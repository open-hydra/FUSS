module FUSS_Lib_Newstate
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Newstate_Conservative, Newstate_Primitive, Newstate_Conservative_Pyrolysis

contains

  ! ===========================================================================
  !  CONSERVATIVE NEWSTATE
  ! ===========================================================================
  subroutine Newstate_Conservative ( domain, irk )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_time_scheme, obj_irs
    use FUSS_Global_m
    use FUSS_Lib_Solid
    use FUSS_Lib_RK
    use FUSS_Lib_IRS
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain 
    integer, intent(in)                   :: irk
    ! Local
    integer :: i, j, k, b
    integer :: n_rk
    logical :: irs_enabled
    real(R8) :: irs_beta

    n_rk        = obj_time_scheme%n_RK
    irs_enabled = obj_irs%enabled
    irs_beta    = obj_irs%beta

    ! ------------------------------------------------------------------
    ! Branch on IRS once outside all loops – avoids repeated evaluation
    ! ------------------------------------------------------------------
    if ( irs_enabled ) then

      ! PHASE 1: compute residuals only (no state update yet)
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Scale_Residual_Cons( domain%blk(b)%R(i,j,k),       &
                                    domain%blk(b)%dtlocal(i,j,k), &
                                    domain%blk(b)%vol(i,j,k),     &
                                    domain%blk(b)%qvol(i,j,k)      )
        enddo; enddo; enddo
        !$omp end do
      enddo

      ! PHASE 2: smooth residuals (serial across blocks by design)
      call Residual_Smoothing( domain, irs_beta )

      ! PHASE 3: update state after smoothing
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Update_State_Cons_IRS( domain%blk(b)%T(i,j,k),     &
                                      domain%blk(b)%TO(i,j,k),    &
                                      domain%blk(b)%matID(i,j,k), &
                                      domain%blk(b)%R(i,j,k),     &
                                      irk, n_rk, b, i, j, k       )
        enddo; enddo; enddo
        !$omp end do
      enddo

    else

      ! No IRS: one pass – compute and update in same loop
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Compute_Cons_NoIRS( domain%blk(b)%T(i,j,k),     &
                                   domain%blk(b)%TO(i,j,k),    &
                                   domain%blk(b)%matID(i,j,k),   &
                                   domain%blk(b)%R(i,j,k),     &
                                   domain%blk(b)%dtlocal(i,j,k), &
                                   domain%blk(b)%vol(i,j,k),     &
                                   domain%blk(b)%qvol(i,j,k),    &
                                   irk, n_rk, b, i, j, k         )
        enddo; enddo; enddo
        !$omp end do
      enddo

    endif

    contains

      ! Scale residual in-place for conservative integration (IRS path)
      subroutine Scale_Residual_Cons( residual, dt, volume, qvol )
        implicit none
        real(R8), intent(inout) :: residual
        real(R8), intent(in)    :: dt, volume, qvol
        residual = ( - residual/volume + qvol) * dt
      end subroutine Scale_Residual_Cons

      ! Apply RK stage + cons2prim after IRS smoothing
      subroutine Update_State_Cons_IRS( T, TO, matID, residual, irk, n_rk, b, i, j, k )
        implicit none
        real(R8), intent(inout) :: T, residual
        real(R8), intent(in)    :: TO, matID
        integer,  intent(in)    :: irk, n_rk, b, i, j, k
        ! Local
        real(R8) :: hO, h

        call co_H( matID, TO, hO )
        call co_H( matID, T, h )

        h = RK_stage ( irk, n_rk, h, hO, residual )

        call co_T( matID, h, T )

        call Check_And_Fix_State( T, b, i, j, k )

      end subroutine Update_State_Cons_IRS

      ! Full compute+update without IRS
      subroutine Compute_Cons_NoIRS ( T, TO, matID, residual, dt, volume, qvol, irk, n_rk, b, i, j, k )
        implicit none
        real(R8), intent(inout) :: T, residual
        real(R8), intent(in)    :: TO, matID, dt, volume, qvol
        integer,  intent(in)    :: irk, n_rk, b, i, j, k
        ! Local
        real(R8) :: hO, h

        residual = ( - residual/volume + qvol) * dt

        call co_H( matID, TO, hO )
        call co_H( matID, T, h )

        h = RK_stage ( irk, n_rk, h, hO, residual )

        call co_T( matID, h, T )

        call Check_And_Fix_State( T, b, i, j, k )
      
      end subroutine Compute_Cons_NoIRS

  end subroutine Newstate_Conservative


  ! ===========================================================================
  !  CONSERVATIVE NEWSTATE WITH PYROLYSIS (plan 06 section 3.4)
  !
  !  Integrates the stored energy E(T, chi) and the reaction advancements chi
  !  with the same RK_stage weights, from the start-of-step state (TO, chiO)
  !  and the stage increments (R folded as the no-IRS residual, dchi from
  !  Pyrolysis_Stage), inverts E -> T per cell and refreshes the property
  !  caches. Pyrolysis_Stage already added the gas-advection flux to R. The
  !  two legacy routines above and below are untouched; IRS is refused.
  ! ===========================================================================
  subroutine Newstate_Conservative_Pyrolysis ( domain, irk )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_time_scheme, obj_pyrolysis, obj_darcy
    use FUSS_Global_m
    use FUSS_Lib_Pyrolysis,  only: Energy_Of, Temperature_Of, Solid_State
    use FUSS_Lib_RK
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain 
    integer, intent(in)                   :: irk
    ! Local
    integer  :: i, j, k, b, m, r, n_rk
    real(R8) :: residual, E_n, E_k, E_new, T_new, tau, eps_pg, p_gas

    ! Type 2 (plan 07): its own routine, so that this one stays plan 06's text
    if ( obj_darcy%enabled ) then
      call Newstate_Conservative_Darcy( domain, irk )
      return
    endif

    n_rk  = obj_time_scheme%n_RK
    p_gas = obj_pyrolysis%p_gas

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(3) schedule(static) private(i,j,k,m,r,residual,E_n,E_k,E_new,T_new,tau,eps_pg)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        m = nint( domain%blk(b)%matID(i,j,k) )
        residual = ( - domain%blk(b)%R(i,j,k) / domain%blk(b)%vol(i,j,k) + domain%blk(b)%qvol(i,j,k) ) &
                 * domain%blk(b)%dtlocal(i,j,k)
        E_n   = Energy_Of( m, domain%blk(b)%TO(i,j,k), domain%blk(b)%chiO(:,i,j,k), p_gas )
        E_k   = Energy_Of( m, domain%blk(b)%T(i,j,k),  domain%blk(b)%chi(:,i,j,k),  p_gas )
        E_new = RK_stage( irk, n_rk, E_k, E_n, residual )
        do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
          domain%blk(b)%chi(r,i,j,k) = RK_stage( irk, n_rk, domain%blk(b)%chi(r,i,j,k), &
                                                 domain%blk(b)%chiO(r,i,j,k), domain%blk(b)%dchi(r,i,j,k) )
        enddo
        call Temperature_Of( m, E_new, domain%blk(b)%chi(:,i,j,k), p_gas, domain%blk(b)%T(i,j,k), T_new, b, i, j, k )
        domain%blk(b)%T(i,j,k) = T_new
        call Check_And_Fix_State( T_new, b, i, j, k )
        call Solid_State( m, T_new, domain%blk(b)%chi(:,i,j,k), domain%blk(b)%rho_s(i,j,k), &
                          domain%blk(b)%cp_s(i,j,k), domain%blk(b)%k_s(i,j,k), tau, eps_pg )
      enddo; enddo; enddo
      !$omp end do
    enddo

  end subroutine Newstate_Conservative_Pyrolysis


  ! ===========================================================================
  !  CONSERVATIVE NEWSTATE WITH THE TYPE-2 GAS CLOSURE (plan 07 section 3.5)
  !
  !  E = H_s(T, tau) + G h_g(p, T) - eps_pg p with the stored gas mass G a state
  !  (FUSS_Lib_Pyrolysis, type-2 section). The pressure solve of the step fixed
  !  G^(n+1) = eta_st p^(n+1) and the outflow that carries the gas out; at every
  !  RK stage the state the energy is inverted at carries the share of that
  !  step the stage has accounted for: with c_s the stage's time fraction
  !  (SSP-RK3: 1, 1/2, 1 -- the RK_stage recursion applied to a unit
  !  increment), G_s = G^n + c_s (G^(n+1) - G^n) and p_s likewise. Stage 1
  !  starts from E^n itself (re-evaluating (T^n, chi^n) at another (p, G)
  !  would add stored energy without a flux and break the budget identity).
  ! ===========================================================================
  subroutine Newstate_Conservative_Darcy ( domain, irk )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_time_scheme, obj_pyrolysis
    use FUSS_Global_m
    use FUSS_Lib_Pyrolysis,  only: Energy_Of2, Temperature_Of2, Solid_State
    use FUSS_Lib_RK
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in)                   :: irk
    ! Local
    integer  :: i, j, k, b, m, r, n_rk, s
    real(R8) :: residual, E_n, E_k, E_new, T_new, tau, eps_pg
    real(R8) :: c_prev, c_cur, p_n, p_1, G_n, G_1

    n_rk = obj_time_scheme%n_RK
    ! stage time fractions: c_cur of this stage's output, c_prev of its input
    c_prev = 0.0d0
    c_cur  = 0.0d0
    do s = 1, irk
      c_prev = c_cur
      c_cur  = RK_stage( s, n_rk, c_cur, 0.0d0, 1.0d0 )
    enddo

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(3) schedule(static) private(i,j,k,m,r,residual,E_n,E_k,E_new,T_new,tau,eps_pg,p_n,p_1,G_n,G_1)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        m = nint( domain%blk(b)%matID(i,j,k) )
        residual = ( - domain%blk(b)%R(i,j,k) / domain%blk(b)%vol(i,j,k) + domain%blk(b)%qvol(i,j,k) ) &
                 * domain%blk(b)%dtlocal(i,j,k)
        p_n = domain%blk(b)%p_pgO(i,j,k)
        p_1 = domain%blk(b)%p_pg(i,j,k)
        G_n = domain%blk(b)%eta_stO(i,j,k) * p_n
        G_1 = domain%blk(b)%eta_st(i,j,k)  * p_1
        E_n = Energy_Of2( m, domain%blk(b)%TO(i,j,k), domain%blk(b)%chiO(:,i,j,k), p_n, G_n )
        if ( irk == 1 ) then
          E_k = E_n
        else
          E_k = Energy_Of2( m, domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), &
                            p_n + c_prev * ( p_1 - p_n ), G_n + c_prev * ( G_1 - G_n ) )
        endif
        E_new = RK_stage( irk, n_rk, E_k, E_n, residual )
        do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
          domain%blk(b)%chi(r,i,j,k) = RK_stage( irk, n_rk, domain%blk(b)%chi(r,i,j,k), &
                                                 domain%blk(b)%chiO(r,i,j,k), domain%blk(b)%dchi(r,i,j,k) )
        enddo
        call Temperature_Of2( m, E_new, domain%blk(b)%chi(:,i,j,k), p_n + c_cur * ( p_1 - p_n ), &
                              G_n + c_cur * ( G_1 - G_n ), domain%blk(b)%T(i,j,k), T_new, b, i, j, k )
        domain%blk(b)%T(i,j,k) = T_new
        call Check_And_Fix_State( T_new, b, i, j, k )
        call Solid_State( m, T_new, domain%blk(b)%chi(:,i,j,k), domain%blk(b)%rho_s(i,j,k), &
                          domain%blk(b)%cp_s(i,j,k), domain%blk(b)%k_s(i,j,k), tau, eps_pg )
      enddo; enddo; enddo
      !$omp end do
    enddo

  end subroutine Newstate_Conservative_Darcy


  ! ===========================================================================
  !  PRIMITIVE NEWSTATE
  ! ===========================================================================
  subroutine Newstate_Primitive ( domain, irk )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_time_scheme, obj_irs
    use FUSS_Global_m
    use FUSS_Lib_Solid
    use FUSS_Lib_RK
    use FUSS_Lib_IRS
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain 
    integer, intent(in)                   :: irk
    ! Local
    integer :: i, j, k, b
    integer :: n_rk
    logical :: irs_enabled
    real(R8) :: irs_beta

    n_rk        = obj_time_scheme%n_RK
    irs_enabled = obj_irs%enabled
    irs_beta    = obj_irs%beta

    ! ------------------------------------------------------------------
    ! Branch on IRS once outside all loops
    ! ------------------------------------------------------------------
    if ( irs_enabled ) then

      ! PHASE 1: convert residual to primitive space only
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Compute_Prim_Residual( domain%blk(b)%T(i,j,k),       &
                                      domain%blk(b)%matID(i,j,k),   &
                                      domain%blk(b)%qvol(i,j,k),    &
                                      domain%blk(b)%R(i,j,k),       &
                                      domain%blk(b)%dtlocal(i,j,k), &
                                      domain%blk(b)%vol(i,j,k)      )
        enddo; enddo; enddo
        !$omp end do
      enddo

      call Residual_Smoothing( domain, irs_beta )

      ! PHASE 3: apply RK stage after smoothing
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          domain%blk(b)%T(i,j,k) = RK_stage( irk, n_rk,               &
                                             domain%blk(b)%T(i,j,k),  &
                                             domain%blk(b)%TO(i,j,k), &
                                             domain%blk(b)%R(i,j,k)   )
          call Check_And_Fix_State( domain%blk(b)%T(i,j,k), b, i, j, k )
        enddo; enddo; enddo
        !$omp end do
      enddo
    
    else

      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Compute_And_Update_Prim( domain%blk(b)%T(i,j,k),       &
                                        domain%blk(b)%TO(i,j,k),      &
                                        domain%blk(b)%matID(i,j,k),   &
                                        domain%blk(b)%qvol(i,j,k),    &
                                        domain%blk(b)%R(i,j,k),       &
                                        domain%blk(b)%dtlocal(i,j,k), &
                                        domain%blk(b)%vol(i,j,k),     &
                                        irk, n_rk, b, i, j, k         )
        enddo; enddo; enddo
        !$omp end do
      enddo

    endif

    contains

      ! Converts conservative residual to primitive residual in-place
      ! (no RK update – used when IRS smoothing follows)
      subroutine Compute_Prim_Residual ( T, matID, qvol, residual, dt, volume )
        implicit none
        real(R8), intent(inout) :: T, residual
        real(R8), intent(in)    :: matID, dt, qvol, volume
        ! Local
        real(R8) :: rho, cp, drhodT, dcpdT
        
        call co_Rho(matID, T, rho)
        call co_Cp(matID, T, cp)
        call co_DrhoDT(matID, T, drhodT)
        call co_DcpDT(matID, T, dcpdT)

        ! Residuo sulle primitive (proprietà variabili)
        residual = ( - residual/volume + qvol) *dt / (T*cp*drhodT + rho*cp + T*rho*dcpdT)

      end subroutine Compute_Prim_Residual

      ! Full primitive residual + RK update in one shot (no IRS)
      subroutine Compute_And_Update_Prim ( T, TO, matID, qvol, residual, dt, volume, irk, n_rk, b, i, j, k )

        implicit none
        real(R8), intent(inout) :: T, residual
        real(R8), intent(in)    :: TO, matID, dt, qvol, volume
        integer,  intent(in)    :: irk, n_rk, b, i, j, k
        
        call Compute_Prim_Residual( T, matID, qvol, residual, dt, volume )
        T = RK_Stage ( irk, n_rk, T, TO, residual )
        call Check_And_Fix_State( T, b, i, j, k )

      end subroutine Compute_And_Update_Prim

  end subroutine Newstate_Primitive


  ! ===========================================================================
  !  STATE SANITY CHECK AND FLOOR ENFORCEMENT
  ! ===========================================================================
  subroutine Check_And_Fix_State ( T, b, i, j, k )
    implicit none
    real(R8), intent(in) :: T
    integer, intent(in)  :: b, i, j, k

    ! Bail out on NaN or negative temperature
    if ( isnan(T) .or. T < 0d0 ) then
      write(*,'(A,4I4)') "Integration failed at b, i, j, k:", b, i, j, k
      write(*,*) T
      write(*,'(A)') " [ERROR] NaN or T<0 detected"
      error stop 3
    endif

  end subroutine Check_And_Fix_State


end module FUSS_Lib_Newstate
