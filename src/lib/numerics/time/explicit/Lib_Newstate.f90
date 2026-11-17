module FUSS_Lib_Newstate
  use FUSS_Mod_MPI, only: fuss_abort_in_thread
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Newstate_Conservative, Newstate_Primitive, Newstate_Conservative_Pyrolysis

  !> IRS fallbacks of this rank (cell updates that took their own residual), reported and reset by
  !> FUSS_postprocess every shell-diter iterations.
  integer, public :: irs_fallbacks = 0

  integer, parameter :: IRS_HALVINGS = 10   ! own update halved down to 2**-10, then not applied

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
          domain%blk(b)%RS0(i,j,k) = domain%blk(b)%R(i,j,k)      ! the cell's own update, kept
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
                                      domain%blk(b)%RS0(i,j,k),   &
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

      ! Apply RK stage + cons2prim after IRS smoothing.
      ! The smoothed update of a cell mixes its neighbours' updates: above the explicit limit it can
      ! overshoot, beyond the temperatures of the boundaries and out of the properties table. If the
      ! new state is not admissible (Admissible_T), the cell takes its own update (residual_own, the
      ! scaled residual before smoothing), halved until it is; with no update at all the stage is a
      ! convex combination of admissible states, which is admissible. An admissible smoothed update
      ! is used as before.
      subroutine Update_State_Cons_IRS( T, TO, matID, residual, residual_own, irk, n_rk, b, i, j, k )
        implicit none
        real(R8), intent(inout) :: T, residual
        real(R8), intent(in)    :: TO, matID, residual_own
        integer,  intent(in)    :: irk, n_rk, b, i, j, k
        ! Local
        real(R8) :: hO, h, Tnew, theta
        integer  :: n

        call co_H( matID, TO, hO )
        call co_H( matID, T, h )

        call co_T( matID, RK_stage ( irk, n_rk, h, hO, residual ), Tnew )

        if ( .not. Admissible_T( Tnew ) ) then
          !$omp atomic
          irs_fallbacks = irs_fallbacks + 1
          theta = 1d0
          do n = 0, IRS_HALVINGS
            call co_T( matID, RK_stage ( irk, n_rk, h, hO, theta * residual_own ), Tnew )
            if ( Admissible_T( Tnew ) ) exit
            theta = 0.5d0 * theta
          enddo
          if ( .not. Admissible_T( Tnew ) ) call co_T( matID, RK_stage ( irk, n_rk, h, hO, 0d0 ), Tnew )
        endif

        T = Tnew
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
          domain%blk(b)%RS0(i,j,k) = domain%blk(b)%R(i,j,k)      ! the cell's own update, kept
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
          call Update_State_Prim_IRS( domain%blk(b)%T(i,j,k),   &
                                      domain%blk(b)%TO(i,j,k),  &
                                      domain%blk(b)%R(i,j,k),   &
                                      domain%blk(b)%RS0(i,j,k), &
                                      irk, n_rk, b, i, j, k     )
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
        real(R8) :: rho, cp

        call co_Rho(matID, T, rho)
        call co_Cp(matID, T, cp)

        ! dT = dh / (dh/dT), with h(T) = int rho cp dT the volumetric enthalpy of the conservative path
        ! (the Energy column of properties.dat): dh/dT = rho cp, also with temperature-dependent properties.
        ! Dividing by d(rho cp T)/dT instead, as before 2026-10-09, added T d(rho cp)/dT to the heat capacity.
        residual = ( - residual/volume + qvol) * dt / ( rho * cp )

      end subroutine Compute_Prim_Residual

      ! Apply the RK stage after IRS smoothing, with the fallback of Update_State_Cons_IRS
      subroutine Update_State_Prim_IRS ( T, TO, residual, residual_own, irk, n_rk, b, i, j, k )
        implicit none
        real(R8), intent(inout) :: T
        real(R8), intent(in)    :: TO, residual, residual_own
        integer,  intent(in)    :: irk, n_rk, b, i, j, k
        ! Local
        real(R8) :: Tnew, theta
        integer  :: n

        Tnew = RK_stage ( irk, n_rk, T, TO, residual )

        if ( .not. Admissible_T( Tnew ) ) then
          !$omp atomic
          irs_fallbacks = irs_fallbacks + 1
          theta = 1d0
          do n = 0, IRS_HALVINGS
            Tnew = RK_stage ( irk, n_rk, T, TO, theta * residual_own )
            if ( Admissible_T( Tnew ) ) exit
            theta = 0.5d0 * theta
          enddo
          if ( .not. Admissible_T( Tnew ) ) Tnew = RK_stage ( irk, n_rk, T, TO, 0d0 )
        endif

        T = Tnew
        call Check_And_Fix_State( T, b, i, j, k )

      end subroutine Update_State_Prim_IRS

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
  !  STATE SANITY CHECK
  ! ===========================================================================
  ! A temperature the properties tables can evaluate: co_Kappa, co_Rho, co_Cp and co_H interpolate
  ! between the entries at int(T) and int(T)+1, so T must lie in [Tmin, Tmax) of the tables.
  logical function Admissible_T ( T )
    use FUSS_Config_Types_m, only: obj_table
    implicit none
    real(R8), intent(in) :: T

    Admissible_T = .false.
    if ( isnan(T) ) return
    if ( T < real(obj_table%Tmin, R8) .or. T >= real(obj_table%Tmax, R8) ) return
    Admissible_T = .true.

  end function Admissible_T


  subroutine Check_And_Fix_State ( T, b, i, j, k )
    use FUSS_Config_Types_m, only: obj_table
    implicit none
    real(R8), intent(in) :: T
    integer, intent(in)  :: b, i, j, k

    ! Stop on NaN or on a temperature outside the properties table, where the lookups would read
    ! past its ends (upstream f3401c5; the check was NaN or T < 0 only). With MPI (plan 11): fuss_abort
    ! on the master thread or one rank keeps the code; a worker thread of a multi-rank run stops its
    ! process (the launcher reports 255). Code 3 = numerical, until plan 14 S6 (D3).
    ! No continuation here: giving this routine the cell's previous temperature (an
    ! intent(inout) T and a T_old argument) changed the code generated for the type-2
    ! RK loop and moved 512 manifest rows by round-off (measured 2026-10-17).
    if ( .not. Admissible_T( T ) ) then
      write(*,'(A,4I5)') ' [ERROR] Integration failed at b, i, j, k:', b, i, j, k
      write(*,'(A,ES14.6,A,I0,A,I0,A)') ' [ERROR] T = ', T, ' K, outside the properties table [', &
        obj_table%Tmin, ', ', obj_table%Tmax, ') K'
      call fuss_abort_in_thread ( 3 )
      error stop 3
    endif

  end subroutine Check_And_Fix_State


end module FUSS_Lib_Newstate
