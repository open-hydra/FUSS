module FUSS_Lib_Newstate
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Newstate_Conservative, Newstate_Primitive

contains

  ! ===========================================================================
  !  CONSERVATIVE NEWSTATE
  ! ===========================================================================
  subroutine Newstate_Conservative ( domain, irk )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_time_scheme, obj_irs, obj_mesh_motion
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

    else if ( .not. obj_mesh_motion%enabled ) then

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

    else

      ! Moving mesh (ALE). Kept as a separate branch, rather than folded into
      ! the expression above, so that the static path's arithmetic is untouched
      ! and reproduces earlier results BIT-identically.
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Compute_Cons_ALE( domain%blk(b), irk, n_rk, b, i, j, k )
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

      !> Arbitrary-Lagrangian-Eulerian update for one cell on a moving mesh.
      !>
      !> Integrating the volumetric enthalpy h [J/m^3] over a control volume
      !> that changes size gives
      !>
      !>   (V h)^{n+1} = (V h)^n + dt*( -R + qvol*V )  +  sum_f h_f dV_swept(f)
      !>
      !> where the last term is the ALE convective flux: enthalpy carried across
      !> the faces as they sweep. Dividing by V^{n+1} and subtracting h gives the
      !> increment used below.
      !>
      !> TWO PROPERTIES THIS FORM HAS, BOTH REQUIRED:
      !>
      !>  1. Free-stream preservation. For a uniform field the geometric term
      !>     h*(V^n - V^{n+1}) and the swept term h*sum_f dV_swept cancel, because
      !>     the discrete GCL guarantees sum_f dV_swept == V^{n+1} - V^n. A uniform
      !>     temperature therefore stays uniform under ARBITRARY mesh motion, to
      !>     round-off. This is what test 5.1 checks, and it is the reason the
      !>     swept volumes must come from MORPH rather than be re-derived here.
      !>
      !>  2. Conservation. The face value h_f is upwinded on the sign of
      !>     dV_swept: a face that sweeps outward (cell gains volume) carries the
      !>     NEIGHBOUR's enthalpy, one that sweeps inward carries this cell's.
      !>     Across a shared face the two cells see equal and opposite dV_swept,
      !>     so they both select the SAME upwind value and their contributions
      !>     cancel exactly in the global sum. Using a centred average here would
      !>     preserve free-stream just as well but would NOT telescope, and the
      !>     conservation audit would fail.
      !>
      !> Reduces exactly to the static formula when V^n == V^{n+1} and
      !> dV_swept == 0.
      subroutine Compute_Cons_ALE ( blk, irk, n_rk, b, i, j, k )
        use FUSS_Advanced_Types_m, only: FUSS_block_type
        implicit none
        type(FUSS_block_type), intent(inout) :: blk
        integer,  intent(in) :: irk, n_rk, b, i, j, k
        ! Local
        integer  :: f
        real(R8) :: hO, h, h_nb, h_up, ale, vnew, vold, dV, res

        vnew = blk%vol(i,j,k)
        vold = blk%vol_old(i,j,k)

        call co_H( blk%matID(i,j,k), blk%TO(i,j,k), hO )
        call co_H( blk%matID(i,j,k), blk%T(i,j,k),  h  )

        ! ALE convective flux, upwinded on the swept-volume sign.
        ale = 0.0d0
        do f = 1, 6
          dV = blk%dV_swept(f,i,j,k)
          if ( dV == 0.0d0 ) cycle
          if ( dV > 0.0d0 ) then
            call neighbour_h ( blk, f, i, j, k, h_nb )
            h_up = h_nb
          else
            h_up = h
          endif
          ale = ale + h_up * dV
        enddo

        ! Increment in h. Written as one quotient so the two geometric terms
        ! cancel as completely as floating point allows for a uniform field.
        res = ( h * ( vold - vnew )                                        &
                + blk%dtlocal(i,j,k) * ( -blk%R(i,j,k) + blk%qvol(i,j,k)*vnew ) &
                + ale ) / vnew

        h = RK_stage ( irk, n_rk, h, hO, res )

        call co_T( blk%matID(i,j,k), h, blk%T(i,j,k) )
        call Check_And_Fix_State( blk%T(i,j,k), b, i, j, k )

      end subroutine Compute_Cons_ALE


      !> Volumetric enthalpy of the cell across face f. Ghost cells carry valid
      !> T and matID (filled every RK stage by Fill_Ghost_Cell / fill_matIDg),
      !> so boundary faces need no special case.
      subroutine neighbour_h ( blk, f, i, j, k, h_nb )
        use FUSS_Advanced_Types_m, only: FUSS_block_type
        implicit none
        type(FUSS_block_type), intent(in)  :: blk
        integer,               intent(in)  :: f, i, j, k
        real(R8),              intent(out) :: h_nb
        ! Local
        integer :: ii, jj, kk

        ii = i; jj = j; kk = k
        select case (f)
          case (1); ii = i - 1
          case (2); ii = i + 1
          case (3); jj = j - 1
          case (4); jj = j + 1
          case (5); kk = k - 1
          case (6); kk = k + 1
        end select

        call co_H( blk%matID(ii,jj,kk), blk%T(ii,jj,kk), h_nb )

      end subroutine neighbour_h


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
  !  PRIMITIVE NEWSTATE
  ! ===========================================================================
  subroutine Newstate_Primitive ( domain, irk )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_time_scheme, obj_irs, obj_mesh_motion
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
      stop "NaN or T<0 detected"
    endif

  end subroutine Check_And_Fix_State


end module FUSS_Lib_Newstate
