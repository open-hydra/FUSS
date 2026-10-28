module FUSS_Lib_Newstate
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Newstate_Conservative, Newstate_Primitive

  !> FAS coarse-level updates halved because the forcing pushed T outside the table (no-IRS
  !> conservative path, which otherwise stops); reported and reset by FUSS_postprocess every
  !> shell-diter iterations.
  integer, public :: fas_update_fallbacks = 0

  !> The FAS guards (the prolongated correction in Lib_Multigrid, the coarse smoothing update here)
  !> halve down to 2**-10, then do not apply. The IRS admissibility fallback that used the same
  !> constant was removed in October 2026 (decision D2 of the hydra sync, as MOSE): a smoothed
  !> update that leaves the properties table stops the run like an explicit one.
  integer, parameter, public :: FAS_HALVINGS = 10
  public :: Admissible_T

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
                                   irk, n_rk, b, i, j, k,        &
                                   domain%fas_forcing            )
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

      ! RK stage after IRS smoothing. A smoothed update that leaves the properties table stops
      ! the run, as an explicit one does (Check_And_Fix_State): the halving fallback of 2026-10-14
      ! was removed in October 2026 (decision D2 of the hydra sync, as MOSE a24a236).
      subroutine Update_State_Cons_IRS( T, TO, matID, residual, irk, n_rk, b, i, j, k )
        implicit none
        real(R8), intent(inout) :: T, residual
        real(R8), intent(in)    :: TO, matID
        integer,  intent(in)    :: irk, n_rk, b, i, j, k
        ! Local
        real(R8) :: hO, h

        call co_H( matID, TO, hO )
        call co_H( matID, T, h )

        call co_T( matID, RK_stage ( irk, n_rk, h, hO, residual ), T )
        call Check_And_Fix_State( T, b, i, j, k )

      end subroutine Update_State_Cons_IRS

      ! Full compute+update without IRS. On a FAS coarse level (guarded = the level carries a
      ! forcing) an update that leaves the properties table is halved until admissible (at most
      ! FAS_HALVINGS times, then not applied) instead of stopping the run: a forced coarse
      ! problem is not a physical time evolution.
      subroutine Compute_Cons_NoIRS ( T, TO, matID, residual, dt, volume, qvol, irk, n_rk, b, i, j, k, guarded )
        implicit none
        real(R8), intent(inout) :: T, residual
        real(R8), intent(in)    :: TO, matID, dt, volume, qvol
        integer,  intent(in)    :: irk, n_rk, b, i, j, k
        logical,  intent(in)    :: guarded
        ! Local
        real(R8) :: hO, h, h0stage, theta
        integer  :: n

        residual = ( - residual/volume + qvol) * dt

        call co_H( matID, TO, hO )
        call co_H( matID, T, h )
        h0stage = h

        h = RK_stage ( irk, n_rk, h, hO, residual )

        call co_T( matID, h, T )

        if ( guarded .and. .not. Admissible_T( T ) ) then
          !$omp atomic
          fas_update_fallbacks = fas_update_fallbacks + 1
          theta = 0.5d0
          do n = 1, FAS_HALVINGS
            call co_T( matID, RK_stage ( irk, n_rk, h0stage, hO, theta * residual ), T )
            if ( Admissible_T( T ) ) exit
            theta = 0.5d0 * theta
          enddo
          ! theta = 0: a convex combination of the admissible stage and time-level states
          if ( .not. Admissible_T( T ) ) call co_T( matID, RK_stage ( irk, n_rk, h0stage, hO, 0d0 ), T )
        endif

        call Check_And_Fix_State( T, b, i, j, k )

      end subroutine Compute_Cons_NoIRS

  end subroutine Newstate_Conservative


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
          call Update_State_Prim_IRS( domain%blk(b)%T(i,j,k),   &
                                      domain%blk(b)%TO(i,j,k),  &
                                      domain%blk(b)%R(i,j,k),   &
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

        ! dT = dh / (dh/dT), h = int rho cp dT (the Energy column): dh/dT = rho cp
        residual = ( - residual/volume + qvol) * dt / ( rho * cp )

      end subroutine Compute_Prim_Residual

      ! Apply the RK stage after IRS smoothing; an update outside the table stops the run, as
      ! in Update_State_Cons_IRS
      subroutine Update_State_Prim_IRS ( T, TO, residual, irk, n_rk, b, i, j, k )
        implicit none
        real(R8), intent(inout) :: T
        real(R8), intent(in)    :: TO, residual
        integer,  intent(in)    :: irk, n_rk, b, i, j, k

        T = RK_stage ( irk, n_rk, T, TO, residual )
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
  ! T the property lookups can evaluate (they read the entries int(T) and int(T)+1): Tmin <= T < Tmax
  logical function Admissible_T ( T )
    use FUSS_Config_Types_m, only: obj_table
    use FUSS_Finite_m, only: Not_Finite
    implicit none
    real(R8), intent(in) :: T

    Admissible_T = .false.
    if ( Not_Finite(T) ) return
    if ( T < real(obj_table%Tmin, R8) .or. T >= real(obj_table%Tmax, R8) ) return
    Admissible_T = .true.

  end function Admissible_T


  subroutine Check_And_Fix_State ( T, b, i, j, k )
    use FUSS_Config_Types_m, only: obj_table
    use FUSS_Mod_MPI, only: stop_from_thread
    implicit none
    real(R8), intent(in) :: T
    integer, intent(in)  :: b, i, j, k

    ! stop_from_thread (exit code 1), not mpi_stop_all or mpi_abort_all: this runs in OpenMP threads, MPI is FUNNELED
    if ( .not. Admissible_T( T ) ) then
      write(*,'(A,4I5)') ' [ERROR] Integration failed at b, i, j, k:', b, i, j, k
      write(*,'(A,ES14.6,A,I0,A,I0,A)') ' [ERROR] T = ', T, ' K, outside the properties table [', &
        obj_table%Tmin, ', ', obj_table%Tmax, ') K'
      call stop_from_thread('NaN or T outside the properties table')
    endif

  end subroutine Check_And_Fix_State


end module FUSS_Lib_Newstate
