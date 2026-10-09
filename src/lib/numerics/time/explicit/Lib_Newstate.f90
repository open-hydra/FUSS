module FUSS_Lib_Newstate
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Newstate_Conservative, Newstate_Primitive

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
    ! past its ends. error stop, not mpi_abort_all: this runs inside OpenMP threads (MPI is
    ! FUNNELED); unlike stop, it also ends the run with a non-zero exit status.
    if ( .not. Admissible_T( T ) ) then
      write(*,'(A,4I5)') ' [ERROR] Integration failed at b, i, j, k:', b, i, j, k
      write(*,'(A,ES14.6,A,I0,A,I0,A)') ' [ERROR] T = ', T, ' K, outside the properties table [', &
        obj_table%Tmin, ', ', obj_table%Tmax, ') K'
      error stop 'NaN or T outside the properties table'
    endif

  end subroutine Check_And_Fix_State


end module FUSS_Lib_Newstate
