module FUSS_Mod_Explicit
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Explicit_Step

contains

  subroutine Explicit_Step ( domain, External_Function )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m
    use FUSS_Global_m
    use FUSS_Mod_dt,         only: Set_Global_dt, Compute_dt
    use FUSS_Mod_Metrics,    only: Update_Mesh
    use FUSS_Mod_ALE_Remap,  only: Remap_State_ALE
    use FUSS_Lib_Ghost,      only: Fill_Ghost_Cell
    use FUSS_Mod_Fluxes,     only: Fluxes
    use FUSS_Mod_BC_Fluxes,  only: BC_Fluxes
    use FUSS_Mod_Newstate,   only: RK_Newstate
    use FUSS_Mod_Diagnostic, only: Compute_Residual
    use FUSS_Mod_MPI, only: is_local_block, mpi_reduce_sum_r8, &
                           mpi_is_root, mpi_bcast_logical, mpi_bcast_integer
    use FUSS_Mod_Timers, only: timer_iter_begin, timer_iter_end, timer_report, &
                               timer_sync_begin, timer_sync_end
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    external :: External_Function
    ! Local
    logical  :: endsim, iosim, endmg, advance_time
    integer  :: i_rk, b, level
    real(R8) :: average

    call timer_iter_begin()

    level = obj_multigrid%MG_level
    obj_multigrid%change_MG = .false.
    domain(level) % iter = domain(level) % iter + 1
    obj_sim_param%iter_from_call = obj_sim_param%iter_from_call + 1
    obj_sim_param%iter_general   = obj_sim_param%iter_general + 1

    ! --- Time step, on the geometry this step starts from.
    advance_time = .false.
    if (obj_sim_param%HYDRA_time_accurate) then
      call Set_Global_dt ( domain(level) )
      advance_time = .true.
    else
      domain(level) % dtglobal = 1d5
      call Compute_dt ( domain(level), obj_time_scheme%vnn, obj_time_scheme%rampa_vnn_iter )  ! Compute local and minimum time step
      if ( obj_time_scheme%time_accurate ) then
        call Set_Global_dt ( domain(level) )  ! Time-accurate: apply global minimum time step
        advance_time = .true.
      endif
    endif

    ! --- Move the mesh over this step, then advance the clock.
    !
    ! This must come AFTER dt is known: the motion law is evaluated at t+dt, and
    ! on the very first step dtglobal is not yet set to anything meaningful.
    ! Placing the update before Compute_dt (as first written) fed the law an
    ! uninitialised step size and the mesh jumped to an arbitrary position.
    !
    ! dt is therefore taken from the geometry at the START of the step, a
    ! first-order lag in the stability limit. It costs accuracy, never
    ! conservation: the swept volumes come from the actual old and new node
    ! positions, not from dt, so the discrete GCL holds regardless.
    !
    ! Outside the OpenMP region below -- the geometry kernels open their own.
    call Update_Mesh ( domain(level), domain(level)%time, domain(level)%dtglobal )

    if ( advance_time ) then
      domain(level) % time = domain(level) % time + domain(level) % dtglobal
    endif

    !$omp parallel

    ! Conservative remap onto the moved mesh -- ONCE per step, before the RK
    ! stages. Ghost temperatures must be current for the upwind lookup at block
    ! interfaces, hence the fill. Both are no-ops when the mesh is static.
    if ( obj_mesh_motion%enabled ) then
      call Fill_Ghost_Cell ( domain(level) )
      call Remap_State_ALE ( domain(level) )
    endif

    call Copy_State ( domain(level) )

    rk: do i_rk = 1, obj_time_scheme%n_RK
        
      call Fill_Ghost_Cell ( domain(level) )               ! Fill ghost cells
      call Fluxes ( domain(level) )                        ! Diffusive fluxes 
      call BC_Fluxes ( domain(level) )                     ! Boundary fluxes

      call External_Function ( domain(level) )             ! External function (e.g. source terms)

      call RK_Newstate ( domain(level), i_rk )             ! State update

    enddo rk
    !$omp end parallel

    ! Save residuals
    if ( (obj_sim_param%iter_from_call == 1) .or. (mod (domain(level) % iter, obj_io%res_diter) == 0) ) then
      obj_sim_param%residuotot = 0.d0
      do b = 1, domain(level) % nb
        if (.not. is_local_block(b)) cycle
        call Compute_Residual ( new=domain(level)%blk(b)%T, &
                                old=domain(level)%blk(b)%TO, &
                                dt=domain(level)%blk(b)%dtlocal, &
                                n=domain(level)%blk(b)%dim, &
                                average=average, &
                                total=obj_sim_param%residuotot )
      enddo
      call timer_sync_begin()
      call mpi_reduce_sum_r8(obj_sim_param%residuotot)
      call timer_sync_end()
      if (mpi_is_root) obj_sim_param%residuotot = sqrt ( obj_sim_param%residuotot ) ! L2 norm time derivative
    endif
    
    ! Determine simulation control flags on root, then broadcast to all ranks
    !
    ! Both flags MUST be initialised: only one of them is assigned below, and
    ! BOTH are read. On a single-grid run `endmg` was never assigned and was
    ! then read in the `elseif` -- so whatever was on the stack decided whether
    ! the step counted as an output step. That is not hypothetical: running the
    ! same case twice with the same binary produced solution files numbered
    ! 1..101 one time and 10..1010 the next, because id_stampa advances once per
    ! TODO==2 event. The solution itself was bit-identical both times, which is
    ! exactly what made it hard to see -- only the file NAMES moved, so it read
    ! as a harness quirk rather than as uninitialised state.
    !
    ! The multigrid direction is worse and has no symptom at all: for level /= 1
    ! it is `endsim` that goes unassigned, and a stray .true. there sets TODO=3
    ! and ends the run early, reported as a normal completion.
    endsim = .false.
    endmg  = .false.

    if (mpi_is_root) then
      if (level == 1) then
        endsim = ( obj_sim_param%iter_from_call >= domain(1) % itermax ) &
            .or. ( obj_sim_param%residuotot <= obj_sim_param%res_threshold ) &
            .or. ( domain(1) % time >= obj_sim_param%time_threshold )
      else
        endmg = ( obj_sim_param%iter_from_call >= domain(level) % itermax ) &
            .or. ( obj_sim_param%residuotot <= obj_sim_param%res_threshold ) &
            .or. ( domain(level) % time >= obj_sim_param%time_threshold )
      endif

      iosim  = ( mod (domain(level) % iter, obj_io%sol_diter) == 0) &
            .or. ( domain(level) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime )

      if ( endsim ) then
        obj_sim_param%TODO = 3
      elseif ( iosim .or. endmg ) then
        obj_sim_param%TODO = 2
        if (endmg .and. obj_multigrid%MGL > 1) then
          obj_multigrid%change_MG = .true.
        endif
      else
        obj_sim_param%TODO = 1
      endif
    end if

    ! Broadcast simulation control from root to all ranks
    call timer_sync_begin()
    call mpi_bcast_integer(obj_sim_param%TODO)
    call mpi_bcast_logical(obj_multigrid%change_MG)
    call timer_sync_end()

    ! Wall-clock report (collective: every rank runs the same iteration counter)
    call timer_iter_end()
    if ( obj_io%timer_diter > 0 ) then
      if ( mod (domain(level) % iter, obj_io%timer_diter) == 0 ) &
        call timer_report ( level, domain(level) % iter )
    endif

  end subroutine Explicit_Step


  subroutine Copy_State ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, j, k, b

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        
        domain % blk(b) % TO (i,j,k) = domain % blk(b) % T (i,j,k)
      
      enddo; enddo; enddo
      !$omp end do
    enddo

  end subroutine Copy_State

end module FUSS_Mod_Explicit