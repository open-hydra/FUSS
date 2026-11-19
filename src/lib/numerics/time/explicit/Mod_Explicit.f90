module FUSS_Mod_Explicit
  use FUSS_Mod_MPI, only: fuss_abort_check
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Explicit_Step, Evaluate_Residual

contains

  subroutine Explicit_Step ( domain, External_Function, inner )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m
    use FUSS_Global_m
    use FUSS_Mod_dt,         only: Set_Global_dt, Compute_dt
    use FUSS_Mod_Metrics,    only: Update_Mesh
    use FUSS_Mod_ALE_Remap,  only: Remap_State_ALE, Remap_Pyrolysis_ALE
    use FUSS_Mod_Recession,  only: Ablation_Update_Env
    use FUSS_Lib_Ghost,      only: Fill_Ghost_Cell, Fill_Ghost_chi
    use FUSS_Mod_Pyrolysis,  only: Pyrolysis_Stage, Pyrolysis_Residual, Pyrolysis_Budget_Rewind
    use FUSS_Mod_Darcy,      only: Solve_Pressure, Darcy_Residual, Darcy_First_Step_Pending, &
                                   Darcy_First_Step_History, Darcy_Restore_State, Darcy_Mass_Correction
    use FUSS_Mod_Fluxes,     only: Fluxes
    use FUSS_Mod_BC_Fluxes,  only: BC_Fluxes
    use FUSS_Mod_Newstate,   only: RK_Newstate
    use FUSS_Mod_Diagnostic, only: Compute_Residual
    use FUSS_IO_BC,          only: has_tdep_bc
    use FUSS_Mod_MPI, only: is_local_block, mpi_reduce_sum_r8_array, &
                           mpi_is_root, mpi_bcast_logical, mpi_bcast_integer
    use FUSS_Mod_Timers, only: timer_iter_begin, timer_iter_end, timer_report, &
                               timer_sync_begin, timer_sync_end
    use FUSS_Hook_Interface_m, only: external_function_if
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    procedure(external_function_if) :: External_Function
    logical, intent(in), optional :: inner
    ! Local
    logical  :: endsim, iosim, endmg, inner_
    integer  :: i_rk, b, level
    real(R8) :: average, t_n, t_np1
    real(R8), allocatable :: res_block(:)   ! residual sum over the cells of each block
    logical  :: advance, converged
    ! Stage times t_n + c dt of euler, RK2 (Heun) and SSP RK3
    real(R8), parameter :: c_stage(3,3) = reshape([0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 1d0, 0.5d0], [3,3], order=[2,1])

    inner_ = .false.
    if (present(inner)) inner_ = inner

    call timer_iter_begin()

    level = obj_multigrid%MG_level
    if (.not. inner_) obj_multigrid%change_MG = .false.
    domain(level) % iter = domain(level) % iter + 1
    if (.not. inner_) then
      obj_sim_param%iter_from_call = obj_sim_param%iter_from_call + 1
      obj_sim_param%iter_general   = obj_sim_param%iter_general + 1
    endif

    ! --- Time step, on the geometry this step starts from. The clock: domain%time runs at the
    ! stage's time inside the RK loop (the 301/302 time tables, Ghost_Wall) and is t^(n+1) after it;
    ! the step's own readers below (mesh motion, env tables, the Darcy log) take t_n / t_np1 explicitly.
    t_n = domain(level) % time
    advance = obj_sim_param%HYDRA_time_accurate .or. obj_time_scheme%time_accurate
    if (obj_sim_param%HYDRA_time_accurate) then
      call Set_Global_dt ( domain(level) )
    else
      domain(level) % dtglobal = 1d5
      call Compute_dt ( domain(level), obj_time_scheme%vnn, obj_time_scheme%rampa_vnn_iter )  ! Compute local and minimum time step
      if ( obj_time_scheme%time_accurate ) call Set_Global_dt ( domain(level) )  ! Time-accurate: global minimum time step
    endif
    t_np1 = t_n
    if ( advance ) t_np1 = t_n + domain(level) % dtglobal

    ! --- Move the mesh over this step.
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
    call Update_Mesh ( domain(level), t_n, domain(level)%dtglobal )

    ! Ablating walls (plan 08): the env values at t^(n+1), before the type-2
    ! pressure solve reads the faces' p_w. Serial, outside the parallel region.
    call Ablation_Update_Env ( domain(level), t_np1 )

    !$omp parallel

    ! T^(n-1) of the Darcy closure, shifted BEFORE the remap (plan 08 section
    ! 2.9 phase 0): the remap moves the histories onto the new mesh, so they
    ! must already hold this step's values. On a static mesh nothing between
    ! here and Copy_State, where the shift used to be, touches TO or T_prev.
    if ( obj_darcy%enabled ) call Darcy_Shift_History ( domain(level) )

    ! Conservative remap onto the moved mesh -- ONCE per step, before the RK
    ! stages. Ghost temperatures must be current for the upwind lookup at block
    ! interfaces, hence the fill (at t_n: the clock has not moved yet). Both are
    ! no-ops when the mesh is static.
    if ( obj_mesh_motion%enabled ) then
      if ( obj_pyrolysis%enabled ) then
        ! the pyrolysis state (plan 08 section 2.9; law = nodeshift): energy,
        ! sub-phase masses and the type-2 gas state, no ghost read
        call Remap_Pyrolysis_ALE ( domain(level) )
      else
        call Fill_Ghost_Cell ( domain(level) )
        call Remap_State_ALE ( domain(level) )
      endif
    endif

    call Copy_State ( domain(level) )

    ! Darcy pressure of this step (plan 07 section 3.3): solved ONCE from the
    ! start-of-step state, before the RK loop; every stage then reads
    ! p^(n+1) as frozen data (PATO solves its mass model once per step and
    ! feeds the energy model explicitly). A new run's first step is run
    ! twice: its mid-step extrapolations (production rate, storage
    ! coefficient) have no history on the first pass; the second pass takes
    ! it from the first pass's end state. The stage loop below is plan 06's,
    ! untouched; the first pass repeats it here, clock included.
    if ( obj_darcy%enabled ) then
      if ( Darcy_First_Step_Pending() ) then
        call Solve_Pressure ( domain(level), t_np1 )
        do i_rk = 1, obj_time_scheme%n_RK
          if (advance) then
            !$omp single
            domain(level) % time = t_n + c_stage(obj_time_scheme%n_RK, i_rk) * domain(level) % dtglobal   ! BC time tables
            !$omp end single
          endif
          call Fill_Ghost_Cell ( domain(level) )
          call Fill_Ghost_chi ( domain(level) )
          call Fluxes ( domain(level) )
          call BC_Fluxes ( domain(level) )
          call Pyrolysis_Stage ( domain(level), i_rk )
          call External_Function ( domain(level) )
          call RK_Newstate ( domain(level), i_rk )
        enddo
        call Darcy_First_Step_History ( domain(level) )
        call Darcy_Restore_State ( domain(level) )
        call Pyrolysis_Budget_Rewind ()
      endif
      call Solve_Pressure ( domain(level), t_np1 )
    endif

    rk: do i_rk = 1, obj_time_scheme%n_RK

      if (advance) then
        !$omp single
        domain(level) % time = t_n + c_stage(obj_time_scheme%n_RK, i_rk) * domain(level) % dtglobal   ! BC time tables
        !$omp end single
      endif
      call Fill_Ghost_Cell ( domain(level) )               ! Fill ghost cells (301/302 time tables at the stage's time)
      if ( obj_pyrolysis%enabled ) call Fill_Ghost_chi ( domain(level) )   ! mirror chi at wall faces (plan 06 section 3.1)
      call Fluxes ( domain(level) )                        ! Diffusive fluxes 
      call BC_Fluxes ( domain(level) )                     ! Boundary fluxes
      if ( obj_pyrolysis%enabled ) call Pyrolysis_Stage ( domain(level), i_rk )   ! kinetics, gas outflow, advection flux, budget

      call External_Function ( domain(level) )             ! External function (e.g. source terms)

      if (domain(level) % fas_forcing) call Add_FAS_Forcing ( domain(level) )   ! FAS coarse level

      call RK_Newstate ( domain(level), i_rk )             ! State update

    enddo rk

    ! Gas-mass closure of the type-2 step (plan 07): the realised minus the
    ! predicted production goes into the stored gas, T re-inverted at constant E
    if ( obj_darcy%enabled ) call Darcy_Mass_Correction ( domain(level) )
    !$omp end parallel
    ! a fatal error that a worker thread recorded in the region (plan 11)
    call fuss_abort_check ()
    if (advance) domain(level) % time = t_n + domain(level) % dtglobal

    if (inner_) then
      call timer_iter_end()
      return
    endif

    ! Save residuals
    if ( (obj_sim_param%iter_from_call == 1) .or. (mod (domain(level) % iter, obj_io%res_diter) == 0) ) then
      ! Sums in a fixed order, the same for any number of threads and ranks: over the cells of each block
      ! (Compute_Residual), then over the blocks in turn on the root; each entry of res_block comes from the one rank
      ! that holds the block, so its sum over the ranks is exact
      obj_sim_param%residuotot = 0.d0
      allocate ( res_block(domain(level)%nb) )
      res_block = 0.d0
      do b = 1, domain(level) % nb
        if (.not. is_local_block(b)) cycle
        call Compute_Residual ( new=domain(level)%blk(b)%T, &
                                old=domain(level)%blk(b)%TO, &
                                dt=domain(level)%blk(b)%dtlocal, &
                                n=domain(level)%blk(b)%dim, &
                                average=average, &
                                total=res_block(b) )
      enddo
      ! With pyrolysis the reaction state is part of the residual (plan 06):
      ! a constant temperature with reactions still advancing is not converged.
      ! Added per block into res_block (plan 14 S5), before the reduce: each entry
      ! still comes from the one rank that holds the block.
      if ( obj_pyrolysis%enabled ) call Pyrolysis_Residual ( domain(level), res_block )
      ! ... and so is the gas pressure under type 2 (plan 07): a stationary
      ! temperature with p still relaxing (the quasi-steady Picard steps of
      ! gates 4.2/4.3) is not converged either.
      if ( obj_darcy%enabled ) call Darcy_Residual ( domain(level), res_block )
      call timer_sync_begin()
      call mpi_reduce_sum_r8_array(res_block, domain(level)%nb)
      call timer_sync_end()
      if (mpi_is_root) then
        do b = 1, domain(level) % nb
          obj_sim_param%residuotot = obj_sim_param%residuotot + res_block(b)
        enddo
        obj_sim_param%residuotot = sqrt ( obj_sim_param%residuotot ) ! L2 norm time derivative
      endif
      deallocate ( res_block )
    endif
    
    ! Determine simulation control flags on root, then broadcast to all ranks. Only one of endsim
    ! (finest level) and endmg (coarse levels) is evaluated below; both are initialised at the top of
    ! the root block (uninitialised, they renumbered the output files run to run and could end a
    ! multigrid run early: plan 10, upstream 83562b5). Non-root ranks never read them.
    if (mpi_is_root) then
      endsim = .false.
      endmg  = .false.
      ! With time tables a time-accurate run goes on to time-threshold even when T stops changing.
      ! A receding surface (law = nodeshift, plan 08) is never steady: the residual measures the
      ! stage update only, not the remap, and a uniform flux-free state has none -- such a run
      ! ends on time or iterations.
      converged = ( obj_sim_param%residuotot <= obj_sim_param%res_threshold ) .and. .not. (advance .and. has_tdep_bc) &
                  .and. trim(obj_mesh_motion%law) /= 'nodeshift'
      if (level == 1) then
        endsim = ( obj_sim_param%iter_from_call >= domain(1) % itermax ) &
            .or. converged &
            .or. ( domain(1) % time >= obj_sim_param%time_threshold )
      else
        endmg = ( obj_sim_param%iter_from_call >= domain(level) % itermax ) &
            .or. converged &
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


  !> Residual of the current state of one level, with no update: ghost fill, interior and boundary
  !> fluxes (and the FAS forcing when the level carries one). Used by the FAS transfers.
  subroutine Evaluate_Residual ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Lib_Ghost,     only: Fill_Ghost_Cell
    use FUSS_Mod_Fluxes,    only: Fluxes
    use FUSS_Mod_BC_Fluxes, only: BC_Fluxes
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain

    !$omp parallel
    call Fill_Ghost_Cell ( domain )
    call Fluxes ( domain )
    call BC_Fluxes ( domain )
    !$omp end parallel
    if (domain % fas_forcing) call Add_FAS_Forcing ( domain )

  end subroutine Evaluate_Residual


  !> R <- R - fas_f: the volume-integrated FAS forcing enters the update as +fas_f/V
  !> (Scale_Residual_Cons applies -R/V + qvol), identically for cons and prim.
  subroutine Add_FAS_Forcing ( domain )
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
        domain % blk(b) % R (i,j,k) = domain % blk(b) % R (i,j,k) - domain % blk(b) % fas_f (i,j,k)
      enddo; enddo; enddo
      !$omp end do
    enddo

  end subroutine Add_FAS_Forcing


  subroutine Copy_State ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_pyrolysis, obj_darcy
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

    ! Start-of-step reaction state, the chi analogue of TO (plan 06 section 3.1)
    if ( obj_pyrolysis%enabled ) then
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3)
        do k = 1, domain % blk(b) % dim(3)
        do j = 1, domain % blk(b) % dim(2)
        do i = 1, domain % blk(b) % dim(1)
          domain % blk(b) % chiO (:,i,j,k) = domain % blk(b) % chi (:,i,j,k)
        enddo; enddo; enddo
        !$omp end do
      enddo
    endif

    ! Start-of-step gas pressure p^n and stored-gas coefficient (plan 07 section 3.3)
    if ( obj_darcy%enabled ) then
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3)
        do k = 1, domain % blk(b) % dim(3)
        do j = 1, domain % blk(b) % dim(2)
        do i = 1, domain % blk(b) % dim(1)
          domain % blk(b) % p_pgO (i,j,k)   = domain % blk(b) % p_pg (i,j,k)
          domain % blk(b) % eta_stO (i,j,k) = domain % blk(b) % eta_st (i,j,k)
        enddo; enddo; enddo
        !$omp end do
      enddo
    endif

  end subroutine Copy_State


  !> T^(n-1) of the Darcy closure (plan 07; plan 08 section 2.9 phase 0),
  !> taken from TO, which still holds the previous step's start temperature --
  !> not on the first step of a run, where TO is unset and T_prev keeps what
  !> setup gave it (T on a new run, the file's on a restart). Called before the
  !> remap, so that a moving mesh remaps this step's history, and before
  !> Copy_State overwrites TO. Inside the parallel region.
  subroutine Darcy_Shift_History ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_sim_param
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, j, k, b

    if ( obj_sim_param % iter_from_call <= 1 ) return   ! the first step of THIS process (a resumed restart continues domain%iter)
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        domain % blk(b) % T_prev (i,j,k) = domain % blk(b) % TO (i,j,k)
      enddo; enddo; enddo
      !$omp end do
    enddo

  end subroutine Darcy_Shift_History

end module FUSS_Mod_Explicit