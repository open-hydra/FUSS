module FUSS_Wrap_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: FUSS_setup

contains

  subroutine FUSS_setup ( simulation )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m
    use FUSS_Global_m
    use FUSS_Read_Ini,             only: Read_Inifile
    use FUSS_Load_Table,           only: Load_Table
    use FUSS_Assign_Setup,         only: Assign_Setup
    use FUSS_IO_Solution,          only: Read_IC, Setup_Output_Solution
    use FUSS_Mod_Allocate_Data,    only: Setup_Data_Structure, deallocate_remote_computation_data
    use FUSS_Mod_Multigrid,        only: Setup_Multigrid, Restriction
    use FUSS_IO_BC,                only: Setup_BC
    use FUSS_Mod_Metrics,          only: Setup_Metrics
    use FUSS_IO_Probes,            only: Setup_Probes, Restart_Probes, Drop_Probe_Rows, Drop_Rows_After
    use FUSS_IO_Wall,              only: Initialize_Wall_File, Set_Wall_Variables, IOwall
    use FUSS_Lib_Ghost,            only: Fill_Ghost_Cell, Fill_matIDg
    use FUSS_Mod_MPI,              only: mpi_is_root, partition_blocks, mpi_stop_all, mpi_stop_if_any, mpi_bcast_logical
    use FUSS_Load_Pyrolysis,       only: Load_Pyrolysis_Tables, Check_Pyrolysis_Tables
    use FUSS_Mod_Pyrolysis,        only: Write_Energy_Budget
    use FUSS_Mod_Darcy,            only: Darcy_Setup_Blocks, Fill_Ghost_p, Darcy_Update_Velocity, Darcy_Print_Setup, &
                                         Darcy_Is_Open
    use FUSS_Mod_Timers,           only: timer_run_begin
    use FUSS_Mod_GhostExchange,    only: build_ghost_schedule, build_local_bc_index, &
                                         allocate_exchange_schedules
    use FUSS_Mod_Recession,        only: Ablation_Update_Env, Setup_Ablation, Setup_Recession
    implicit none
    type(FUSS_simulation_type), intent(inout) :: simulation
    ! Local
    integer :: m, ios
    logical :: setup_error   ! an error found by a check of rank 0 (Stop_Together)
    logical :: resume        ! a restart goes on with the counters of its solution

    !! ------------------------------------------------------
    !! ------------------------------------------------------
    ! Print Header
    if (mpi_is_root) call Print_Header ()
    
    ! Load tables (thermo, transport)
    call Load_Table ()
    setup_error = .false.
    if (mpi_is_root) call Check_Table ( setup_error )
    call Stop_Together ( setup_error, 'material properties error (see above)' )

    ! Read input.ini
    call Read_Inifile ()
    ! Flags DERIVED from the ini must exist before Check_Input, or a guard that
    ! tests them is a comment for the input that sets them indirectly.
    ! Example: `irs-beta > 0` enables IRS, but the promotion lived in
    ! Assign_Setup, which runs after Check_Input, so `irs-beta = 0.5` with no
    ! `irs` key passed the mesh-motion IRS guard and smoothed the remapped
    ! residual anyway. Guard demo: run_guard_demos.sh irs_beta.
    ! irs-beta alone switches IRS (upstream 3447dd9: the `irs` key is gone)
    obj_irs%enabled = ( obj_irs%beta > 0d0 )
    ! Same reason, plan 07: the Darcy closure is selected by a [FUSS-Pyrolysis]
    ! key; every guard, the table checks and Assign_Setup read this flag.
    obj_darcy%enabled = obj_pyrolysis%enabled .and. ( trim(obj_pyrolysis%gas_closure) == 'type2' )
    setup_error = .false.
    if (mpi_is_root) call Check_Input ( setup_error )
    call Stop_Together ( setup_error, 'input file error (see above)' )
    call Check_FAS ()

    ! Pyrolysis material files (plan 06 section 3.2): loaded on every rank once the input is validated (a file every
    ! rank reads stops every rank together, mpi_stop_all), value-checked on rank 0
    if (obj_pyrolysis%enabled) call Load_Pyrolysis_Tables ()
    setup_error = .false.
    if (mpi_is_root .and. obj_pyrolysis%enabled) call Check_Pyrolysis_Tables ( setup_error )
    call Stop_Together ( setup_error, 'pyrolysis material error (see above)' )

    ! Allocate container objects Domain and IOfield.
    allocate ( simulation%domain  ( obj_multigrid%MGL ) )
    allocate ( simulation%IOfield ( obj_multigrid%MGL ) )
    simulation%domain(1)%itermax = obj_sim_param%iter_threshold

    ! Assign setup
    call Assign_Setup ()

    ! Read file for initial solution.
    call Read_IC ( simulation%IOfield(1) )
    ! The ranks decide together on an error of the read (Check_IC)
    call Check_IC ()
    ! A restart goes on from the time of the solution it reads: the physical time of a time-accurate solution, minus the
    ! global iteration count of a steady one. A solution that holds no time (a vtk file without TIME, a Tecplot ASCII
    ! file without SOLUTIONTIME) cannot tell where to go on: stop. Each rank reads the solution itself (with one folder
    ! per rank, its own file), so the ranks decide together: the first rank that found no time writes the message
    call mpi_stop_if_any ( .not. obj_sim_param%newrun .and. .not. obj_io%restart%time_found, &
      'restart from a solution without a time (see above)', &
      ' [ERROR] Restart from '//trim(obj_io%nameinit)//': the solution holds no time (vtk TIME, Tecplot SOLUTIONTIME), '// &
      'and a restart goes on from the time of its solution. To start from this solution, use it as the initial '// &
      'condition: INPUT/ic.* with newrun = true.' )

    ! Allocation of data structures and copy solution from IOfield.
    call Setup_Data_Structure ( simulation%domain(1), simulation%IOfield(1) )
    ! Every value that the initial condition gives the domain is finite, before the state is tested and the coarse
    ! levels are built from it
    call Check_IC_Values ( simulation%domain(1) )
    call Check_State ( simulation%domain(1) )

    ! With multigrid, allocate Grid and IOfield 2,...,MGL
    if ( obj_multigrid%MGL > 1 ) then
      call Setup_Multigrid ( simulation )
    end if

    ! Read Boundary Conditions file
    call Setup_BC ( simulation%domain )
    ! The checks of rank 0 on the boundary conditions, in the order of a run that stops at the first refusal: a later
    ! check reads what an earlier one refused (plan 14 S6, D22)
    setup_error = .false.
    if (mpi_is_root) then
      call Check_BC ( setup_error )
      ! Periodic faces (201): the combinations not built yet stop here (plan 14 S9d Part A, interim refusals R1-R3)
      if ( .not. setup_error ) call Check_Periodic_Compatibility ( simulation%domain(1), setup_error )
      if ( .not. setup_error .and. obj_pyrolysis%enabled ) &
        call Check_Pyrolysis_BC_Topology ( simulation%domain(1), setup_error )
      ! Ablating walls (plan 08): their rules hold with or without pyrolysis, so they are checked apart
      if ( .not. setup_error .and. obj_ablation%enabled ) call Check_Ablation ( simulation%domain(1), setup_error )
      if ( .not. setup_error .and. ( obj_ablation%enabled .or. trim(obj_mesh_motion%law) == 'nodeshift' ) ) &
        call Check_Ablation_BC_Topology ( simulation%domain(1), setup_error )
    end if
    call Stop_Together ( setup_error, 'boundary condition error (see above)' )
    ! A restart reads the ablating faces' state back; law = nodeshift: the receding faces, the surface-node map
    call Setup_Ablation ( simulation%domain(1), obj_sim_param%newrun )
    call Setup_Recession ( simulation%domain(1) )
    
    ! MPI: partition blocks once on finest-level cell counts (restriction and
    ! prolongation are block-local, so a block must live on the same rank at
    ! every level), then build per-level ghost/chimera communication schedules.
    ! The rank of every block comes from the map that decomposition-map names in INPUT/ (with the phase prefix of
    ! bc.txt), in the format of the map of the ATLAS mesh decomposer; with none, the default, from the LPT, and no file is
    ! read.
    if ( trim(obj_io%decomposition_map) == 'none' ) then
      call partition_blocks(simulation%domain(1)%nb, &
        [(product(simulation%domain(1)%blk(ios)%dim(1:3)), ios=1, simulation%domain(1)%nb)])
    else
      call partition_blocks(simulation%domain(1)%nb, &
        [(product(simulation%domain(1)%blk(ios)%dim(1:3)), ios=1, simulation%domain(1)%nb)], &
        reshape([(simulation%domain(1)%blk(ios)%dim(1:3), ios=1, simulation%domain(1)%nb)], &
                [3, simulation%domain(1)%nb]), &
        'INPUT/'//trim(FUSS_phase_prefix)//trim(obj_io%decomposition_map))
    end if

    call allocate_exchange_schedules(obj_multigrid%MGL)
    do m = 1, obj_multigrid%MGL
      simulation%domain(m)%mg_level = m
      call build_ghost_schedule(simulation%domain(m))
      call build_local_bc_index(simulation%domain(m))
      call fill_matIDg (simulation%domain(m))
    end do

    ! Darcy closure (plan 07 section 3.3): the per-block sweep direction and
    ! the sweep work arrays, once the BC types and the block ownership are known.
    if (obj_darcy%enabled) call Darcy_Setup_Blocks ( simulation%domain(1) )

    ! Read grid file and setup metrics in each grid level.
    do m = 1, obj_multigrid%MGL
      call Setup_Metrics ( simulation%domain(m) )
    end do
    
    ! Setup probes location, output files, and variables to be printed (before the
    ! nodes of remote blocks are freed: placing a probe searches every block)
    call Setup_Probes ( simulation%domain(1), obj_sim_param%newrun )

    ! Free heavy arrays on non-local blocks to save memory
    do m = 1, obj_multigrid%MGL
      call deallocate_remote_computation_data(simulation%domain(m))
    end do

 
    ! Solution setup
    call Setup_Output_Solution ( simulation%IOfield )

    ! Initialize the wall files, one per multigrid level (wall, wall-level<n>): the variables on
    ! every rank (they size the wall values' reduction, plan 11), the files and patch meshes on root
    allocate ( IOwall(obj_multigrid%MGL) )
    if (obj_io % write_wall) call Set_Wall_Variables ()
    if (mpi_is_root .and. obj_io % write_wall) then
      do m = 1, obj_multigrid%MGL
        call Initialize_Wall_File ( simulation%domain(m), simulation%IOfield(1)%tec%extension)
      enddo
    endif

    ! Print simulation onto the logfile/shell.
    if (mpi_is_root) call Print_Shell_Info ()

    ! Residual history file
    if (mpi_is_root) then
      open(newunit=obj_io%unitRES,file='OUTPUT/'//trim(FUSS_phase_prefix)//'residual-history.dat',status='unknown',form='formatted')
      if ( .not. obj_sim_param%newrun ) then
        ios = 0
        do while ( ios == 0 )
          read (obj_io%unitRES, *, iostat=ios)
        enddo
        backspace (obj_io%unitRES)
      endif
    end if

    ! Setting the format for writing residual_history file
    write(obj_io%unitRES_format,'(A11,I0,A7)') '(I8,E20.10,', 1, 'E20.10)'

    ! Initialize run variables depending on scheme
    select case (trim(obj_time_scheme%solver_type))

      case ('euler', 'RK2', 'RK3')

        obj_multigrid%MG_level = obj_multigrid%MGL

        do m = 1, obj_multigrid%MGL
          simulation%domain(m)%iter = 0

          if (obj_sim_param%newrun) then
            if (obj_time_scheme%time_accurate) then
              simulation%domain(m)%time = 0.0d0
            else
              simulation%domain(m)%time = -1.0d0
            endif
            obj_sim_param%iter_general = 0
          else
            ! The restart solution is read on the finest level: its time and iteration count hold for every level (the
            ! coarse IOfields are built from the fine one and hold no time)
            if (obj_time_scheme%time_accurate) then
              if (obj_io%restart%steady) then
                ! A transient started from a steady solution starts at t = 0
                simulation%domain(m)%time = 0.0d0
              else
                simulation%domain(m)%time = simulation%IOfield(1)%solutiontime
              endif
              obj_sim_param%iter_general = 0
            else
              ! The global iteration count of the solution: a steady one holds it as its time (Read_vtk_tec), a
              ! transient one holds its physical time, of which the integer part is taken
              if (obj_io%restart%steady) then
                obj_sim_param%iter_general = nint(simulation%IOfield(1)%solutiontime)
              else
                obj_sim_param%iter_general = int(simulation%IOfield(1)%solutiontime)
              endif
              simulation%domain(m)%time = -1.d0
            endif
          endif
          obj_sim_param%time_from_call = simulation%domain(m)%time 
          obj_sim_param%iter_from_call = 0
        enddo

        ! The counters stored with a solution, in every format: the restart goes on as the uninterrupted run would (vnn
        ! ramp, output cadence and numbering, residual history, probes, grid level), if that run was of the same kind,
        ! steady or time-accurate. Otherwise, or with a solution that holds no counters (written before them), the run
        ! starts its counters as a new one, numbers its outputs from the files in OUTPUT/ (Wrap_Postprocess), and
        ! starts the multigrid from the coarsest level.
        resume = .false.
        if ( .not. obj_sim_param%newrun ) &
          resume = obj_io%restart%found .and. ( obj_io%restart%steady .neqv. obj_time_scheme%time_accurate )
        obj_io%restart%resumed = resume
        if ( resume ) then
          obj_sim_param%iter_general   = obj_io%restart%cycle
          obj_multigrid%MG_level       = max(1, min(obj_multigrid%MGL, obj_io%restart%grid_level))
          ! The solution read is that of the finest level, and so is its iteration count
          simulation%domain(1)%iter    = obj_io%restart%iteration
          obj_sim_param%time_from_call = obj_io%restart%last_output_time
          obj_io%output_number         = obj_io%restart%output_number
          ! The unit of res-threshold, when the solution holds it: the first iteration of the restart keeps it
          obj_sim_param%res_reference  = obj_io%restart%res_reference
          obj_sim_param%res_restored   = obj_io%restart%res_found
          ! A run that went on past the solution read left in residual-history.dat and in the probe files the rows of
          ! the iterations after it, which the restart writes again: drop them
          if (mpi_is_root) then
            close(obj_io%unitRES)
            call Drop_Rows_After ( 'OUTPUT/'//trim(FUSS_phase_prefix)//'residual-history.dat', 3, &
                                   obj_io%restart%cycle, simulation%domain(1)%time )
            open(newunit=obj_io%unitRES,file='OUTPUT/'//trim(FUSS_phase_prefix)//'residual-history.dat', &
                 status='old',position='append',form='formatted')
          endif
          call Drop_Probe_Rows ( obj_io%restart%cycle, simulation%domain(1)%time )
        endif
        ! Without it (a solution written before it was stored, or a time-accurate one), a steady restart takes the
        ! residual of its own first iteration, as a new run (as MOSE)
        if ( .not. ( obj_sim_param%newrun .or. obj_time_scheme%time_accurate .or. obj_sim_param%res_restored ) ) &
          obj_io%warning_message = '[WARNING] The restart solution holds no residual reference of a steady run: '// &
            'res-threshold is relative to the residual of the first iteration of the restart.'
        ! Time-based probe rows fall on the multiples of their dtime: a restart that goes on in time writes the next one
        ! not reached yet; a new transient starts from the first, as a new run
        if ( .not. obj_sim_param%newrun ) &
          call Restart_Probes ( simulation%domain(1)%time, &
                                resume = .not. ( obj_io%restart%steady .and. obj_time_scheme%time_accurate ) )
        call Fill_Ghost_Cell ( simulation%domain(1) )
        ! Ablating walls (plan 08): the env values at the start time, so the
        ! first step and a restart see the table's p_w
        call Ablation_Update_Env ( simulation%domain(1), simulation%domain(1)%time, setup = .true. )
        ! Gas-pressure ghosts and the t = 0 cell velocity (plan 07): serial
        ! here, inside the parallel region of every later step.
        if (obj_darcy%enabled) then
          call Fill_Ghost_p ( simulation%domain(1) )
          call Darcy_Update_Velocity ( simulation%domain(1) )
        endif

        ! Interpolate initial solution on coarser domains
        if ( obj_multigrid%MGL > 1 ) then
          do m = 2, obj_multigrid%MGL
            call Restriction ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
            call Fill_Ghost_Cell ( simulation%domain(m) )
          end do
        endif

    end select

    ! Pyrolysis energy budget at t = t0, after the caches are built and the
    ! clock is set (plan 06 section 3.6); every rank (its block sums go to root)
    if (obj_pyrolysis%enabled) call Write_Energy_Budget ( simulation%domain(1) )

    ! Print warnings onto the logfile/shell.
    if (mpi_is_root) call Print_Warnings ()

    ! If errors are found, print error messages and stop the simulation (every rank together).
    setup_error = .false.
    if (mpi_is_root) call Stop_Simulation ( setup_error )
    call Stop_Together ( setup_error, 'setup error (see above)' )

    ! Set-up is over: start the solver clock
    call timer_run_begin()

  contains

    subroutine Print_Header() 
      
      write(*,*)
      write(*,'(A89)') '********************************************************************'
      write(*,'(A89)') '**                                                                **'
      write(*,'(A89)') '**  FFFFFFFFFFFF   UU         UU    SSSSSSSSSSSS    SSSSSSSSSSSS  **'
      write(*,'(A89)') '**  FF             UU         UU    SS              SS            **'
      write(*,'(A89)') '**  FF             UU         UU    SS              SS            **'
      write(*,'(A89)') '**  FFFFFFFF       UU         UU    SSSSSSSSSSSS    SSSSSSSSSSSS  **'
      write(*,'(A89)') '**  FF             UU         UU              SS              SS  **'
      write(*,'(A89)') '**  FF              UU       UU               SS              SS  **'
      write(*,'(A89)') '**  FF                UUUUUUU       SSSSSSSSSSSS    SSSSSSSSSSSS  **'   
      write(*,'(A89)') '**                                                                **'
      write(*,'(A89)') '**                Fourier Unsteady Solid Solver                   **'
      write(*,'(A89)') '**                                                                **'
      write(*,'(A89)') '********************************************************************'
      write(*,*)

    end subroutine Print_Header


    subroutine Check_Table ( has_error )
      implicit none
      logical, intent(out) :: has_error

      has_error = .false.

      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Loading'
      write(*,'(A)') ' ========================================================================================='

      ! Properties
      if (index(obj_table%error_message,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Properties', 'FAIL'
        write(*,'(4X,A)') trim(obj_table%error_message)
        has_error = .true.
      else
        write(*,'(A,T35,A)') '   Properties', 'OK'
      endif

    end subroutine Check_Table


    !> An error of the read of the initial condition (Read_IC, the reader writes it in obj_io%error_message) stops the
    !> run. Every rank reads the file, and a rank can fail to read what rank 0 reads: the ranks take the decision
    !> together (the largest of a flag); rank 0 writes the line of the set-up summary with the verdict of all ranks, the
    !> first rank that found an error writes it (with its rank when the run has more than one), and every rank stops
    !> with mpi_stop_all, exit code 1.
    subroutine Check_IC ()
      use FUSS_Mod_MPI, only: mpi_rank_, mpi_allreduce_max_r8, mpi_allreduce_min_r8, mpi_size_
      implicit none
      real(R8) :: flag, flag_all, first

      flag = 0.0_R8
      if ( index(obj_io%error_message,'ERROR') > 0 ) flag = 1.0_R8
      call mpi_allreduce_max_r8 ( flag, flag_all )
      if ( mpi_is_root ) then
        if ( flag_all > 0.0_R8 ) then
          write(*,'(A,T35,A)') '   Initial conditions', 'FAIL'
        else
          write(*,'(A,T35,A)') '   Initial conditions', 'OK'
        endif
      endif
      if ( flag_all > 0.0_R8 ) then
        call mpi_allreduce_min_r8 ( merge(real(mpi_rank_, R8), huge(1.0_R8), flag > 0.0_R8), first )
        if ( flag > 0.0_R8 .and. nint(first) == mpi_rank_ ) then
          if ( mpi_size_ > 1 ) then
            write(*,'(4X,A,I0,A)') trim(obj_io%error_message)//' (rank ', mpi_rank_, ')'
          else
            write(*,'(4X,A)') trim(obj_io%error_message)
          endif
        endif
        call mpi_stop_all ( 'initial condition error (see above)' )
      endif

    end subroutine Check_IC


    !> The initial condition, or the solution that a restart reads, gives the domain the coordinates of its nodes and the
    !> temperature, the material and the heat source of its cells, and a value there that is not finite (NaN, infinity)
    !> stops the run at set-up, with the place of the value. Not_Finite reads the bits of each value, an operation on
    !> integers that raises no floating-point exception: a test by comparison would stop a build that traps the invalid
    !> operations (the DEBUG build) with SIGFPE at the first such operation on the value, and -ffast-math cannot fold the
    !> test away. Each rank checks the values it holds; the ranks decide together (mpi_stop_if_any), the first rank that
    !> found a value writes it, and every rank stops: exit code 1, in every build.
    subroutine Check_IC_Values ( domain )
      use FUSS_Advanced_Types_m
      use FUSS_Finite_m, only: Not_Finite
      use FUSS_Mod_MPI,  only: mpi_stop_if_any
      implicit none
      type(FUSS_domain_type), intent(in) :: domain
      ! Local
      character(len=1), parameter :: axis(3) = ['x', 'y', 'z']
      character(len=hlen) :: what
      integer :: b, i, j, k, d

      what = ''
      blocks: do b = 1, domain%nb
        do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2); do i = 0, domain%blk(b)%dim(1)
          do d = 1, 3
            if ( Not_Finite( domain%blk(b)%node(i,j,k)%c(d) ) ) then
              write(what,'(3A,I0,A,2(I0,A),I0)') 'coordinate ', axis(d), ' of block ', b, ', node i, j, k = ', i, ', ', j, ', ', k
              exit blocks
            endif
          enddo
        enddo; enddo; enddo
        do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
          if ( Not_Finite( domain%blk(b)%T(i,j,k) ) ) then
            write(what,'(A,I0,A,2(I0,A),I0)') 'T of block ', b, ', cell i, j, k = ', i, ', ', j, ', ', k
            exit blocks
          elseif ( Not_Finite( domain%blk(b)%matID(i,j,k) ) ) then
            write(what,'(A,I0,A,2(I0,A),I0)') 'matID of block ', b, ', cell i, j, k = ', i, ', ', j, ', ', k
            exit blocks
          elseif ( Not_Finite( domain%blk(b)%qvol(i,j,k) ) ) then
            write(what,'(A,I0,A,2(I0,A),I0)') 'qvol of block ', b, ', cell i, j, k = ', i, ', ', j, ', ', k
            exit blocks
          endif
        enddo; enddo; enddo
      enddo blocks

      call mpi_stop_if_any ( len_trim(what) > 0, 'initial condition error (see above)', &
        '[ERROR] '//trim(obj_io%nameinit)//' holds a value that is not finite (NaN or infinity): '//trim(what) )

    end subroutine Check_IC_Values


    !> The initial field must be admissible: every matID a table index, every T inside the table. Each rank checks the
    !> field it read; the ranks decide together (mpi_stop_if_any), the first rank that found a cell writes it.
    subroutine Check_State ( domain )
      use FUSS_Advanced_Types_m
      use FUSS_Config_Types_m, only: obj_table
      use FUSS_Mod_MPI, only: mpi_stop_if_any
      use ir_precision, only: str
      implicit none
      type(FUSS_domain_type), intent(in) :: domain
      ! Local
      integer :: b, i, j, k, m
      real(R8) :: T
      logical :: found
      character(len=hlen) :: detail

      found = .false.
      detail = ''
      blocks: do b = 1, domain%nb
        do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
          m = nint ( domain%blk(b)%matID(i,j,k) )
          if (m < 1 .or. m > obj_table%n) then
            write(detail,'(A,4I6,A)') ' [ERROR] initial condition: matID '//trim(str(.true.,m))// &
              ' at b, i, j, k:', b, i, j, k, ' is not a material of properties.dat (1..'//trim(str(.true.,obj_table%n))//')'
            found = .true.
            exit blocks
          endif
          T = domain%blk(b)%T(i,j,k)
          if (.not. (T >= obj_table%Tmin .and. T < real(obj_table%Tmax,R8))) then
            write(detail,'(A,4I6,A,ES12.4,A)') ' [ERROR] initial condition at b, i, j, k:', b, i, j, k, &
              '  T =', T, ' K, outside the properties table ['//trim(str(.true.,obj_table%Tmin))//', '// &
              trim(str(.true.,obj_table%Tmax))//') K'
            found = .true.
            exit blocks
          endif
        enddo; enddo; enddo
      enddo blocks

      call mpi_stop_if_any ( found, 'initial condition outside the properties table', detail )

    end subroutine Check_State


    subroutine Check_Input ( has_error )
      use FUSS_Input_Registry
      use FUSS_Read_Ini, only: input_warnings, unknown_options_message, long_values_message
      implicit none
      logical, intent(out) :: has_error
      character(len=hlen) :: out
      integer :: r

      has_error = .false.
      ! Order (plan 14 S6, D-S6-2): Scan_Ini's errors first (levels that are not an integer, a probe numbered past a gap; upstream
      ! 189ba21): a probe past a gap is also an option nothing registered, and its own message is the useful one. Then ours'
      ! two checks (plan 10 D1, plan 14 D7), before the registry as in plan 10: a misspelt key keeps its default and passes
      ! the registry, so the unknown option is the cause to name. Then the registry (upstream's rules on every key).
      out = ''
      if ( allocated(reg%load_error) ) then
        out = reg%load_error
      else if ( len_trim(unknown_options_message) > 0 ) then
        has_error = .true.
        write(*,'(A,T35,A)') '   Input file', 'FAIL'
        write(*,'(4X,A)') trim(unknown_options_message)
        write(*,'(4X,A)') 'This part of input.ini is not being applied. A misspelt key keeps'
        write(*,'(4X,A)') 'its default, and a repeated section header is ignored entirely --'
        write(*,'(4X,A)') 'in both cases the run would otherwise look completely normal.'
      else if ( len_trim(long_values_message) > 0 ) then
        has_error = .true.
        write(*,'(A,T35,A)') '   Input file', 'FAIL'
        write(*,'(4X,A)') trim(long_values_message)
        write(*,'(4X,A)') 'The value would be cut to the limit without a message.'
      else
        out = Validate_Registry()
        ! decomposition-map names a file of INPUT/: the folder, and the phase prefix of bc.txt, are put before it
        if ( index(out,'ERROR') == 0 .and. trim(obj_io%decomposition_map) /= 'none' .and. &
             index(obj_io%decomposition_map, '/') > 0 ) &
          out = '[ERROR] input.ini, [FUSS-IO] decomposition-map = '//trim(obj_io%decomposition_map)// &
                ': the value is the name of a file in INPUT/, without a folder'
      endif
      if ( index(out,'ERROR') > 0 ) then
        has_error = .true.
        write(*,'(A,T35,A)') '   Input file', 'FAIL'
        write(*,'(4X,A)') trim(out)
      else if ( .not. has_error ) then
        write(*,'(A,T35,A)') '   Input file', 'OK'
      endif
      if (allocated(input_warnings)) then
        do r = 1, size(input_warnings)
          write(*,'(4X,A)') 'WARNING: '//trim(input_warnings(r))
        enddo
      endif

      ! Ours' rules between keys (plan 05/06/07/08): each runs only on an input that passed what precedes it, as when
      ! each of them stopped the run at its end
      if ( .not. has_error ) call Check_Input_Consistency ( has_error )
      if ( .not. has_error ) call Check_Mesh_Motion_Compatibility ( has_error )
      if ( .not. has_error ) call Check_Pyrolysis_Compatibility ( has_error )

    end subroutine Check_Input


    !> Keys that are each valid but contradict one another, or the build.
    subroutine Check_Input_Consistency ( has_error )
      implicit none
      logical, intent(out) :: has_error
      logical :: bad

      bad = .false.

#ifdef FUSS_NO_TECIO
      ! Without TecIO, ORION cannot write or read a binary Tecplot file and ends the
      ! run there with a bare stop: exit status 0, after the setup.
      if ( index(obj_io%sol_format, 'binary') > 0 .and. index(obj_io%sol_format, 'tecplot') > 0 ) then
        write(*,'(A)') ' [ERROR] sol-format = tecplot binary needs TecIO, and this FUSS is built without it'
        write(*,'(A)') '         (USE_TECIO = false). Use tecplot ascii, or rebuild with USE_TECIO = true.'
        bad = .true.
      endif
      if ( index(obj_io%ic_format, 'binary') > 0 .and. index(obj_io%ic_format, 'tecplot') > 0 ) then
        write(*,'(A)') ' [ERROR] ic-format = tecplot binary needs TecIO, and this FUSS is built without it'
        write(*,'(A)') '         (USE_TECIO = false). Use tecplot ascii, or rebuild with USE_TECIO = true.'
        bad = .true.
      endif
#endif

      has_error = bad

    end subroutine Check_Input_Consistency


    !> Reject combinations of moving mesh with convergence-acceleration features
    !> whose interaction has not been analysed.
    !>
    !> These are enforced as ERRORS rather than documented as caveats on
    !> purpose. Both would otherwise produce a plausible-looking answer that is
    !> quietly wrong, which is the worst possible failure mode for a solver:
    !>
    !>   * IRS smooths the residual AFTER the volumetric source has been folded
    !>     in and scaled by dt (Lib_Newstate.f90, Scale_Residual_Cons ->
    !>     Residual_Smoothing). With a moving mesh the residual also carries the
    !>     ALE swept-volume flux, and Jacobi-smoothing that term across cells
    !>     destroys the discrete GCL balance it was constructed to satisfy.
    !>
    !>   * Multigrid builds the metrics of every level once, at setup
    !>     (Wrap_Setup -> Setup_Metrics over all MGL levels). There is no
    !>     coarse-level mesh-motion law, so coarse levels would keep the initial
    !>     geometry while the fine level moves.
    !>
    !> Lifting either restriction is a separate piece of work, not a flag flip.
    subroutine Check_Mesh_Motion_Compatibility ( has_error )
      implicit none
      logical, intent(out) :: has_error
      logical :: bad

      ! Test `law` directly rather than obj_mesh_motion%enabled: Check_Input
      ! runs immediately after Read_Inifile, whereas `enabled` is derived later
      ! in Assign_Setup. Reading the flag here would silently test an unset
      ! value and the guard would never fire -- which is exactly what happened
      ! the first time this was written.
      has_error = .false.
      if ( trim(obj_mesh_motion%law) == 'static' ) return

      bad = .false.

      if ( obj_irs%enabled ) then
        write(*,'(A)') ' [ERROR] mesh motion is not compatible with implicit residual smoothing.'
        write(*,'(A)') '         IRS smooths the residual after the ALE swept-volume flux has'
        write(*,'(A)') '         been added, which breaks the discrete GCL balance.'
        write(*,'(A)') '         Set irs = false, or law = static.'
        bad = .true.
      endif

      if ( trim(obj_time_scheme%integration_variables) /= 'cons' ) then
        write(*,'(A)') ' [ERROR] mesh motion requires integration-variables = cons.'
        write(*,'(A)') '         The primitive path integrates T with a rho*cp coefficient and'
        write(*,'(A)') '         carries no ALE swept-volume term, so it would silently ignore'
        write(*,'(A)') '         the volume change instead of conserving energy across it.'
        bad = .true.
      endif

      ! A moving mesh in a non-time-accurate run is not merely inaccurate, it is
      ! undefined. Mod_Explicit advances domain%time only when one of the two
      ! time-accurate flags is set; otherwise time stays pinned at -1.0 and the
      ! motion law is evaluated at (-1 + dtglobal) every iteration. dtglobal is
      ! itself the evolving local-time-stepping minimum, so the mesh does not
      ! freeze -- it wobbles as convergence proceeds, and the swept volumes are
      ! whatever that wobble produced. The GCL still holds, so nothing would
      ! flag it.
      if ( .not. obj_time_scheme%time_accurate .and. &
           .not. obj_sim_param%HYDRA_time_accurate ) then
        write(*,'(A)') ' [ERROR] mesh motion requires a time-accurate run.'
        write(*,'(A)') '         With local time stepping the simulation clock is never advanced,'
        write(*,'(A)') '         so the motion law would be evaluated at a meaningless time that'
        write(*,'(A)') '         drifts with the local dt.'
        write(*,'(A)') '         Set time-accurate = true, or law = static.'
        bad = .true.
      endif

      ! The oscillatory prescribed law anchors its displacement to a REFERENCE
      ! mesh, which it captures from whatever geometry it is first handed. On a
      ! restart that geometry is the already-moved mesh out of the solution
      ! file, so the displacement would be applied a second time on top of
      ! itself and the run would continue from a different -- and physically
      ! meaningless -- geometry. Nothing downstream could detect it: the GCL
      ! would still hold, because the mesh is self-consistent, just wrong.
      !
      ! Restoring it properly means persisting the reference mesh, which is
      ! effort spent on a verification-only law. The laws that Phase 3 actually
      ! needs (translation, and node-shifting driven by the wall recession rate)
      ! displace INCREMENTALLY from the current mesh and so restart correctly
      ! with no extra state at all. Refuse the one case that cannot.
      if ( trim(obj_mesh_motion%law) == 'prescribed' .and. .not. obj_sim_param%newrun ) then
        write(*,'(A)') ' [ERROR] the prescribed mesh-motion law cannot be restarted.'
        write(*,'(A)') '         It displaces from a reference mesh captured at setup, which on a'
        write(*,'(A)') '         restart is the already-moved geometry, so the displacement would'
        write(*,'(A)') '         be applied twice.'
        write(*,'(A)') '         Use law = translation, which displaces incrementally, or start a'
        write(*,'(A)') '         new run.'
        bad = .true.
      endif

      if ( obj_multigrid%MGL > 1 ) then
        write(*,'(A)') ' [ERROR] mesh motion is not compatible with multigrid.'
        write(*,'(A)') '         Coarse-level metrics are built once at setup and there is no'
        write(*,'(A)') '         coarse-level motion law, so coarse grids would keep the'
        write(*,'(A)') '         initial geometry while the fine grid moves.'
        write(*,'(A)') '         Use a single grid level, or law = static.'
        bad = .true.
      endif

      has_error = bad

    end subroutine Check_Mesh_Motion_Compatibility


    !> Refuse the combinations pyrolysis (plan 06) has no treatment for. This
    !> is the single, per-phase-extensible home for every pyrolysis / Darcy /
    !> ablation setup guard: later phases add their refusals here, each with a
    !> run_guard_demos.sh entry that shows it firing. Deliberately NOT inside
    !> Check_Mesh_Motion_Compatibility, which returns early on a static mesh.
    subroutine Check_Pyrolysis_Compatibility ( has_error )
      use FUSS_Input_Registry, only: Key_Is_Set
      use FUSS_Register_Darcy, only: Parse_Open_Codes
      implicit none
      logical, intent(out) :: has_error
      logical :: bad, ok
      integer :: codes(8), ncodes
      character(len=64) :: badtok

      bad = .false.

      ! 12. [FUSS-Darcy] is read only under gas-closure = type2 (plan 07 §3.8);
      !     present with another closure, or without [FUSS-Pyrolysis], nothing
      !     would apply it -- the silent-non-application the input layer refuses.
      if ( obj_darcy%section_present .and. .not. obj_darcy%enabled ) then
        write(*,'(A)') ' [ERROR] [FUSS-Darcy] is read only with [FUSS-Pyrolysis] gas-closure = type2.'
        write(*,'(A)') '         Remove the section, or select the type-2 gas closure.'
        bad = .true.
      endif

      if ( .not. obj_pyrolysis%enabled ) then
        has_error = bad
        return
      endif

      ! 1. IRS. The promotion of irs-beta > 0 runs before Check_Input; the
      !    second clause is a belt against it moving back.
      if ( obj_irs%enabled .or. obj_irs%beta > 0d0 ) then
        write(*,'(A)') ' [ERROR] pyrolysis is not compatible with implicit residual smoothing.'
        write(*,'(A)') '         The smoothed residual would no longer match the kinetics increments'
        write(*,'(A)') '         it is integrated with. Set irs = false and irs-beta = 0.'
        bad = .true.
      endif

      ! 2. Multigrid: there is no coarse-level reaction state.
      if ( obj_multigrid%MGL > 1 ) then
        write(*,'(A)') ' [ERROR] pyrolysis is not compatible with multigrid.'
        write(*,'(A)') '         The reaction advancements are not restricted or prolongated.'
        write(*,'(A)') '         Use a single grid level.'
        bad = .true.
      endif

      ! 3. The stored energy E(T, chi) is the integrated variable (plan 06 R1).
      if ( trim(obj_time_scheme%integration_variables) /= 'cons' ) then
        write(*,'(A)') ' [ERROR] pyrolysis requires integration-variables = cons.'
        write(*,'(A)') '         The primitive path integrates T with a rho*cp coefficient that'
        write(*,'(A)') '         assumes rho = rho(T); with pyrolysis rho depends on the history.'
        bad = .true.
      endif

      ! 4. (the type-2 refusal of phase 2 was lifted by plan 07)

      ! 5. The kinetics integrate physical time.
      if ( .not. obj_time_scheme%time_accurate .and. .not. obj_sim_param%HYDRA_time_accurate ) then
        write(*,'(A)') ' [ERROR] pyrolysis requires a time-accurate run.'
        write(*,'(A)') '         The Arrhenius kinetics integrate physical time; local time stepping'
        write(*,'(A)') '         would advance the reactions of every cell by a different, meaningless dt.'
        write(*,'(A)') '         Set time-accurate = true.'
        bad = .true.
      endif

      ! 6. The pyrolysis state is remapped through the swept volumes only for
      !    the receding mesh of an ablating wall (plan 08 section 2.9).
      if ( trim(obj_mesh_motion%law) /= 'static' .and. trim(obj_mesh_motion%law) /= 'nodeshift' ) then
        write(*,'(A)') ' [ERROR] pyrolysis on a moving mesh needs law = nodeshift.'
        write(*,'(A)') '         Only the receding mesh of an ablating wall (306) remaps the pyrolysis state.'
        bad = .true.
      endif

      ! 7. The type-1 closure needs the direction the gas leaves through
      !    (type 2 takes the direction from the pressure field).
      if ( .not. obj_darcy%enabled .and. trim(obj_pyrolysis%outflow_direction) == 'none' ) then
        write(*,'(A)') ' [ERROR] pyrolysis: outflow-direction must be set (+i, -i, +j, -j, +k or -k).'
        write(*,'(A)') '         The type-1 closure integrates the gas production along that direction.'
        bad = .true.
      endif

      ! 8. (lifted by plan 11 stages 3a and 3b: pyrolysis runs on any rank count;
      !    the ablating wall keeps its own refusal in Check_Ablation until stage 4)

      ! Darcy closure, ini-level rules (plan 07 §3.8).
      if ( obj_darcy%enabled ) then
        ! 13. p-gas prescribes the gas pressure of type 1; under type 2 the
        !     pressure is solved, so an explicit p-gas contradicts the choice.
        if ( Key_Is_Set( trim(codename)//'-Pyrolysis', 'p-gas' ) ) then
          write(*,'(A)') ' [ERROR] [FUSS-Pyrolysis] p-gas is refused with gas-closure = type2: the gas pressure is solved.'
          write(*,'(A)') '         Remove p-gas (use [FUSS-Darcy] p-init / p-wall).'
          bad = .true.
        endif
        ! 14. The Jacobi damping must not exceed 1 (the registry checks one-sided rules only).
        if ( obj_darcy%pseudo_vnn > 1.0d0 ) then
          write(*,'(A,ES10.3)') ' [ERROR] [FUSS-Darcy] pseudo-vnn must be in (0, 1]; given ', obj_darcy%pseudo_vnn
          bad = .true.
        endif
        ! 15. Every open-bc-types token is an openable wall code.
        call Parse_Open_Codes( obj_darcy%open_bc_types, codes, ncodes, ok, badtok )
        if ( .not. ok ) then
          write(*,'(A)') ' [ERROR] [FUSS-Darcy] open-bc-types: bad token "'//trim(badtok)//'"'
          write(*,'(A)') '         Allowed: blank-separated codes from {301, 302, 303, 304, 305, 306, 103}.'
          bad = .true.
        endif
        ! 16. The grid-following material frame is a placeholder key.
        if ( obj_darcy%axes_follow_grid ) then
          write(*,'(A)') ' [ERROR] [FUSS-Darcy] axes-follow-grid = true is not built: the permeability tensor is given in the mesh frame.'
          bad = .true.
        endif
      endif

      has_error = bad

    end subroutine Check_Pyrolysis_Compatibility


    !> bc.txt-level rules of plan 06 sections 2.1 and 3.7 (root, after
    !> Check_BC): one block, no 101/102/103/104 face (chi ghosts are plan 08) and
    !> a wall BC (301-305) on the outflow end of the block.
    subroutine Check_Pyrolysis_BC_Topology ( domain, has_error )
      implicit none
      type(FUSS_domain_type), intent(in) :: domain
      logical, intent(out) :: has_error
      logical :: bad, kinetics_on
      integer :: i, fo, r, nopen

      bad = .false.

      ! Darcy closure, bc-level rules (plan 07 §3.4, §3.8). The chimera
      ! refusal must outlive plan 08's lifting of the single-block rule below.
      if ( obj_darcy%enabled ) then
        do i = 1, domain%nbound
          if ( domain%bc(i)%type == 102 .or. domain%bc(i)%type == 104 ) then
            write(*,'(A)') ' [ERROR] Darcy: a chimera (102, 104) face has no pressure condition under gas-closure = type2.'
            write(*,'(A,I0,A,I0,A,I0)') '         bc.txt record ', i, ' is of type ', domain%bc(i)%type, ' on face ', domain%bc(i)%f
            bad = .true.
            exit
          endif
        enddo
        nopen = 0
        do i = 1, domain%nbound
          if ( Darcy_Is_Open( domain%bc(i)%type ) ) nopen = nopen + 1
        enddo
        kinetics_on = ( obj_darcy%pi_source > 0.0d0 )
        do r = 1, obj_pyrolysis%nsub
          if ( obj_pyrolysis%reac(r)%A > 0.0d0 ) kinetics_on = .true.
        enddo
        if ( nopen == 0 .and. kinetics_on ) then
          write(*,'(A)') ' [ERROR] Darcy: no face is open to the gas while the pyrolysis gas is produced.'
          write(*,'(A)') '         The gas has nowhere to go: open a wall face through [FUSS-Darcy] open-bc-types.'
          bad = .true.
        endif
        if ( nopen == 0 .and. obj_darcy%quasi_steady ) then
          write(*,'(A)') ' [ERROR] Darcy: a sealed domain has no unique quasi-steady pressure (singular line systems).'
          write(*,'(A)') '         Open a face through [FUSS-Darcy] open-bc-types, or set quasi-steady = false.'
          bad = .true.
        endif
      endif

      ! Guard 9 (plan 08 section 2.9): block connections (101) under type 2
      ! only -- the type-1 outflow line integral restarts at every block, so a
      ! second block's gas would be lost; chimera (102, 104) has no chi or p
      ! donor path, and a coupled face (103) no chi exchange (the coupled
      ! ablating wall is 307).
      if ( .not. obj_darcy%enabled .and. domain%nb > 1 ) then
        write(*,'(A)') ' [ERROR] pyrolysis: gas-closure = type1 is single-block (its outflow line integral does not '// &
          'cross block faces); use type2 for a multi-block mesh.'
        write(*,'(A,I0,A)') '         The mesh has ', domain%nb, ' blocks.'
        bad = .true.
      endif
      do i = 1, domain%nbound
        select case ( domain%bc(i)%type )
          case (101)
            if ( obj_darcy%enabled ) cycle
            if ( domain%bc(i)%per /= 0 ) then
              ! a periodic face (201) is read as a 101 with per /= 0 (IO_BC, upstream 477e50b)
              write(*,'(A)') ' [ERROR] pyrolysis: gas-closure = type1 does not take a periodic face (201): FUSS reads it '// &
                'as a block connection (101), and the outflow line integral does not cross block faces.'
              write(*,'(A,I0,A,I0,A,I0)') '         bc.txt record ', i, ' is a periodic face (201) on face ', &
                domain%bc(i)%f, ' of block ', domain%bc(i)%b
            else
              write(*,'(A)') ' [ERROR] pyrolysis: gas-closure = type1 is single-block (its outflow line integral does not '// &
                'cross block faces); use type2 for a multi-block mesh.'
              write(*,'(A,I0,A,I0)') '         bc.txt record ', i, ' is a block connection (101) on face ', domain%bc(i)%f
            endif
            bad = .true.
            exit
          case (102, 103, 104)
            write(*,'(A)') ' [ERROR] pyrolysis: 102, 103 and 104 faces are refused (a 101 face needs gas-closure = type2).'
            write(*,'(A,I0,A,I0,A,I0)') '         bc.txt record ', i, ' is of type ', domain%bc(i)%type, &
              ' on face ', domain%bc(i)%f
            bad = .true.
            exit
        end select
      enddo

      select case ( trim(obj_pyrolysis%outflow_direction) )
        case ('+i') ; fo = 2
        case ('-i') ; fo = 1
        case ('+j') ; fo = 4
        case ('-j') ; fo = 3
        case ('+k') ; fo = 6
        case ('-k') ; fo = 5
        case default; fo = 0
      end select
      if ( obj_darcy%enabled ) fo = 0     ! guard 10 is the type-1 line integral's
      if ( fo > 0 ) then
        do i = 1, domain%nbound
          if ( domain%bc(i)%b /= 1 .or. domain%bc(i)%f /= fo ) cycle
          if ( domain%bc(i)%type < 301 .or. domain%bc(i)%type > 307 ) then
            write(*,'(A)') ' [ERROR] pyrolysis: outflow face must be a wall BC (301-307).'
            write(*,'(A,I0,A,I0,A)') '         Face ', fo, ' of block 1 carries type ', domain%bc(i)%type, &
              ': a symmetry plane or a null face cannot pass the pyrolysis gas.'
            bad = .true.
            exit
          endif
        enddo
      endif

      has_error = bad

    end subroutine Check_Pyrolysis_BC_Topology


    !> Periodic faces (201, upstream 477e50b): the set-up refusals of the combinations whose periodic treatment is not
    !> built yet (plan 14 S9d Part A, interim; root, after Check_BC). Each refusal is one block between two marker lines,
    !> R1 to R3, which Part B removes with its guard demos (commits B1 to B3). A 201 is read as a 101 with per = -1
    !> (IO_BC); Setup_Periodic numbers the pairs later, inside Setup_Metrics, so a 101 with per /= 0 is the test here.
    subroutine Check_Periodic_Compatibility ( domain, has_error )
      implicit none
      type(FUSS_domain_type), intent(in) :: domain
      logical, intent(out) :: has_error
      integer :: i, ip

      has_error = .false.
      ip = 0
      do i = 1, domain%nbound
        if ( domain%bc(i)%type == 101 .and. domain%bc(i)%per /= 0 ) then
          ip = i
          exit
        endif
      enddo
      if ( ip == 0 ) return

      ! S9d interim refusal R1 (Part A; removed by Part B commit B1)
      ! law = translation or prescribed: not supported yet in this commit. The rigid motion x -> R x + t of each pair is
      ! fitted once (Setup_Periodic) and the per-step refresh of the ghost metrics (Update_Mesh) keeps it, which is right
      ! for a motion that commutes with the pair (d = R d_s at the image nodes: any translation of a translational pair,
      ! a translation along the axis of a rotational one); Part B allows those and refuses the others, which make the
      ! problem itself non-periodic.
      if ( trim(obj_mesh_motion%law) == 'translation' .or. trim(obj_mesh_motion%law) == 'prescribed' ) then
        write(*,'(A)') ' [ERROR] mesh motion: law = '//trim(obj_mesh_motion%law)//' does not yet move a periodic face (201).'
        write(*,'(A,I0,A,I0,A,I0,A)') '         bc.txt record ', ip, ' is a periodic face (201) on face ', &
          domain%bc(ip)%f, ' of block ', domain%bc(ip)%b, '. Nothing checks yet that a step keeps the two faces'
        write(*,'(A)') '         of a pair images of each other by the rigid motion fitted at set-up.'
        has_error = .true.
      endif
      ! end of S9d interim refusal R1

      ! S9d interim refusal R2 (Part A; removed by Part B commit B2)
      ! law = nodeshift: not supported yet in this commit. The surface-node map joins block-boundary nodes by equal
      ! coordinates (Setup_Recession), so the
      ! two faces of a pair, one period apart, would recede apart, and Assert_Interface_Nodes (Mod_Recession), which
      ! checks every 101 face for bitwise shared nodes, would stop the first step ("the mesh tore at a block interface").
      if ( trim(obj_mesh_motion%law) == 'nodeshift' ) then
        write(*,'(A)') ' [ERROR] mesh motion: law = nodeshift does not yet take a periodic face (201).'
        write(*,'(A,I0,A,I0,A,I0,A)') '         bc.txt record ', ip, ' is a periodic face (201) on face ', &
          domain%bc(ip)%f, ' of block ', domain%bc(ip)%b, '. The receding nodes are not joined across a pair'
        write(*,'(A)') '         (one period apart), and its interface check would take the pair for a torn interface.'
        has_error = .true.
      endif
      ! end of S9d interim refusal R2

      ! S9d interim refusal R3 (Part A; removed by Part B commit B3)
      ! gas-closure = type2: not supported yet in this commit. Connection_Terms (Mod_Darcy) takes the gas permeability of
      ! the cell across a pair unrotated
      ! (right for a translational pair, or when K_v, K_c and B are invariant under the pair's rotation), and no
      ! periodic Darcy case has been run.
      if ( obj_darcy%enabled ) then
        write(*,'(A)') ' [ERROR] Darcy: a periodic face (201) is not yet supported under gas-closure = type2.'
        write(*,'(A,I0,A,I0,A,I0,A)') '         bc.txt record ', ip, ' is a periodic face (201) on face ', &
          domain%bc(ip)%f, ' of block ', domain%bc(ip)%b, '. The gas permeability of the cell across the pair'
        write(*,'(A)') '         is taken unrotated (wrong for a rotational pair unless K_v, K_c and B are invariant'
        write(*,'(A)') '         under its rotation), and no periodic Darcy case has been run.'
        has_error = .true.
      endif
      ! end of S9d interim refusal R3

    end subroutine Check_Periodic_Compatibility


    !> bc.txt-level rules of the ablating wall 306 (plan 08 sections 2.1, 2.3,
    !> 2.10), with or without pyrolysis (root, after Check_BC).
    subroutine Check_Ablation ( domain, has_error )
      implicit none
      type(FUSS_domain_type), intent(in) :: domain
      logical, intent(out) :: has_error
      logical :: bad
      integer :: i, r, m, fo
      real(R8) :: v(8)

      bad = .false.
      ! (the MPI refusal is lifted: plan 11 stage 4)
      if ( mod( obj_ablation%smoothing_width, 2 ) == 0 ) then
        write(*,'(A,I0,A)') ' [ERROR] [FUSS-Ablation] smoothing-width = ', obj_ablation%smoothing_width, &
          ' must be odd: the stencil is centred on its node (1 = no smoothing).'
        bad = .true.
      endif
      ! multigrid: refused where the first 306 record is read (Read_BCfile)
      if ( obj_mesh_motion%enabled ) then
        select case ( trim(obj_mesh_motion%law) )
          case ( 'translation', 'prescribed' )
            write(*,'(A)') ' [ERROR] the ablating wall (306) cannot ride a '//trim(obj_mesh_motion%law)// &
              ' mesh: its recession moves the mesh (law = nodeshift, plan 08), or nothing does (law = static).'
            bad = .true.
        end select
      endif

      fo = 0
      if ( obj_pyrolysis%enabled .and. .not. obj_darcy%enabled ) then
        select case ( trim(obj_pyrolysis%outflow_direction) )
          case ('+i') ; fo = 2
          case ('-i') ; fo = 1
          case ('+j') ; fo = 4
          case ('-j') ; fo = 3
          case ('+k') ; fo = 6
          case ('-k') ; fo = 5
        end select
      endif

      do i = 1, domain%nbound
        if ( domain%bc(i)%type /= 306 .and. domain%bc(i)%type /= 307 ) cycle
        associate ( bc => domain%bc(i) )
          ! the type-1 line integral carries the gas to the outflow face only
          if ( fo > 0 .and. ( bc%b /= 1 .or. bc%f /= fo ) ) then
            write(*,'(A,I0,A,I0,A)') ' [ERROR] bc.txt record ', i, ': under gas-closure = type1 a ', bc%type, &
              ' face must be the outflow face (outflow-direction): the gas it blows is the line integral that ends there.'
            bad = .true.
          endif
          ! a coupled ablating wall runs with the chemistry on (plan 08 section 2.5 ii);
          ! its exchange values have no table here (hydra, or Setup_Ablation's coupled-env-file)
          if ( bc%type == 307 ) then
            if ( .not. obj_pyrolysis%enabled ) then
              write(*,'(A,I0,A)') ' [ERROR] bc.txt record ', i, ': a coupled ablating wall (307) runs with the '// &
                'chemistry on, which needs [FUSS-Pyrolysis] (its gas and its B'' table).'
              bad = .true.
            endif
            if ( bad ) exit
            cycle
          endif
          m = 0
          if ( obj_pyrolysis%enabled ) m = nint( domain%blk(bc%b)%matID(bc%i,bc%j,bc%k) )
          do r = 1, bc%env%nrow
            v(1) = bc%env%time(r)
            v(2:8) = bc%env%var(:,r)
            if ( v(2) <= 0.0d0 .or. v(3) < 0.0d0 .or. v(7) < 0.0d0 .or. v(6) < 0.0d0 .or. v(6) > 1.0d0 &
                 .or. v(8) <= 0.0d0 ) then
              write(*,'(A,I0,A,I0,A)') ' [ERROR] bc.txt record ', i, ': env row ', r, &
                ' needs p_w > 0, rhoeUeCH0 >= 0, hconv >= 0, 0 <= chem <= 1 and T_edge > 0'
              bad = .true.
            endif
            if ( v(6) > 0.0d0 .and. .not. obj_pyrolysis%enabled ) then
              write(*,'(A,I0,A,I0,A)') ' [ERROR] bc.txt record ', i, ': env row ', r, ' has chem > 0, but a 306 face '// &
                'without [FUSS-Pyrolysis] (kinematic mode) has no pyrolysis gas and no B'' table: every row needs chem = 0.'
              bad = .true.
            endif
            if ( v(6) > 0.0d0 .and. bc%ibprime == 0 ) then
              write(*,'(A,I0,A,I0,A)') ' [ERROR] bc.txt record ', i, ': env row ', r, &
                ' has chem > 0 and the face has no B'' table (none).'
              bad = .true.
            endif
            if ( m > 0 .and. obj_darcy%enabled ) then
              associate ( mt => obj_pyrolysis%mat(m) )
                if ( v(2) < mt%p_lev(1) .or. v(2) > mt%p_lev(mt%Np) ) then
                  write(*,'(A,I0,A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A)') ' [ERROR] bc.txt record ', i, ': env row ', r, &
                    ' p_w = ', v(2), ' Pa lies outside the gas table of the face''s material [', mt%p_lev(1), ', ', &
                    mt%p_lev(mt%Np), '] Pa'
                  bad = .true.
                endif
              end associate
            endif
            if ( bad ) exit
          enddo
        end associate
        if ( bad ) exit
      enddo

      has_error = bad

    end subroutine Check_Ablation


    !> bc.txt-level rules of a RECEDING mesh (law = nodeshift, plan 08 section
    !> 2.3), with or without pyrolysis (root, after Check_Ablation). The node-
    !> map rules (a) and (b) need the surface-node map and are Setup_Recession's.
    subroutine Check_Ablation_BC_Topology ( domain, has_error )
      implicit none
      type(FUSS_domain_type), intent(in) :: domain
      logical, intent(out) :: has_error
      logical :: bad
      integer :: i, b, f, n306
      integer, allocatable :: face_of(:)
      real(R8) :: m0

      bad = .false.
      has_error = .false.
      if ( trim(obj_mesh_motion%law) /= 'nodeshift' ) return

      n306 = count( [( domain%bc(i)%type == 306, i = 1, domain%nbound )] )
      if ( n306 == 0 ) then
        write(*,'(A)') ' [ERROR] law = nodeshift moves the ablating faces (306), and bc.txt has none.'
        has_error = .true.
        return
      endif

      ! one receding face per block; its opposite face is not 306 as well
      allocate( face_of(domain%nb) )
      face_of = 0
      do i = 1, domain%nbound
        if ( domain%bc(i)%type /= 306 ) cycle
        b = domain%bc(i)%b ; f = domain%bc(i)%f
        if ( face_of(b) == 0 ) then
          face_of(b) = f
        elseif ( face_of(b) /= f ) then
          if ( ( face_of(b) + 1 ) / 2 == ( f + 1 ) / 2 ) then
            write(*,'(A,I0,A,I0,A,I0,A)') ' [ERROR] block ', b, ': its 306 faces ', face_of(b), ' and ', f, &
              ' are opposite faces; under law = nodeshift the face opposite the receding one stays put.'
          else
            write(*,'(A,I0,A,I0,A,I0,A)') ' [ERROR] block ', b, ' carries 306 faces on two face ids (', face_of(b), &
              ' and ', f, '); under law = nodeshift a block recedes from one face only.'
          endif
          bad = .true.
          exit
        endif
      enddo

      ! a receding block: no coupled face (MOSE does not move), one material, no qvol
      do i = 1, domain%nbound
        b = domain%bc(i)%b
        if ( face_of(b) == 0 ) cycle
        select case ( domain%bc(i)%type )
          case ( 103, 104, 307 )
            write(*,'(A,I0,A,I0,A,I0,A)') ' [ERROR] block ', b, ' recedes (law = nodeshift) and owns a coupled face ', &
              domain%bc(i)%type, ' (bc.txt record ', i, '): the coupled solver''s mesh does not move.'
            bad = .true.
            exit
        end select
      enddo
      do b = 1, domain%nb
        if ( face_of(b) == 0 ) cycle
        associate ( blk => domain%blk(b) )
          m0 = blk%matID(1,1,1)
          if ( any( blk%matID(1:blk%dim(1),1:blk%dim(2),1:blk%dim(3)) /= m0 ) ) then
            write(*,'(A,I0,A)') ' [ERROR] block ', b, ' recedes (law = nodeshift) and holds more than one material: '// &
              'matID is not remapped, so a receding block must be single-material.'
            bad = .true.
          endif
          if ( any( blk%qvol(1:blk%dim(1),1:blk%dim(2),1:blk%dim(3)) /= 0.0d0 ) ) then
            write(*,'(A,I0,A)') ' [ERROR] block ', b, ' recedes (law = nodeshift) and has a volumetric source qvol: '// &
              'qvol stays attached to the cells and is not remapped.'
            bad = .true.
          endif
        end associate
      enddo

      has_error = bad

    end subroutine Check_Ablation_BC_Topology


    subroutine Check_FAS()
      implicit none

      if (trim(obj_multigrid%cycle) == 'none') return
      if (obj_multigrid%MGL < 2) &
        call mpi_stop_all('[FUSS-Multigrid] cycle needs levels >= 2')
      if (obj_time_scheme%time_accurate .or. obj_sim_param%HYDRA_time_accurate) &
        call mpi_stop_all('[FUSS-Multigrid] cycle works on steady runs only (time-accurate = false)')
      if (obj_sim_param%HYDRA_MG) &
        call mpi_stop_all('[FUSS-Multigrid] cycle is not available under the hydra coupling')
      if (trim(obj_time_scheme%integration_variables) == 'prim') &
        call mpi_stop_all('[FUSS-Multigrid] cycle needs integration-variables = cons')
      if (obj_multigrid%pre_iter + obj_multigrid%post_iter < 1) &
        call mpi_stop_all('[FUSS-Multigrid] pre-iter + post-iter must be at least 1')

    end subroutine Check_FAS


    subroutine Check_BC ( has_error )
      implicit none
      logical, intent(out) :: has_error

      has_error = index(obj_io_bc%error_message,'ERROR') > 0
      if ( has_error ) then
        write(*,'(A,T35,A)') '   Boundary conditions', 'FAIL'
        write(*,'(4X,A)') trim(obj_io_bc%error_message)
      else
        write(*,'(A,T35,A)') '   Boundary conditions', 'OK'
      endif
      write(*,'(A)') ' ========================================================================================='

    end subroutine Check_BC


    subroutine Print_Warnings
      implicit none

      ! IO warnings
      if (index(obj_io%warning_message,'WARNING')>0) write(*,'(A)') obj_io%warning_message
      ! Material table warnings
      if (index(obj_table%warning_message,'WARNING')>0) write(*,'(A)') trim(obj_table%warning_message)
      ! Numerical scheme warnings
      if (index(obj_time_scheme%warning_message,'WARNING')>0) write(*,'(A)') obj_time_scheme%warning_message
      if (index(obj_irs%warning_message,'WARNING')>0) write(*,'(A)') obj_irs%warning_message

    end subroutine Print_Warnings


    subroutine Stop_Simulation ( has_error )
      implicit none
      logical, intent(out) :: has_error

      has_error = .false.

      ! Numerical scheme errors
      if (index(obj_time_scheme%error_message,'ERROR')>0) then
        write(*,'(A)') obj_time_scheme%error_message;  has_error = .true.
      endif
      if (index(obj_irs%error_message,'ERROR')>0) then
        write(*,'(A)') obj_irs%error_message;          has_error = .true.
      endif

    end subroutine Stop_Simulation


    !> A check of the set-up that runs on rank 0 alone writes what it finds there. Rank 0 broadcasts whether it found an
    !> error, and on an error every rank stops together (mpi_stop_all): rank 0 has written its message before any rank
    !> ends, and the run ends with exit code 1. With MPI_ABORT from rank 0 alone, which makes only a best attempt to stop
    !> the processes (MPI-4.1, section 11.4.2), the launcher can end while rank 0 is writing, and the other ranks go on
    !> meanwhile.
    subroutine Stop_Together ( found, why )
      implicit none
      logical, intent(inout)       :: found
      character(len=*), intent(in) :: why

      call mpi_bcast_logical ( found )
      if ( found ) call mpi_stop_all ( why )

    end subroutine Stop_Together


    subroutine Print_Shell_Info()
      use IR_Precision,   only: str  
      use FUSS_IO_BC,     only: Print_BC_Summary
      use FUSS_IO_Probes, only: nprobes
      use FUSS_Lib_Metrics, only: delthe, skew_factor_max
      implicit none
      ! Local
      character(llen) :: eosword
      character(len=16) :: angle 
      integer :: b, total_cells, k

      ! ----- Domain topology -----
      total_cells = 0
      do b = 1, simulation%domain(1)%nb
        total_cells = total_cells + &
          product(simulation%domain(1)%blk(b)%dim(1:3))
      end do

      write(*,*)
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Set-up'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Domain'
      write(*,'(A,T35,I0)') '   Blocks', simulation%domain(1)%nb
      write(*,'(A,T35,I0)') '   Cells', total_cells
      write(*,'(A,T35,I0)') '   Boundary faces', simulation%domain(1)%nbound
      ! Mesh type as Check_Mesh_Type takes it from the nodes of every block (delthe = 0: plane;
      ! FUSS has no type 200, so a 3-D mesh has no angle) and the largest time-step skew factor of the cells (1 unless 3-D
      ! cells are skewed in all pairs)
      write(angle,'(F16.6)') delthe * 180d0 / acos ( -1d0 )
      select case ( ndir )
        case ( 1 )
          write(*,'(A,T35,A)') '   Mesh type', '1-D'
        case ( 2 )
          if ( delthe >= 0d0 .and. delthe <= 0d0 ) then
            write(*,'(A,T35,A)') '   Mesh type', '2-D, plane'
          else
            write(*,'(A,T35,A)') '   Mesh type', '2-D, axisymmetric wedge of '//trim(adjustl(angle))//' degrees'
          end if
        case default
          write(*,'(A,T35,A)') '   Mesh type', '3-D'
      end select
      write(*,'(A,T35,F0.6)') '   Time-step skew factor (max)', skew_factor_max
      if ( obj_multigrid%MGL > 1 ) &
        write(*,'(A,T35,I0)') '   Multigrid levels', obj_multigrid%MGL

      ! ----- Boundary conditions -----
      call Print_BC_Summary ()

      ! ----- IO -----
      write(*,*)
      write(*,'(A)') ' Input/Output'
      write(*,'(A,T35,A)') '   Initial conditions file', trim(obj_io%nameinit)
      write(*,'(A,T35,A)') '   Solution format', trim(obj_io%sol_format)
      if ( .not. obj_sim_param%newrun ) then
        write(*,'(A,T35,A)') '   Restart from', trim(obj_io%nameinit)
        if ( obj_io%restart%found .and. ( obj_io%restart%steady .neqv. obj_time_scheme%time_accurate ) ) then
          write(*,'(A,T35,A)') '   Restart counters', 'restored from the solution'
        elseif ( obj_io%restart%found ) then
          write(*,'(A,T35,A)') '   Restart counters', 'not restored (the solution is of the other kind of run)'
        else
          write(*,'(A,T35,A)') '   Restart counters', 'not in the solution: counted as a new run'
        endif
      endif
      write(*,'(A,T35,A)') '   Probes number', str(.true.,nprobes)
      write(*,'(A)') ' ========================================================================================='

      ! ----- Physical model -----
      write(*,*)
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Physical model'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Solid model'
      write(*,'(A,T35,A)') ' Equation', trim(obj_sim_param%description)
      if ( obj_pyrolysis%enabled ) then
        write(*,'(A)') ' Pyrolysis'
        write(*,'(A,T35,A)') '   Gas closure', trim(obj_pyrolysis%gas_closure)
        if ( .not. obj_darcy%enabled ) then
          write(*,'(A,T35,A)') '   Outflow direction', trim(obj_pyrolysis%outflow_direction)
          write(*,'(A,T35,ES12.4)') '   Gas pressure [Pa]', obj_pyrolysis%p_gas
        endif
        write(*,'(A,T35,I0)') '   Materials', obj_pyrolysis%nmat
        write(*,'(A,T35,I0)') '   Reactions (all materials)', obj_pyrolysis%nsub
        if ( obj_darcy%enabled ) call Darcy_Print_Setup ( simulation%domain(1) )
      endif
      write(*,'(A)') ' ========================================================================================='

      ! ----- Numerical scheme -----
      write(*,*)
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Numerical scheme'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A,T35,A)') ' Space'
      write(*,'(A,T35,A)') '   Diffusive face metric', trim(obj_space_scheme%diffusive_metric)
      write(*,*)
      write(*,'(A,T35,A)') ' Time'
      write(*,'(A,T35,A)') '   Scheme', trim(obj_time_scheme%description)
      select case (obj_time_scheme%integration_variables)
        case ('prim')
          write(*,'(A,T35,A)') '   Integration variables', 'Primitive'
        case default
          write(*,'(A,T35,A)') '   Integration variables', 'Conservative'
      end select
      if (len_trim(obj_irs%description)>0) &
        write(*,'(A,T35,A)') '   Implicit residual smoothing', trim(obj_irs%description)
      write(*,'(A)') ' ========================================================================================='
      write(*,*)

    end subroutine Print_Shell_Info

  end subroutine FUSS_setup

end module FUSS_Wrap_Setup