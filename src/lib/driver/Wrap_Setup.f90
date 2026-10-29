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
    use FUSS_IO_Wall,              only: Initialize_Wall_File, IOwall
    use FUSS_Lib_Ghost,            only: Fill_Ghost_Cell, Fill_matIDg
    use FUSS_Mod_MPI,              only: mpi_is_root, partition_blocks, mpi_stop_all, mpi_stop_if_any, mpi_bcast_logical
    use FUSS_Mod_Timers,           only: timer_run_begin
    use FUSS_Mod_GhostExchange,    only: build_ghost_schedule, build_local_bc_index, &
                                         allocate_exchange_schedules
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
    setup_error = .false.
    if (mpi_is_root) call Check_Input ( setup_error )
    call Stop_Together ( setup_error, 'input file error (see above)' )
    call Check_FAS ()

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
    setup_error = .false.
    if (mpi_is_root) call Check_BC ( setup_error )
    call Stop_Together ( setup_error, 'boundary condition error (see above)' )
    
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

    ! Initialize the wall files, one per multigrid level (wall, wall-level<n>)
    allocate ( IOwall(obj_multigrid%MGL) )
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

        ! Interpolate initial solution on coarser domains
        if ( obj_multigrid%MGL > 1 ) then
          do m = 2, obj_multigrid%MGL
            call Restriction ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
            call Fill_Ghost_Cell ( simulation%domain(m) )
          end do
        endif

    end select

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
      use FUSS_Read_Ini, only: input_warnings
      implicit none
      logical, intent(out) :: has_error
      character(len=hlen) :: out
      integer :: r

      out = Validate_Registry()
      ! decomposition-map names a file of INPUT/: the folder, and the phase prefix of bc.txt, are put before it
      if ( index(out,'ERROR') == 0 .and. trim(obj_io%decomposition_map) /= 'none' .and. &
           index(obj_io%decomposition_map, '/') > 0 ) &
        out = '[ERROR] input.ini, [FUSS-IO] decomposition-map = '//trim(obj_io%decomposition_map)// &
              ': the value is the name of a file in INPUT/, without a folder'
      has_error = index(out,'ERROR') > 0
      if ( has_error ) then
        write(*,'(A,T35,A)') '   Input file', 'FAIL'
        write(*,'(4X,A)') trim(out)
      else
        write(*,'(A,T35,A)') '   Input file', 'OK'
      endif
      if (allocated(input_warnings)) then
        do r = 1, size(input_warnings)
          write(*,'(4X,A)') 'WARNING: '//trim(input_warnings(r))
        enddo
      endif

    end subroutine Check_Input


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