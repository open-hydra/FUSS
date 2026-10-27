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
    use FUSS_IO_Probes,            only: Setup_Probes
    use FUSS_IO_Wall,              only: Initialize_Wall_File, IOwall
    use FUSS_Lib_Ghost,            only: Fill_Ghost_Cell, Fill_matIDg
    use FUSS_Mod_MPI,              only: mpi_is_root, partition_blocks, mpi_stop_all, mpi_bcast_logical
    use FUSS_Mod_Timers,           only: timer_run_begin
    use FUSS_Mod_GhostExchange,    only: build_ghost_schedule, build_local_bc_index, &
                                         allocate_exchange_schedules
    implicit none
    type(FUSS_simulation_type), intent(inout) :: simulation
    ! Local
    integer :: m, ios
    logical :: setup_error   ! an error found by a check of rank 0 (Stop_Together)

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

    ! Allocation of data structures and copy solution from IOfield.
    call Setup_Data_Structure ( simulation%domain(1), simulation%IOfield(1) )
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
    call partition_blocks(simulation%domain(1)%nb, &
      [(product(simulation%domain(1)%blk(ios)%dim(1:3)), ios=1, simulation%domain(1)%nb)])

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
            ! The restart file is the finest level's (IOfield(1)): the coarse IOfield carry no time
            if (obj_time_scheme%time_accurate) then
              simulation%domain(m)%time = simulation%IOfield(1)%solutiontime
              obj_sim_param%iter_general = 0
            else
              ! ORION a978930 reads SOLUTIONTIME with its sign: a negative time is a steady solution stamped with
              ! minus its iteration count (Wrap_Postprocess.f90:65); a positive one (a file written before the
              ! repin, |count| with STRANDID = 0) is read as before. Stopgap until the restart counters replace it.
              if (simulation%IOfield(1)%solutiontime < 0.0d0) then
                obj_sim_param%iter_general = nint(-simulation%IOfield(1)%solutiontime)
              else
                obj_sim_param%iter_general = int(simulation%IOfield(1)%solutiontime)
              endif
              simulation%domain(m)%time = -1.d0
            endif
          endif
          obj_sim_param%time_from_call = simulation%domain(m)%time 
          obj_sim_param%iter_from_call = 0
        enddo
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
      implicit none
      ! Local
      character(llen) :: eosword
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
      if ( obj_multigrid%MGL > 1 ) &
        write(*,'(A,T35,I0)') '   Multigrid levels', obj_multigrid%MGL

      ! ----- Boundary conditions -----
      call Print_BC_Summary ()

      ! ----- IO -----
      write(*,*)
      write(*,'(A)') ' Input/Output'
      write(*,'(A,T35,A)') '   Initial conditions file', trim(obj_io%nameinit)
      write(*,'(A,T35,A)') '   Solution format', trim(obj_io%sol_format)
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
      write(*,'(A,T35,A)') '   Diffusive face metric', 'inverse-mean'
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