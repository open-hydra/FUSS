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
    use FUSS_IO_Wall,              only: Initialize_Wall_File
    use FUSS_Lib_Ghost,            only: Fill_Ghost_Cell, Fill_matIDg
    use FUSS_Mod_MPI,              only: mpi_is_root, partition_blocks
    use FUSS_Mod_Timers,           only: timer_run_begin
    use FUSS_Mod_GhostExchange,    only: build_ghost_schedule, build_local_bc_index
    implicit none
    type(FUSS_simulation_type), intent(inout) :: simulation
    ! Local
    integer :: m, ios

    !! ------------------------------------------------------
    !! ------------------------------------------------------
    ! Print Header
    if (mpi_is_root) call Print_Header ()
    
    ! Load tables (thermo, transport)
    call Load_Table ()
    if (mpi_is_root) call Check_Table ()

    ! Read input.ini
    call Read_Inifile ()
    ! Flags DERIVED from the ini must exist before Check_Input, or a guard that
    ! tests them is a comment for the input that sets them indirectly.
    ! Example: `irs-beta > 0` enables IRS, but the promotion lived in
    ! Assign_Setup, which runs after Check_Input, so `irs-beta = 0.5` with no
    ! `irs` key passed the mesh-motion IRS guard and smoothed the remapped
    ! residual anyway. Guard demo: run_guard_demos.sh irs_beta.
    if (obj_irs%beta > 0d0) obj_irs%enabled = .true.
    if (mpi_is_root) call Check_Input ()

    ! Allocate container objects Domain and IOfield.
    allocate ( simulation%domain  ( obj_multigrid%MGL ) )
    allocate ( simulation%IOfield ( obj_multigrid%MGL ) )
    simulation%domain(1)%itermax = obj_sim_param%iter_threshold

    ! Assign setup
    call Assign_Setup ()

    ! Read file for initial solution.
    call Read_IC ( simulation%IOfield(1) )
    if (mpi_is_root) call Check_IC ()

    ! Allocation of data structures and copy solution from IOfield.
    call Setup_Data_Structure ( simulation%domain(1), simulation%IOfield(1) )

    ! With multigrid, allocate Grid and IOfield 2,...,MGL
    if ( obj_multigrid%MGL > 1 ) then
      call Setup_Multigrid ( simulation )
    end if

    ! Read Boundary Conditions file
    call Setup_BC ( simulation%domain )
    if (mpi_is_root) call Check_BC ()
    
    ! MPI: partition blocks and build ghost cell communication schedule
    do m = 1, obj_multigrid%MGL
      call partition_blocks(simulation%domain(m)%nb, &
        [(product(simulation%domain(m)%blk(ios)%dim(1:3)), ios=1, simulation%domain(m)%nb)])
      call build_ghost_schedule(simulation%domain(m))
      call build_local_bc_index(simulation%domain(m))
      call fill_matIDg (simulation%domain(m))
    end do

    ! Read grid file and setup metrics in each grid level.
    do m = 1, obj_multigrid%MGL
      call Setup_Metrics ( simulation%domain(m) )
    end do
    
    ! Free heavy arrays on non-local blocks to save memory
    do m = 1, obj_multigrid%MGL
      call deallocate_remote_computation_data(simulation%domain(m))
    end do

    ! Setup probes location, output files, and variables to be printed.
    call Setup_Probes ( simulation%domain(1), obj_sim_param%newrun )
 
    ! Solution setup
    call Setup_Output_Solution ( simulation%IOfield )

    ! Initialize Wall file
    if (mpi_is_root .and. obj_io % write_wall) &
      call Initialize_Wall_File ( simulation%domain(1), simulation%IOfield(1)%tec%extension)

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
            if (obj_time_scheme%time_accurate) then
              simulation%domain(m)%time = simulation%IOfield(m)%solutiontime
              obj_sim_param%iter_general = 0
            else
              obj_sim_param%iter_general = int(simulation%IOfield(m)%solutiontime)
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

    ! If errors are found, print error messages and stop the simulation.
    if (mpi_is_root) call Stop_Simulation()

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


    subroutine Check_Table()
      implicit none
      logical :: has_error

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

      if (has_error) error stop 1

    end subroutine Check_Table


    subroutine Check_IC()
      implicit none

      if (index(obj_io%error_message,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Initial conditions', 'FAIL'
        write(*,'(4X,A)') trim(obj_io%error_message)
        error stop 1
      else
        write(*,'(A,T35,A)') '   Initial conditions', 'OK'
      endif

    end subroutine Check_IC


    subroutine Check_Input()
      use FUSS_Input_Registry
      use FUSS_Read_Ini, only: unknown_options_message
      implicit none
      character(len=hlen) :: out

      ! An option in a FUSS-owned section that nothing registered is not being
      ! applied, whatever the user believed. Reported BEFORE the registry
      ! validation, because a misspelt key makes the validation pass -- the
      ! parameter simply keeps its default and looks perfectly valid.
      if ( len_trim(unknown_options_message) > 0 ) then
        write(*,'(A,T35,A)') '   Input file', 'FAIL'
        write(*,'(4X,A)') trim(unknown_options_message)
        write(*,'(4X,A)') 'This part of input.ini is not being applied. A misspelt key keeps'
        write(*,'(4X,A)') 'its default, and a repeated section header is ignored entirely --'
        write(*,'(4X,A)') 'in both cases the run would otherwise look completely normal.'
        error stop 1
      endif

      out = Validate_Registry()
      if (index(out,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Input file', 'FAIL'
        write(*,'(4X,A)') trim(out)
        error stop 1
      else
        write(*,'(A,T35,A)') '   Input file', 'OK'
      endif

      call Check_Mesh_Motion_Compatibility()

    end subroutine Check_Input


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
    subroutine Check_Mesh_Motion_Compatibility()
      implicit none
      logical :: bad

      ! Test `law` directly rather than obj_mesh_motion%enabled: Check_Input
      ! runs immediately after Read_Inifile, whereas `enabled` is derived later
      ! in Assign_Setup. Reading the flag here would silently test an unset
      ! value and the guard would never fire -- which is exactly what happened
      ! the first time this was written.
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

      if ( bad ) error stop 1

    end subroutine Check_Mesh_Motion_Compatibility


    subroutine Check_BC()
      implicit none

      if (index(obj_io_bc%error_message,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Boundary conditions', 'FAIL'
        write(*,'(4X,A)') trim(obj_io_bc%error_message)
        error stop 1
      else
        write(*,'(A,T35,A)') '   Boundary conditions', 'OK'
      endif
      write(*,'(A)') ' ========================================================================================='

    end subroutine Check_BC


    subroutine Print_Warnings
      implicit none

      ! IO warnings
      if (index(obj_io%warning_message,'WARNING')>0) write(*,'(A)') obj_io%warning_message
      ! Numerical scheme warnings
      if (index(obj_time_scheme%warning_message,'WARNING')>0) write(*,'(A)') obj_time_scheme%warning_message
      if (index(obj_irs%warning_message,'WARNING')>0) write(*,'(A)') obj_irs%warning_message

    end subroutine Print_Warnings


    subroutine Stop_Simulation()
      implicit none
      logical :: has_error

      has_error = .false.

      ! Numerical scheme errors
      if (index(obj_time_scheme%error_message,'ERROR')>0) then
        write(*,'(A)') obj_time_scheme%error_message;  has_error = .true.
      endif
      if (index(obj_irs%error_message,'ERROR')>0) then
        write(*,'(A)') obj_irs%error_message;          has_error = .true.
      endif

      if (has_error) error stop 1

    end subroutine Stop_Simulation


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
      write(*,'(A,T35,A)') ' Time'
      write(*,'(A,T35,A)') '   Scheme', trim(obj_time_scheme%description)
      write(*,'(A,T35,A)') '   Integration variables', trim(obj_time_scheme%integration_variables)
      if (len_trim(obj_irs%description)>0) &
        write(*,'(A,T35,A)') '   Implicit residual smoothing', trim(obj_irs%description)
      write(*,'(A)') ' ========================================================================================='
      write(*,*)

    end subroutine Print_Shell_Info

  end subroutine FUSS_setup

end module FUSS_Wrap_Setup