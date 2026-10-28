module FUSS_Wrap_Postprocess
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: FUSS_postprocess

  logical, private :: id_stampa_set = .false.   ! numbering of the outputs settled (restored by a restart, or scanned)
  logical, private :: irs_runtime_warned = .false.
  ! FAS setup state: the runtime re-read may not switch cycles on (coarse-level arrays are
  ! allocated at setup only) nor leave a cycling root with no smoothing iterations
  logical, private :: fas_setup_saved = .false.
  logical, private :: fas_cycle_warned = .false., fas_iters_warned = .false.
  character(len=8), private :: fas_cycle_setup = 'none'
  integer, private :: fas_pre_setup = 2, fas_post_setup = 1, fas_coarse_setup = 10
  character(81), private :: RANS_shell_format  =    "('FUSS | Iter =', i9, ' | Global iter =', i9, ' | Energy  residual =', E13.6)"
  character(69), private :: URANS_shell_format =    "('FUSS | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"
  character(96), private :: RANS_shell_format_MG  = "('FUSS Grid Level', i2, ' | Iter =', i9, ' | Global iter =', i9, ' | Energy  residual =', E13.6)"
  character(84), private :: URANS_shell_format_MG = "('FUSS Grid Level', i2, ' | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"

contains

  subroutine FUSS_postprocess ( simulation )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m
    use FUSS_Global_m
    use FUSS_Parameters_m
    use FUSS_IO_Wall,        only: IOwall, Write_Wall_Solution
    use FUSS_IO_Solution,    only: Write_Solution
    use FUSS_Mod_Diagnostic, only: Write_Diagnostic
    use FUSS_IO_Probes,      only: Write_Probes_Data, Close_Probes
    use FUSS_Mod_Multigrid,  only: Prolongation
    use FUSS_Read_Ini,       only: Read_Inifile_Runtime
    use FUSS_Mod_MPI,        only: mpi_is_root, mpi_reduce_sum_r8
    use FUSS_Lib_Newstate,   only: fas_update_fallbacks
    use FUSS_Lib_Multigrid,  only: fas_fallbacks
    use FUSS_Mod_Timers,     only: timer_summary
    use IR_precision
    implicit none
    type(FUSS_simulation_type), intent(inout) :: simulation
    ! Local
    character(len=llen) :: solfile, dgsfile, wallfile, mgsol, mgwall, numsol, numwall, numdgs
    logical  :: iter_based
    integer  :: m, level
    real(R8) :: nfall

    ! FAS keys as validated at setup, before the first runtime re-read can change them
    if (.not. fas_setup_saved) then
      fas_setup_saved  = .true.
      fas_cycle_setup  = obj_multigrid%cycle
      fas_pre_setup    = obj_multigrid%pre_iter
      fas_post_setup   = obj_multigrid%post_iter
      fas_coarse_setup = obj_multigrid%coarse_iter
    endif

    level = obj_multigrid%MG_level

    ! Solution time of level 1 (also the prolongated files) and of the level being computed
    if ( obj_time_scheme%time_accurate ) then
      simulation%IOfield(1)%solutiontime = simulation%domain(1)%time
      IOwall(1)%solutiontime = simulation%domain(1)%time
      simulation%IOfield(level)%solutiontime = simulation%domain(level)%time
      IOwall(level)%solutiontime = simulation%domain(level)%time
    else
      simulation%IOfield(1)%solutiontime = -real(obj_sim_param%iter_general,8)
      IOwall(1)%solutiontime = -real(obj_sim_param%iter_general,8)
      simulation%IOfield(level)%solutiontime = -real(obj_sim_param%iter_general,8)
      IOwall(level)%solutiontime = -real(obj_sim_param%iter_general,8)
    endif

    if (level == 1) then
      solfile  = '/'//trim(FUSS_phase_prefix)//'field'
      wallfile = '/'//trim(FUSS_phase_prefix)//'wall'
      dgsfile  = '/'//trim(FUSS_phase_prefix)//'diagnostic'
    else
      solfile  = '/'//trim(FUSS_phase_prefix)//'field-level'//trim(str(.true.,level))
      wallfile = '/'//trim(FUSS_phase_prefix)//'wall-level'//trim(str(.true.,level))
      dgsfile  = '/'//trim(FUSS_phase_prefix)//'diagnostic-level'//trim(str(.true.,level))
      mgsol    = '/'//trim(FUSS_phase_prefix)//'field-prolongated'
      mgwall   = '/'//trim(FUSS_phase_prefix)//'wall-prolongated'
    endif

    ! END OF SIMULATION
    if ( obj_sim_param%TODO == 3 ) then

      if (mpi_is_root) then
        ! Residual history
        if ( mod (simulation%domain(1) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(1)%time, obj_sim_param%residuotot
          flush(obj_io%unitRES)
        endif

        ! Shell
        if ( obj_time_scheme%time_accurate ) then
          if (obj_multigrid%MGL > 1 ) then
            write (*,URANS_shell_format_MG) level, obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
          else
            write (*,URANS_shell_format) obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
          endif
        else
          if (obj_multigrid%MGL > 1 ) then
            write (*,RANS_shell_format_MG) level, obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot
          else
            write (*,RANS_shell_format) obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot
          endif
        endif
      end if

      ! Last probe sample (collective: the owner rank writes)
      if (level == 1) call Write_Probes_Data( obj_sim_param%iter_general, simulation%domain(1)%time )

      ! The solution due at the last iteration (sol-diter, else sol-dtime), as at any other output iteration (TODO == 2): it
      ! takes the next output number and, when due by time only, sets the time of the last output, which the final solution
      ! below then holds for a restart. With sol-overwrite = true the final solution is that output, written once.
      ! Collective: every rank evaluates the same condition. A run ends on the finest level only (Explicit_Step sets
      ! TODO = 3 at level 1)
      iter_based = ( mod (simulation%domain(1) % iter, obj_io%sol_diter) == 0 )
      if ( level == 1 .and. ( iter_based .or. &
           simulation%domain(1) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime ) ) then
        if (.not. id_stampa_set) call Settle_Numbering ()
        obj_io%output_number = obj_io%output_number + 1
        if ( .not. iter_based ) obj_sim_param%time_from_call = simulation%domain(1) % time
        if ( .not. obj_io%sol_overwrite ) then
          if (mpi_is_root) then
            if ( iter_based ) then
              write(*,*) ' ... writing iter-based solution'
            else
              write(*,*) ' ... writing time-based solution'
            endif
          endif
          numsol  = trim(solfile)//trim(str(.true.,obj_io%output_number))
          numwall = trim(wallfile)//trim(str(.true.,obj_io%output_number))
          numdgs  = trim(dgsfile)//trim(str(.true.,obj_io%output_number))
          call Write_Solution ( simulation%domain(1), simulation%IOfield(1), numsol )
          if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(1), numwall )
          if (.not. obj_time_scheme%time_accurate) &
            call Write_Diagnostic ( simulation%domain(1), simulation%IOfield(1), numdgs )
        endif
      endif

      ! Timing summary (collective — every rank contributes its accumulators)
      call timer_summary()

      ! Write the level being computed (collective)
      call Write_Solution ( simulation%domain(level), simulation%IOfield(level), solfile )
      if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(level), wallfile )
      if (level == 1 .and. .not. obj_time_scheme%time_accurate) &
        call Write_Diagnostic ( simulation%domain(1), simulation%IOfield(1), dgsfile )

      ! The history files hold their last rows: close them before the ranks finalize MPI (a unit left open is closed by
      ! the Fortran run-time library only when the process ends, after MPI_FINALIZE)
      call Close_Probes()
      if (mpi_is_root) close(obj_io%unitRES)
    

    ! INTERMEDIATE SOLUTION EVALUATION
    elseif ( obj_sim_param%TODO == 2 ) then

      if (level == 1) then
        ! Restart with numbered snapshots: the number goes on from the solution read (obj_io%output_number, restored with
        ! the counters). A solution without counters (written before them) gives no number: the legacy path scans
        ! OUTPUT/ for the last field<n> and continues after it
        if (.not. id_stampa_set) call Settle_Numbering ()
        obj_io%output_number = obj_io%output_number + 1
        if (.not.obj_io%sol_overwrite) solfile  = trim(solfile)//trim(str(.true.,obj_io%output_number))
        if (.not.obj_io%sol_overwrite) wallfile = trim(wallfile)//trim(str(.true.,obj_io%output_number))
        if (.not.obj_io%sol_overwrite) dgsfile  = trim(dgsfile)//trim(str(.true.,obj_io%output_number))
      endif

      ! Residuals history
      if (mpi_is_root) then
        if ( mod(simulation%domain(level) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(level)%time, obj_sim_param%residuotot
          flush(obj_io%unitRES)
        endif
      end if

      ! Probes
      if (level == 1) call Write_Probes_Data( obj_sim_param%iter_general, simulation%domain(1)%time)

      ! Shell (root only)
      if (mpi_is_root) then
        if ( mod (simulation%domain(level) % iter, obj_io%shell_diter) == 0d0 ) then
          if (obj_time_scheme%time_accurate) then
            if (obj_multigrid%MGL > 1 ) then
              write (*,URANS_shell_format_MG) level, obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            else
              write (*,URANS_shell_format) obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            endif
          else
            if (obj_multigrid%MGL > 1 ) then
              write (*,RANS_shell_format_MG) level, obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot
            else
              write (*,RANS_shell_format) obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot
            endif
          endif
        endif
      end if
      
      ! Solution output (collective — all ranks participate in gather)
      ! Iter-based solution
      if ( mod(simulation%domain(level) % iter, obj_io%sol_diter) == 0d0 ) then
        if (mpi_is_root) write(*,*) ' ... writing iter-based solution'
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
        if (obj_multigrid%MGL > 1 .and. level > 1 .and. prolongated_due()) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(1), mgwall )
          obj_multigrid%prolongated_last = simulation%domain(level) % iter
        endif
        if (level == 1 .and. .not. obj_time_scheme%time_accurate) &
          call Write_Diagnostic ( simulation%domain(level), simulation%IOfield(level), dgsfile )
      ! Time-based solution
      elseif ( simulation%domain(level) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime ) then
        obj_sim_param%time_from_call = simulation%domain(level) % time
        if (mpi_is_root) write(*,*) ' ... writing time-based solution'
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
        if (obj_multigrid%MGL > 1 .and. level > 1 .and. prolongated_due()) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(1), mgwall )
          obj_multigrid%prolongated_last = simulation%domain(level) % iter
        endif
        if (level == 1 .and. .not. obj_time_scheme%time_accurate) &
          call Write_Diagnostic ( simulation%domain(level), simulation%IOfield(level), dgsfile )
      endif

      if (obj_multigrid%change_MG) then
        if (mpi_is_root) write(*,*) ' ... writing solution of grid-level ', level
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
        if (obj_multigrid%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(1), mgwall )
        endif
      endif

      if (obj_multigrid%change_MG) then
        obj_multigrid%MG_level = obj_multigrid%MG_level - 1
        obj_multigrid%prolongated_last = 0
        obj_sim_param%iter_from_call = 0
        obj_sim_param%time_from_call = 0.0
      endif


    ! AUXILIARY OPERATIONS AT EACH ITERATION
    elseif ( obj_sim_param%TODO == 1 ) then
      
      if (mpi_is_root) then
        if ( mod (simulation%domain(level) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(level)%time, obj_sim_param%residuotot
          flush(obj_io%unitRES)
        end if
      end if

      ! Probes
      if (level == 1) call Write_Probes_Data( obj_sim_param%iter_general, simulation%domain(1)%time )
      
      ! Shell
      if (mpi_is_root) then
        if ( mod (simulation%domain(level) % iter, obj_io%shell_diter) == 0d0 ) then
          if (obj_time_scheme%time_accurate) then
            if (obj_multigrid%MGL > 1 ) then
              write (*,URANS_shell_format_MG) level, obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            else
              write (*,URANS_shell_format) obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            endif
          else
            if (obj_multigrid%MGL > 1 ) then
              write (*,RANS_shell_format_MG) level, obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot
            else
              write (*,RANS_shell_format) obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot
            endif
          endif
        endif
      end if

    endif

    ! FAS corrections cut by the admissibility guard (Lib_Multigrid) since the last report, summed
    ! over the ranks (the IRS fallback report went with the fallback, October 2026)
    if ( trim(obj_multigrid%cycle) /= 'none' ) then
      if ( mod(simulation%domain(level) % iter, obj_io%shell_diter) == 0d0 ) then
        nfall = real(fas_fallbacks, R8)
        call mpi_reduce_sum_r8( nfall )
        if ( mpi_is_root .and. nfall > 0.5d0 ) write(*,'(A,I0,A)') ' FAS: ', nint(nfall), &
          ' prolongated corrections halved or dropped since the last report'
        fas_fallbacks = 0
        ! Coarse-level updates halved because the forcing left the table (Lib_Newstate)
        nfall = real(fas_update_fallbacks, R8)
        call mpi_reduce_sum_r8( nfall )
        if ( mpi_is_root .and. nfall > 0.5d0 ) write(*,'(A,I0,A)') ' FAS: ', nint(nfall), &
          ' coarse smoothing updates halved since the last report'
        fas_update_fallbacks = 0
      endif
    endif

    ! Update input data
    if ( mod (simulation%domain(level) % iter, obj_io%ini_diter) == 0d0 ) then
      call Read_Inifile_Runtime( simulation%domain(level)%iter )
      ! level<n>-iter / iter-threshold take effect at runtime
      if ( obj_multigrid%MGL > 1 ) then
        do m = 1, obj_multigrid%MGL
          simulation%domain(m)%itermax = obj_multigrid%iter_threshold(m)
        enddo
      else
        simulation%domain(1)%itermax = obj_sim_param%iter_threshold
      endif
      ! IRS follows irs-beta, but its arrays exist only if IRS was on at setup
      obj_irs%enabled = ( obj_irs%beta > 0d0 )
      if ( obj_irs%enabled .and. .not. obj_irs%at_setup ) then
        if ( mpi_is_root .and. .not. irs_runtime_warned ) write(*,'(A)') &
          ' WARNING: IRS cannot be switched on at runtime (start with irs-beta > 0): irs-beta ignored'
        irs_runtime_warned = .true.
        obj_irs%enabled = .false.
      endif
      ! FAS coarse-level arrays exist only if a cycle was on at setup (same pattern as IRS)
      if ( trim(fas_cycle_setup) == 'none' .and. trim(obj_multigrid%cycle) /= 'none' ) then
        if ( mpi_is_root .and. .not. fas_cycle_warned ) write(*,'(A)') &
          ' WARNING: [FUSS-Multigrid] cycle cannot be switched on at runtime (start with cycle = v or w): kept none'
        fas_cycle_warned = .true.
        obj_multigrid%cycle = fas_cycle_setup
      endif
      ! A cycling run keeps the setup guards: steady, cons, and a root that iterates
      if ( trim(obj_multigrid%cycle) /= 'none' ) then
        if ( obj_multigrid%pre_iter + obj_multigrid%post_iter < 1 .or. obj_multigrid%coarse_iter < 1 .or. &
             obj_time_scheme%time_accurate .or. trim(obj_time_scheme%integration_variables) == 'prim' ) then
          if ( mpi_is_root .and. .not. fas_iters_warned ) write(*,'(A)') &
            ' WARNING: runtime edit breaks the FAS guards (pre-iter + post-iter >= 1, coarse-iter >= 1,'// &
            ' steady, cons): setup values kept'
          fas_iters_warned = .true.
          obj_multigrid%pre_iter    = fas_pre_setup
          obj_multigrid%post_iter   = fas_post_setup
          obj_multigrid%coarse_iter = fas_coarse_setup
          obj_time_scheme%time_accurate = .false.
          obj_time_scheme%integration_variables = 'cons'
        endif
      endif
    end if
    
    
    ! PRINT FOR HYDRA-COUPLING PROBLEMS
    if ( obj_sim_param%HYDRA_postprocess ) then
      if (mpi_is_root) write(*,*) ' ... writing HYDRA-coupling solution'
      call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
      if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
      if (obj_multigrid%MGL > 1 .and. level > 1) then
        if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
        do m = obj_multigrid%MG_level, 2, -1
          call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
        enddo
        call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
        if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(1), mgwall )
      endif
    endif

  contains

    ! The numbering of the outputs of a restart that gets no counter from its solution: scan OUTPUT/ (legacy path)
    subroutine Settle_Numbering ()
      id_stampa_set = .true.
      if (.not. obj_sim_param%newrun .and. .not. obj_io%sol_overwrite .and. .not. obj_io%restart%resumed) &
        call Resume_File_Counter ( trim(simulation%IOfield(1)%tec%extension) )
    end subroutine Settle_Numbering

    ! Periodic prolongated writes: every field write, or prolongated-diter level iterations apart
    logical function prolongated_due ()
      prolongated_due = obj_multigrid%prolongated_diter == 0 .or. &
        simulation%domain(level) % iter - obj_multigrid%prolongated_last >= obj_multigrid%prolongated_diter
    end function prolongated_due

  end subroutine FUSS_postprocess

  !> Legacy path, for a restart from a solution that holds no counters (written before them): the highest n with
  !> OUTPUT/field<n><ext> present, so that a restart with numbered snapshots appends instead of overwriting. A solution with
  !> counters gives the number itself (obj_io%output_number).
  subroutine Resume_File_Counter ( ext )
    use FUSS_Global_m, only: FUSS_phase_prefix
    use FUSS_Config_Types_m, only: obj_io
    use FUSS_Mod_MPI, only: mpi_is_root
    use ir_precision, only: str
    implicit none
    character(len=*), intent(in) :: ext
    ! Local
    logical :: present

    do
      inquire ( file='OUTPUT/'//trim(FUSS_phase_prefix)//'field'//trim(str(.true.,obj_io%output_number+1))//trim(ext), &
                exist=present )
      if (.not. present) exit
      obj_io%output_number = obj_io%output_number + 1
    enddo
    if (obj_io%output_number > 0 .and. mpi_is_root) &
      write(*,'(A)') ' Numbered output continues from '//trim(str(.true.,obj_io%output_number+1))

  end subroutine Resume_File_Counter

end module FUSS_Wrap_Postprocess
