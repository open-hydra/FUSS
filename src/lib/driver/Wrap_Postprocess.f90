module FUSS_Wrap_Postprocess
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: FUSS_postprocess

  integer, private :: id_stampa=0
  logical, private :: id_stampa_set = .false.   ! counter resumed from OUTPUT/ on a restart
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
    use FUSS_IO_Probes,      only: Write_Probes_Data
    use FUSS_Mod_Multigrid,  only: Prolongation
    use FUSS_Read_Ini,       only: Read_Inifile_Runtime
    use FUSS_Mod_MPI,        only: mpi_is_root, mpi_reduce_sum_r8
    use FUSS_Lib_Newstate,   only: irs_fallbacks, fas_update_fallbacks
    use FUSS_Lib_Multigrid,  only: fas_fallbacks
    use FUSS_Mod_Timers,     only: timer_summary
    use IR_precision
    implicit none
    type(FUSS_simulation_type), intent(inout) :: simulation
    ! Local
    character(len=llen) :: solfile, dgsfile, wallfile, mgsol, mgwall
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

      ! Timing summary (collective — every rank contributes its accumulators)
      call timer_summary()

      ! Write the level being computed (collective)
      call Write_Solution ( simulation%domain(level), simulation%IOfield(level), solfile )
      if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(level), wallfile )
      if (level == 1 .and. .not. obj_time_scheme%time_accurate) &
        call Write_Diagnostic ( simulation%domain(1), simulation%IOfield(1), dgsfile )
    

    ! INTERMEDIATE SOLUTION EVALUATION
    elseif ( obj_sim_param%TODO == 2 ) then

      if (level == 1) then
        ! Restart with numbered snapshots: continue after the last field<n> in OUTPUT/
        if (.not. id_stampa_set) then
          id_stampa_set = .true.
          if (.not. obj_sim_param%newrun .and. .not. obj_io%sol_overwrite) &
            call Resume_File_Counter ( trim(simulation%IOfield(1)%tec%extension) )
        endif
        id_stampa = id_stampa + 1
        if (.not.obj_io%sol_overwrite) solfile  = trim(solfile)//trim(str(.true.,id_stampa))
        if (.not.obj_io%sol_overwrite) wallfile = trim(wallfile)//trim(str(.true.,id_stampa))
        if (.not.obj_io%sol_overwrite) dgsfile  = trim(dgsfile)//trim(str(.true.,id_stampa))
      endif

      ! Residuals history
      if (mpi_is_root) then
        if ( mod(simulation%domain(level) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(level)%time, obj_sim_param%residuotot
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

    ! IRS fallbacks (Lib_Newstate) since the last report, summed over the ranks
    if ( obj_irs%at_setup .and. mod(simulation%domain(level) % iter, obj_io%shell_diter) == 0d0 ) then
      nfall = real(irs_fallbacks, R8)
      call mpi_reduce_sum_r8( nfall )
      if ( mpi_is_root .and. nfall > 0.5d0 ) write(*,'(A,I0,A,I0)') ' IRS: fallback to the own update in ', &
        nint(nfall), ' cell updates since the last report, level ', level
      irs_fallbacks = 0
    endif

    ! FAS corrections cut by the admissibility guard (Lib_Multigrid), same cadence
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
      call Read_Inifile_Runtime()
      ! level<n>-iter / iter-threshold take effect at runtime
      if ( obj_multigrid%MGL > 1 ) then
        do m = 1, obj_multigrid%MGL
          simulation%domain(m)%itermax = obj_multigrid%iter_threshold(m)
        enddo
      else
        simulation%domain(1)%itermax = obj_sim_param%iter_threshold
      endif
      ! IRS arrays exist only if IRS was on at setup
      if ( obj_irs%enabled .and. .not. obj_irs%at_setup ) then
        if ( mpi_is_root .and. .not. irs_runtime_warned ) write(*,'(A)') &
          ' WARNING: IRS cannot be switched on at runtime (start with irs-beta > 0): irs ignored'
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

    ! Periodic prolongated writes: every field write, or prolongated-diter level iterations apart
    logical function prolongated_due ()
      prolongated_due = obj_multigrid%prolongated_diter == 0 .or. &
        simulation%domain(level) % iter - obj_multigrid%prolongated_last >= obj_multigrid%prolongated_diter
    end function prolongated_due

  end subroutine FUSS_postprocess

  !> Highest n with OUTPUT/field<n><ext> present, so a restart appends instead of overwriting
  subroutine Resume_File_Counter ( ext )
    use FUSS_Global_m, only: FUSS_phase_prefix
    use ir_precision, only: str
    implicit none
    character(len=*), intent(in) :: ext
    ! Local
    logical :: present

    do
      inquire ( file='OUTPUT/'//trim(FUSS_phase_prefix)//'field'//trim(str(.true.,id_stampa+1))//trim(ext), &
                exist=present )
      if (.not. present) exit
      id_stampa = id_stampa + 1
    enddo
    if (id_stampa > 0) write(*,'(A)') ' Numbered output continues from '//trim(str(.true.,id_stampa+1))

  end subroutine Resume_File_Counter

end module FUSS_Wrap_Postprocess