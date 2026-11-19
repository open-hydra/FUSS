module FUSS_Read_Ini
  
  implicit none

  ! Warnings on keys that are set but have no effect, printed by Check_Input
  character(len=256), allocatable :: input_warnings(:)
  !> Set by Read_Inifile when input.ini carries options nothing reads, or a repeated section (plan 10 D1). Reported by
  !> Check_Input (Wrap_Setup), which refuses the input: every such option is listed in one pass.
  character(len=1024), public :: unknown_options_message = ''
  !> Set by Read_Inifile when a string value is longer than its variable (it would be cut without a message). Reported
  !> by Check_Input.
  character(len=1024), public :: long_values_message = ''

contains

  subroutine Read_Inifile()
    use Finer,               only: file_ini
    use FUSS_Read_Sim_Param, only: Register_Sim_Param
    use FUSS_Read_IO,        only: Register_IO_Fields, Register_Probes
    use FUSS_Read_Numerics,  only: Register_Numerics
    use FUSS_Read_Pyrolysis, only: Register_Pyrolysis
    use FUSS_Register_Darcy, only: Register_Darcy
    use FUSS_Register_Ablation, only: Register_Ablation
    use FUSS_Backend_INI,    only: Load_Ini, Scan_Ini, Snapshot_Ini_Text, Check_Unknown_Options, Check_Value_Lengths
    use FUSS_Input_Registry
    implicit none
    ! Local
    type(file_ini) :: fini
    integer :: nprobes, nmgl
    character(len=16), allocatable :: probes_name(:)


    ! Load input.ini: a file that cannot be read stops the run
    call Load_Input_File(fini, 'Input file')

    ! Scan input.ini for unknown number of probes and multigrid levels
    call Scan_Ini(fini, nprobes, probes_name, nmgl)

    ! Build registry entries
    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Probes(nprobes, probes_name)
    call Register_Numerics(nmgl)
    call Register_Pyrolysis()
    call Register_Darcy()
    call Register_Ablation()

    ! Registry is built, now load the values from the ini file
    call Load_Ini(fini)
    ! The text of every key at set-up, which each reading during the run is compared with (Reload_Runtime_Ini). Here and
    ! not in Load_Ini, which every reading during the run calls too (Read_Inifile_Runtime): there it would overwrite the
    ! set-up text before the comparison, and no change would ever be reported
    call Snapshot_Ini_Text(fini)

    ! The keys of input.ini that the reader does not give
    call Warn_Unread_Keys(fini)

    ! Keys set to no effect
    call Find_Input_Warnings(fini, nmgl)

    ! Anything in a FUSS-owned section that nothing registered is a typo or a stale key, and either way it is NOT
    ! being applied (plan 10 D1, plan 14 D7). Recorded here, reported by Check_Input.
    call Check_Unknown_Options(fini, unknown_options_message)
    call Check_Value_Lengths(fini, long_values_message)

  end subroutine Read_Inifile


  !> Keys of input.ini that are set but have no effect: iter-threshold when the run has more than
  !> one multigrid level (each level stops at its own level<m>-iter)
  subroutine Find_Input_Warnings(fini, nmgl)
    use Finer, only: file_ini
    implicit none
    type(file_ini), intent(in) :: fini
    integer, intent(in) :: nmgl
    ! Local
    character(len=256) :: sec, line

    if (allocated(input_warnings)) deallocate(input_warnings)
    allocate(input_warnings(0))
    sec = ''
    if (nmgl > 1 .and. fini%has_option(option_name='iter-threshold', section_name=sec)) then
      write(line,'(A,I0,A)') '['//trim(sec)//'] iter-threshold is ignored with ', nmgl, &
        ' multigrid levels: each level stops at its level<m>-iter'
      input_warnings = [ input_warnings, line ]
    endif

  end subroutine Find_Input_Warnings


  !> Prints a warning for every line "key = value" of input.ini whose key the INI reader (FiNeR) does not give in its
  !> section: that key is not read, and its parameter keeps its default. Comment lines (they start with !, ; or #) are
  !> skipped, as the reader skips them since FiNeR 2.1.0, also when they hold "=" or a bracket.
  subroutine Warn_Unread_Keys(fini)
    use Finer,        only: file_ini
    use FUSS_Mod_MPI, only: mpi_is_root
    implicit none
    type(file_ini), intent(in) :: fini
    ! Local
    character(len=1024) :: line, section, key, val
    integer :: unit, ios, n, eq, error

    if (.not. mpi_is_root) return
    open(newunit=unit, file='input.ini', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    section = ''
    n = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      n = n + 1
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (scan(line(1:1), '!;#') == 1) cycle
      if (line(1:1) == '[') then
        if (index(line, ']') > 2) section = adjustl(line(2:index(line, ']')-1))
        cycle
      end if
      eq = index(line, '=')
      if (eq <= 1) cycle
      key = adjustl(line(1:eq-1))
      call fini%get(section_name=trim(section), option_name=trim(key), val=val, error=error)
      if (error /= 0) write(*,'(A,I0,A)') '[WARNING] input.ini, line ', n, ': '//trim(key)//' of ['//trim(section)// &
        '] is not read (the INI reader does not give it): its parameter keeps its default'
    end do
    close(unit)

  end subroutine Warn_Unread_Keys


  !> input.ini read again during the run, every ini-diter iterations (FUSS_postprocess, at iteration iter).
  !> The reading (plan 14 D8-A2): rank 0 reads the file and passes its text to every rank (Load_Input_File); the whole
  !> reading keeps to the rules of the registry, as at set-up, and an error of it does what ini-error says (warn, the
  !> default: the reading is not taken, Write_Not_Taken; stop: every rank stops, exit code 1). The value of ini-error in
  !> force when the reading starts decides.
  !> What a reading without error changes (plan 10 D3, user decision 2026-09-14): only the keys of RUNTIME_MUTABLE
  !> (Backend_INI), each checked again on its own (Reload_Runtime_Ini); a changed value of any other key is IGNORED, the
  !> run goes on with the value it started with, so the set-up guards (mesh motion, pyrolysis, Darcy, ablation) cannot be
  !> bypassed mid-run. Every ignored change, and every reading not taken, is a line of
  !> OUTPUT/<prefix>runtime-ini-ignored.txt (rank 0), created only then: its existence is the signal.
  !> Every rank checks and applies the same text, so every rank takes the same decisions (no broadcast of values).
  subroutine Read_Inifile_Runtime(iter)
    use Finer,               only: file_ini
    use FUSS_Backend_INI,    only: Load_Ini, Reload_Runtime_Ini
    use FUSS_Input_Registry
    use FUSS_Config_Types_m, only: obj_io
    use FUSS_Mod_MPI,        only: mpi_rank_, mpi_is_root, mpi_allreduce_max_r8, mpi_allreduce_min_r8, mpi_stop_all
    use iso_fortran_env,     only: R8 => real64
    implicit none
    integer, intent(in) :: iter   !< The iteration at which input.ini is read again
    ! Local
    type(file_ini) :: fini
    character(len=1024) :: out
    character(len=4096) :: report
    real(R8) :: flag, flag_all, first
    logical  :: warn, taken

    ! The value of ini-error in force when the reading starts
    warn = ( trim(obj_io%ini_error) == 'warn' )

    ! Load input.ini: a file that cannot be read stops the run, or with warn leaves the reading not taken
    call Load_Input_File(fini, 'Input file, read again', warn, iter, taken)
    if (.not. taken) then
      if (mpi_is_root) call Append_Ignored(iter, 'reading not taken: input.ini cannot be read (see the log)')
      return
    end if

    ! The whole reading against the rules of the registry: loaded into the registry, checked, and the values before the
    ! reading put back whatever the verdict (Save_Values / Restore_Values); the keys it may change are applied below
    call Save_Values()
    call Load_Ini(fini)
    out = Validate_Registry()
    call Restore_Values()

    flag = 0.0_R8
    if (index(out,'ERROR') > 0) flag = 1.0_R8
    call mpi_allreduce_max_r8(flag, flag_all)
    if (flag_all > 0.0_R8) then
      call mpi_allreduce_min_r8(merge(real(mpi_rank_, R8), huge(1.0_R8), flag > 0.0_R8), first)
      if (warn) then
        if (flag > 0.0_R8 .and. nint(first) == mpi_rank_) call Write_Not_Taken('Input file, read again', iter, out)
        if (mpi_is_root) call Append_Ignored(iter, 'reading not taken: '//Without_Error_Tag(out))
        return
      end if
      if (flag > 0.0_R8 .and. nint(first) == mpi_rank_) then
        write(*,'(A,T35,A)') '   Input file, read again', 'FAIL'
        write(*,'(4X,A)') trim(out)
      end if
      if (flag > 0.0_R8) call mpi_stop_all('input file error, input.ini read again during the run (see above)')
      call mpi_stop_all('stopped with the ranks that found an error in input.ini read again')
    end if

    ! A reading without error: the runtime-mutable keys whose text changed are applied, every other change is ignored
    call Reload_Runtime_Ini(fini, iter, report)
    if (mpi_is_root .and. len_trim(report) > 0) call Append_Ignored(-1, report)

  end subroutine Read_Inifile_Runtime


  !> A record of OUTPUT/<prefix>runtime-ini-ignored.txt (rank 0 only), opened in append mode, so that the file exists only
  !> once something was ignored or not taken (plan 10 D3). iter >= 0: one line 'iteration <iter>:  <text>'; iter < 0:
  !> text already holds its lines (Reload_Runtime_Ini's report).
  subroutine Append_Ignored(iter, text)
    use FUSS_Global_m, only: FUSS_phase_prefix
    implicit none
    integer,          intent(in) :: iter
    character(len=*), intent(in) :: text
    ! Local
    integer :: u
    character(len=32) :: head

    open(newunit=u, file='OUTPUT/'//trim(FUSS_phase_prefix)//'runtime-ini-ignored.txt', &
         status='unknown', position='append', action='write')
    if (iter >= 0) then
      write(head,'(A,I0,A)') 'iteration ', iter, ':  '
      write(u,'(A)') trim(head)//'  '//trim(text)
    else
      write(u,'(A)') trim(text)
    end if
    close(u)

  end subroutine Append_Ignored


  !> An error of the registry without its leading '[ERROR] ' (the reading not taken is a warning, not an error)
  function Without_Error_Tag(text) result(r)
    implicit none
    character(len=*), intent(in) :: text
    character(len=:), allocatable :: r
    r = trim(adjustl(text))
    if (index(r, '[ERROR] ') == 1) r = r(9:)
  end function Without_Error_Tag


  !> Loads input.ini into fini, and stops the run if it cannot be read: a file that cannot be opened (missing, or without
  !> read permission), a file with no line that is not blank, a file whose read ends with an I/O error, or an error of
  !> the INI reader (FiNeR load gives an error code other than 0: a file with no section, a line that it cannot parse;
  !> FiNeR names such a line on stderr, without its number). Without this stop, the parameters would keep their defaults
  !> at set-up and the values of the previous reading during the run, and an empty file would end the run in the reader.
  !> During the run with ini-error = warn (warn present and true) such a file leaves the reading not taken instead: a
  !> warning (Write_Not_Taken) and taken false, so that the caller keeps every value it had (Read_Inifile_Runtime).
  !> Rank 0 alone reads the file, as FiNeR reads a file (string%read_file of StringiFor), and passes its text to every
  !> rank (mpi_bcast_chars); every rank, rank 0 too, gives that text to FiNeR (load from a source), so every rank takes
  !> the same values. If each rank read input.ini on its own, a file rewritten while the ranks read it during the run
  !> could give them different values, and a value that decides a collective call (sol-diter, res-diter, ini-diter)
  !> would make them call their collectives in different orders (forbidden by MPI-4.1, section 6.14), and hang. A file
  !> that rank 0 cannot read stops every rank, or leaves the reading not taken on every rank: rank 0 writes why and
  !> passes the outcome (mpi_bcast_integer), then every rank stops with mpi_stop_all (exit code 1), or returns; an error
  !> of the INI reader, the same on every rank, is decided together too (the largest of a flag), and the first rank that
  !> found it writes it.
  subroutine Load_Input_File(fini, label, warn, iter, taken)
    use Finer,           only: file_ini
    use stringifor,      only: string
    use FUSS_Mod_MPI,    only: mpi_is_root, mpi_rank_, mpi_bcast_integer, mpi_bcast_chars, mpi_allreduce_max_r8, &
                               mpi_allreduce_min_r8, mpi_stop_all
    use iso_fortran_env, only: R8 => real64
    implicit none
    type(file_ini),    intent(inout) :: fini
    character(len=*),  intent(in)    :: label   !< Line of the set-up summary: Input file, or Input file, read again
    logical, optional, intent(in)    :: warn    !< A reading during the run with ini-error = warn: an error, not taken
    integer, optional, intent(in)    :: iter    !< The iteration of a reading during the run, for the warning
    logical, optional, intent(out)   :: taken   !< The file is loaded (false: a reading that is not taken, with warn)
    ! Local
    character(len=512)  :: reason, msg
    character(len=1024) :: line
    character(len=:), allocatable :: text      !< The text of input.ini that rank 0 reads
    type(string) :: source
    integer  :: unit, ios, error, n, state, it
    real(R8) :: flag, flag_all, first
    logical  :: not_fatal

    not_fatal = .false.
    if (present(warn)) not_fatal = warn
    it = 0
    if (present(iter)) it = iter
    if (present(taken)) taken = .true.

    ! Rank 0 reads input.ini. state: 0 read, 1 cannot be opened, 2 empty, 3 the read of the text failed
    reason = ''
    state = 0
    if (mpi_is_root) then
      open(newunit=unit, file='input.ini', status='old', action='read', iostat=ios, iomsg=msg)
      if (ios /= 0) then
        state = 1
        reason = 'cannot be opened: '//trim(msg)
      else
        n = 0
        do
          read(unit, '(A)', iostat=ios) line
          if (ios /= 0) exit
          if (len_trim(line) > 0) n = n + 1
        end do
        close(unit)
        if (n == 0) then
          state = 2
          reason = 'is empty'
        end if
      end if
      if (state == 0) then
        ! As FiNeR reads a file (file_ini%load with a file name): read_file, then the parse of its text, below. An error
        ! of read_file is an I/O error of the read (StringiFor), not an error of the INI reader
        call source%read_file(file='input.ini', iostat=ios)
        if (ios > 0) then
          state = 3
          write(reason, '(A,I0)') 'cannot be read: I/O error ', ios
        else
          text = source%chars()
        end if
      end if
    end if
    call mpi_bcast_integer(state)
    if (state /= 0) then
      if (not_fatal) then
        if (mpi_is_root) call Write_Not_Taken(label, it, 'input.ini '//trim(reason))
        if (present(taken)) taken = .false.
        return
      end if
      if (mpi_is_root) then
        write(*,'(A,T35,A)') '   '//trim(label), 'FAIL'
        write(*,'(4X,A)') '[ERROR] input.ini '//trim(reason)
      end if
      call mpi_stop_all('input.ini cannot be read (see above)')
    end if

    ! Every rank parses the text that rank 0 read
    call mpi_bcast_chars(text)
    call fini%load(source=text, error=error)
    if (error /= 0) write(reason, '(A,I0,A)') 'cannot be read by the INI reader (FiNeR error ', error, &
      '; a line that FiNeR could not parse is named above, if any)'

    flag = 0.0_R8
    if (error /= 0) flag = 1.0_R8
    call mpi_allreduce_max_r8(flag, flag_all)
    if (flag_all > 0.0_R8) then
      call mpi_allreduce_min_r8(merge(real(mpi_rank_, R8), huge(1.0_R8), flag > 0.0_R8), first)
      if (not_fatal) then
        if (flag > 0.0_R8 .and. nint(first) == mpi_rank_) call Write_Not_Taken(label, it, 'input.ini '//trim(reason))
        if (present(taken)) taken = .false.
        return
      end if
      if (flag > 0.0_R8 .and. nint(first) == mpi_rank_) then
        write(*,'(A,T35,A)') '   '//trim(label), 'FAIL'
        write(*,'(4X,A)') '[ERROR] input.ini '//trim(reason)
      end if
      call mpi_stop_all('input.ini cannot be read (see above)')
    end if

  end subroutine Load_Input_File


  !> The warning of a reading of input.ini during the run that is not taken (ini-error = warn): the line of the set-up
  !> summary with NOT TAKEN, the iteration and the error, without the [ERROR] in front of it (the run goes on), and that
  !> no value of the reading is taken. Written out at once (flush), as the message of a stop is (mpi_stop_all).
  subroutine Write_Not_Taken(label, iter, error)
    use iso_fortran_env, only: output_unit
    implicit none
    character(len=*), intent(in) :: label   !< Line of the set-up summary: Input file, read again
    integer,          intent(in) :: iter    !< The iteration of the reading
    character(len=*), intent(in) :: error   !< The error, as a stop gives it
    ! Local
    character(len=:), allocatable :: reason

    reason = trim(adjustl(error))
    if (index(reason, '[ERROR] ') == 1) reason = reason(9:)
    write(*,'(A,T35,A)') '   '//trim(label), 'NOT TAKEN'
    write(*,'(4X,A,I0,A)') '[WARNING] input.ini read again at iteration ', iter, ' is not taken: '//reason
    write(*,'(4X,A)') 'No value of this reading is taken: the run goes on with the values it had before the reading.'
    flush(output_unit)

  end subroutine Write_Not_Taken

end module FUSS_Read_Ini
