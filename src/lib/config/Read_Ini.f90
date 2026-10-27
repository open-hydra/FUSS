module FUSS_Read_Ini
  
  implicit none

  ! Warnings on keys that are set but have no effect, printed by Check_Input
  character(len=256), allocatable :: input_warnings(:)

contains

  subroutine Read_Inifile()
    use Finer,               only: file_ini
    use FUSS_Read_Sim_Param, only: Register_Sim_Param
    use FUSS_Read_IO,        only: Register_IO_Fields, Register_Probes
    use FUSS_Read_Numerics,  only: Register_Numerics
    use FUSS_Backend_INI,    only: Load_Ini, Scan_Ini
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

    ! Registry is built, now load the values from the ini file
    call Load_Ini(fini)

    ! Keys set to no effect
    call Find_Input_Warnings(fini, nmgl)

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


  subroutine Read_Inifile_Runtime()
    use Finer,               only: file_ini
    use FUSS_Backend_INI,    only: Load_Ini
    use FUSS_Input_Registry
    implicit none
    ! Local
    type(file_ini) :: fini
    character(len=1024) :: out

    ! Load input.ini: a file that cannot be read stops the run
    call Load_Input_File(fini, 'Input file, read again')

    ! Registry is built, now load the values from the ini file
    call Load_Ini(fini)

    ! Validate registry: a bad runtime edit is reported, the value stays as read
    out = Validate_Registry()
    if (index(out,'ERROR')>0) write(*,'(A)') ' [WARNING] runtime input.ini: '//trim(out)

  end subroutine Read_Inifile_Runtime

  !> Loads input.ini into fini, and stops the run if it cannot be read: a file that cannot be opened (missing, or without
  !> read permission), a file with no line that is not blank, a file whose read ends with an I/O error, or an error of
  !> the INI reader (FiNeR load gives an error code other than 0: a file with no section, a line that it cannot parse;
  !> FiNeR names such a line on stderr, without its number). Without this stop, the parameters would keep their defaults
  !> at set-up and the values of the previous reading during the run, and an empty file would end the run in the reader.
  !> Rank 0 alone reads the file, as FiNeR reads a file (string%read_file of StringiFor), and passes its text to every
  !> rank (mpi_bcast_chars); every rank, rank 0 too, gives that text to FiNeR (load from a source), so every rank takes
  !> the same values. If each rank read input.ini on its own, a file rewritten while the ranks read it during the run
  !> could give them different values, and a value that decides a collective call (sol-diter, res-diter, ini-diter)
  !> would make them call their collectives in different orders (forbidden by MPI-4.1, section 6.14), and hang. A file
  !> that rank 0 cannot read stops every rank: rank 0 writes why and
  !> passes the outcome (mpi_bcast_integer), then every rank stops with mpi_stop_all (exit code 1); an error
  !> of the INI reader, the same on every rank, is decided together too (the largest of a flag), and the first rank that
  !> found it writes it.
  subroutine Load_Input_File(fini, label)
    use Finer,           only: file_ini
    use stringifor,      only: string
    use FUSS_Mod_MPI,    only: mpi_is_root, mpi_rank_, mpi_bcast_integer, mpi_bcast_chars, mpi_allreduce_max_r8, &
                               mpi_allreduce_min_r8, mpi_stop_all
    use iso_fortran_env, only: R8 => real64
    implicit none
    type(file_ini),    intent(inout) :: fini
    character(len=*),  intent(in)    :: label   !< Line of the set-up summary: Input file, or Input file, read again
    ! Local
    character(len=512)  :: reason, msg
    character(len=1024) :: line
    character(len=:), allocatable :: text      !< The text of input.ini that rank 0 reads
    type(string) :: source
    integer  :: unit, ios, error, n, state
    real(R8) :: flag, flag_all, first

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
      if (flag > 0.0_R8 .and. nint(first) == mpi_rank_) then
        write(*,'(A,T35,A)') '   '//trim(label), 'FAIL'
        write(*,'(4X,A)') '[ERROR] input.ini '//trim(reason)
      end if
      call mpi_stop_all('input.ini cannot be read (see above)')
    end if

  end subroutine Load_Input_File

end module FUSS_Read_Ini
