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


    ! Load input.ini
    call fini%load(filename='input.ini')

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

    ! Load input.ini
    call fini%load(filename='input.ini')

    ! Registry is built, now load the values from the ini file
    call Load_Ini(fini)

    ! Validate registry: a bad runtime edit is reported, the value stays as read
    out = Validate_Registry()
    if (index(out,'ERROR')>0) write(*,'(A)') ' [WARNING] runtime input.ini: '//trim(out)

  end subroutine Read_Inifile_Runtime


end module FUSS_Read_Ini