module FUSS_Read_Ini
  
  implicit none

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

  end subroutine Read_Inifile


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

    ! Validate registry
    out = Validate_Registry()

  end subroutine Read_Inifile_Runtime


end module FUSS_Read_Ini