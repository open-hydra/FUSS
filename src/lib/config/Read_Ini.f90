module FUSS_Read_Ini

  implicit none

  !> Set by Read_Inifile when input.ini carries options nothing reads. Reported
  !> by Check_Input in Wrap_Setup; kept here rather than stopping on the spot so
  !> that the input report lists every problem in one pass instead of one per
  !> run.
  character(len=1024), public :: unknown_options_message = ''

contains

  subroutine Read_Inifile()
    use Finer,               only: file_ini
    use FUSS_Read_Sim_Param, only: Register_Sim_Param
    use FUSS_Read_IO,        only: Register_IO_Fields, Register_Probes
    use FUSS_Read_Numerics,  only: Register_Numerics
    use FUSS_Backend_INI,    only: Load_Ini, Scan_Ini, Check_Unknown_Options
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

    ! Anything in a FUSS-owned section that nothing registered is a typo or a
    ! stale key, and either way it is NOT being applied. Recorded here, reported
    ! by Check_Input.
    call Check_Unknown_Options(fini, unknown_options_message)

  end subroutine Read_Inifile


  !> Re-read input.ini while the run is in progress.
  !>
  !> Only runtime-mutable parameters are applied; a change to anything else is
  !> ignored and recorded. See Reload_Runtime_Ini for the policy and why it also
  !> closes the mid-run route around the setup-time guards.
  subroutine Read_Inifile_Runtime( iter )
    use Finer,               only: file_ini
    use FUSS_Backend_INI,    only: Reload_Runtime_Ini
    use FUSS_Global_m,       only: FUSS_phase_prefix
    use FUSS_Input_Registry
    implicit none
    integer, intent(in) :: iter
    ! Local
    type(file_ini) :: fini
    character(len=4096) :: report
    integer :: u

    call fini%load(filename='input.ini')

    call Reload_Runtime_Ini(fini, iter, report)

    ! The file is created ONLY when something was actually ignored, so its mere
    ! existence after a run is the signal that input.ini was edited mid-flight
    ! and part of the edit did not take. A line in a log nobody reads would not
    ! carry that.
    if ( len_trim(report) > 0 ) then
      open(newunit=u, file='OUTPUT/'//trim(FUSS_phase_prefix)//'runtime-ini-ignored.txt', &
           status='unknown', position='append', action='write')
      write(u,'(A)') trim(report)
      close(u)
    endif

  end subroutine Read_Inifile_Runtime


end module FUSS_Read_Ini