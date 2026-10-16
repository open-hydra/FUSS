module FUSS_Read_Ini

  implicit none

  !> Set by Read_Inifile when input.ini carries options nothing reads. Reported
  !> by Check_Input in Wrap_Setup; kept here rather than stopping on the spot so
  !> that the input report lists every problem in one pass instead of one per
  !> run.
  character(len=1024), public :: unknown_options_message = ''
  !> Set by Read_Inifile when a string value is longer than its variable (it would be
  !> cut without a message). Reported by Check_Input.
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
    use FUSS_Backend_INI,    only: Load_Ini, Scan_Ini, Check_Unknown_Options, Check_Value_Lengths
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
    call Register_Pyrolysis()
    call Register_Darcy()
    call Register_Ablation()

    ! Registry is built, now load the values from the ini file
    call Load_Ini(fini)

    ! Anything in a FUSS-owned section that nothing registered is a typo or a
    ! stale key, and either way it is NOT being applied. Recorded here, reported
    ! by Check_Input.
    call Check_Unknown_Options(fini, unknown_options_message)
    call Check_Value_Lengths(fini, long_values_message)

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
    use FUSS_Config_Types_m, only: obj_sim_param, obj_io
    use FUSS_Mod_MPI,        only: mpi_is_root, mpi_bcast_integer, mpi_bcast_real
    implicit none
    integer, intent(in) :: iter
    ! Local
    type(file_ini) :: fini
    character(len=4096) :: report
    integer :: u

    ! Under MPI only root reads and validates the file, and the seven runtime-
    ! mutable values are broadcast (plan 11 section 2.4). Two of them decide on
    ! every rank whether a collective runs (res-diter the residual reduction,
    ! ini-diter this routine and its broadcasts), so ranks that saw a mid-run edit
    ! at different moments would deadlock. The report is written by root and read
    ! by no other rank, so it is not broadcast.
    if ( mpi_is_root ) then
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
    endif

    ! RUNTIME_MUTABLE (Backend_INI.f90), in its order
    call mpi_bcast_real    ( obj_sim_param%res_threshold  )
    call mpi_bcast_real    ( obj_sim_param%time_threshold )
    call mpi_bcast_integer ( obj_io%sol_diter   )
    call mpi_bcast_real    ( obj_io%sol_dtime   )
    call mpi_bcast_integer ( obj_io%res_diter   )
    call mpi_bcast_integer ( obj_io%shell_diter )
    call mpi_bcast_integer ( obj_io%ini_diter   )

  end subroutine Read_Inifile_Runtime


end module FUSS_Read_Ini