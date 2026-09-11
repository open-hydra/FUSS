module FUSS_Read_Numerics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Config_Types_m
  use FUSS_Input_Registry
  
  implicit none
  private
  public :: Register_Numerics

contains

  subroutine Register_Numerics (nmgl)
    use FUSS_Parameters_m
    use IR_precision
    implicit none
    ! Local
    integer, intent(in) :: nmgl
    character(len=llen) :: section

    section = trim(codename)//'-Numerics'

    !! ------------------------------------------------------
    !! Time Scheme ------------------------------------------
    !! ------------------------------------------------------
    obj_time_scheme%warning_message = 'none'
    obj_time_scheme%error_message   = 'none'
    obj_time_scheme%description     = 'none'

    ! Solver-type
    call reg%add( trim(section), 'time-scheme', obj_time_scheme%solver_type, 'euler', 'Time integration solver', 'euler, RK2, RK3', .false. )

    ! Stability coefficients and related options
    call reg%add( trim(section), 'vnn', obj_time_scheme%vnn, '0.3', 'VNN parameter', '> 0', .false. )
    call reg%add( trim(section), 'vnn-rise-threshold', obj_time_scheme%rampa_vnn_iter, '0', 'VNN rise threshold', '>= 0', .false. )

    ! Time-accurate switch
    call reg%add( trim(section), 'time-accurate', obj_time_scheme%time_accurate, '.false.', 'Time accurate switch', 'logical', .true. )

    ! New state
    call reg%add( trim(section), 'integration-variables', obj_time_scheme%integration_variables, 'cons', 'Integration variables (cons/prim)', 'cons ,  prim', .false. )


    ! Implicit residual smoothing --------------------------
    obj_irs%description     = 'none'
    obj_irs%warning_message = 'none'
    obj_irs%error_message   = 'none'
    call reg%add( trim(section), 'irs', obj_irs%enabled, '.false.', 'Implicit Residual Smoothing', 'logical', .false. )
    call reg%add( trim(section), 'irs-beta', obj_irs%beta, '0.0', 'IRS beta parameter', '>= 0', .false. )

    !! ------------------------------------------------------
    !! ------------------------------------------------------


    !! ------------------------------------------------------
    !! Space Scheme ------------------------------------------
    !! ------------------------------------------------------

    ! Multigrid levels --------------------------------------
    call Register_Multigrid_Levels(nmgl)

    ! Mesh motion (ALE) -------------------------------------
    call Register_Mesh_Motion()

    !! ------------------------------------------------------
    !! ------------------------------------------------------

  end subroutine Register_Numerics


  !> Moving-mesh options. The default `static` law is exactly inert, so adding
  !> these parameters cannot change any existing case.
  subroutine Register_Mesh_Motion ()
    use FUSS_Config_Types_m
    use FUSS_Parameters_m
    implicit none
    character(len=llen) :: section

    section = trim(codename)//'-MeshMotion'

    obj_mesh_motion%warning_message = 'none'
    obj_mesh_motion%error_message   = 'none'
    obj_mesh_motion%description     = 'none'

    ! The allowed list is load-bearing, not documentation: `enabled` is derived
    ! as (law /= 'static'), while the adapter selects the law by name with a
    ! static default. A misspelt law would therefore switch the ALE path ON and
    ! then move nothing -- a silently wrong run. Validate_Registry rejects any
    ! name not listed here, so keep this in sync with Adapter_Update_Mesh.
    call reg%add( trim(section), 'law', obj_mesh_motion%law, 'static', &
                  'Mesh motion law', 'static ,  prescribed ,  translation', .false. )

    call reg%add( trim(section), 'vel', obj_mesh_motion%vel, '0.0', &
                  'Rigid-translation velocity [m/s]', 'real', .false. )

    call reg%add( trim(section), 'amp', obj_mesh_motion%amp, '0.0', &
                  'Prescribed-motion displacement amplitude per coordinate [m]', '>= 0', .false. )
    call reg%add( trim(section), 'kx', obj_mesh_motion%kx, '0.0', &
                  'Prescribed-motion wavenumbers multiplying x [1/m]', 'real', .false. )
    call reg%add( trim(section), 'ky', obj_mesh_motion%ky, '0.0', &
                  'Prescribed-motion wavenumbers multiplying y [1/m]', 'real', .false. )
    call reg%add( trim(section), 'kz', obj_mesh_motion%kz, '0.0', &
                  'Prescribed-motion wavenumbers multiplying z [1/m]', 'real', .false. )
    call reg%add( trim(section), 'omega', obj_mesh_motion%omega, '0.0', &
                  'Prescribed-motion angular frequency [1/s]', 'real', .false. )

    ! The discrete GCL residual is asserted every step when the mesh moves; a
    ! violation means the geometry and the state update disagree about how much
    ! volume was swept, which corrupts the solution silently.
    call reg%add( trim(section), 'gcl-tolerance', obj_mesh_motion%gcl_tol, '1.0e-10', &
                  'Max relative discrete-GCL residual before the run is stopped', '> 0', .false. )

    ! Boundary taper: only interior nodes move. Used by the conservation audit,
    ! where a fixed domain boundary makes exact energy conservation the expected
    ! result rather than an approximation.
    call reg%add( trim(section), 'taper-to-boundary', obj_mesh_motion%taper, '.false.', &
                  'Taper prescribed motion to zero on the block boundary', 'logical', .false. )

  end subroutine Register_Mesh_Motion


  subroutine Register_Multigrid_Levels(nmgl)
    use FUSS_Config_Types_m
    use FUSS_Parameters_m
    use IR_precision
    implicit none
    integer, intent(in) :: nmgl
    integer :: m
    character(len=llen) :: option

    ! Allocate iteration threshold array for each level
    obj_multigrid%MGL = max(nmgl, 1)
    allocate( obj_multigrid%iter_threshold(obj_multigrid%MGL) )
    obj_multigrid%iter_threshold = 0

    ! Register per-level iteration parameters (level-2-iter, level-3-iter, ...)
    if (.not. obj_sim_param%HYDRA_MG) then 
      do m = 1, obj_multigrid%MGL
        write(option,'(A5,I0,A5)') 'level', m, '-iter'
        call reg%add( trim(codename)//'-Multigrid', trim(option), obj_multigrid%iter_threshold(m), '0', 'Iterations for multigrid level '//trim(str(.true.,m)), '>= 0', .false. )
      enddo
    else
      do m = 1, obj_multigrid%MGL
        write(option,'(A5,I0,A5)') 'level', m, '-iter'
        call reg%add( 'HYDRA-Multigrid', trim(option), obj_multigrid%iter_threshold(m), '0', 'Iterations for multigrid level '//trim(str(.true.,m)), '>= 0', .false. )
      enddo
    endif

  end subroutine Register_Multigrid_Levels  

end module FUSS_Read_Numerics