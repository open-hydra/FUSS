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
    call reg%add( trim(section), 'vnn', obj_time_scheme%vnn, '0.15', &
      'Diffusion number: dt = vnn rho cp / sum_d k_d |grad xi_d|^2 (limits 0.5 euler/RK2, 0.628 RK3; ' // &
      'MOSE''s convention; until October 2026 twice this value)', '> 0', .false. )
    call reg%add( trim(section), 'vnn-rise-threshold', obj_time_scheme%rampa_vnn_iter, '0', 'VNN rise threshold', '>= 0', .false. )

    ! Time-accurate switch
    call reg%add( trim(section), 'time-accurate', obj_time_scheme%time_accurate, '.false.', 'Time accurate switch', 'logical', .true. )

    ! New state
    call reg%add( trim(section), 'integration-variables', obj_time_scheme%integration_variables, 'cons', 'Integration variables (cons/prim)', 'cons ,  prim', .false. )


    ! Implicit residual smoothing --------------------------
    obj_irs%description     = 'none'
    obj_irs%warning_message = 'none'
    obj_irs%error_message   = 'none'
    call reg%add( trim(section), 'irs-beta', obj_irs%beta, '0.0', &
                  'Implicit residual smoothing coefficient: 0 = off, > 0 switches IRS on (only at set-up; '// &
                  'the value is re-read at runtime)', '>= 0', .false. )

    !! ------------------------------------------------------
    !! ------------------------------------------------------

    ! Space scheme ------------------------------------------
    call reg%add( trim(section), 'diffusive-metric', obj_space_scheme%diffusive_metric, 'stencil', &
                  'Face metric of the diffusive fluxes: stencil takes the face Jacobian from the same differences of '// &
                  'the cell centres as the face gradient, exact for a linear field on any mesh except at the chimera '// &
                  'faces (102), which take the two-sided balance with the cell metrics with every value; mean takes the '// &
                  'mean of the two cell metrics at the interior faces (the block connections the two-sided balance, the '// &
                  'walls the cell metric) and overstates the gradient by (1+r)^2/(4r) at a cell-size ratio r; read at set-up only', &
                  'stencil, mean', .false. )

    !! ------------------------------------------------------
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
                  'Mesh motion law (nodeshift: the ablating faces 306 recede, plan 08)', &
                  'static ,  prescribed ,  translation ,  nodeshift', .false. )

    call reg%add( trim(section), 'vel', obj_mesh_motion%vel, '0.0', &
                  'Rigid-translation velocity [m/s]', '', .false. )

    call reg%add( trim(section), 'amp', obj_mesh_motion%amp, '0.0', &
                  'Prescribed-motion displacement amplitude per coordinate [m]', '>= 0', .false. )
    call reg%add( trim(section), 'kx', obj_mesh_motion%kx, '0.0', &
                  'Prescribed-motion wavenumbers multiplying x [1/m]', '', .false. )
    call reg%add( trim(section), 'ky', obj_mesh_motion%ky, '0.0', &
                  'Prescribed-motion wavenumbers multiplying y [1/m]', '', .false. )
    call reg%add( trim(section), 'kz', obj_mesh_motion%kz, '0.0', &
                  'Prescribed-motion wavenumbers multiplying z [1/m]', '', .false. )
    call reg%add( trim(section), 'omega', obj_mesh_motion%omega, '0.0', &
                  'Prescribed-motion angular frequency [1/s]', '', .false. )

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

    ! TEST ONLY red control of the receding conjugate case (plan 08 section 2.3
    ! R5): the ghost metrics of the 101 faces stay those of the first step, so
    ! the interface-flux pairing must fail.
    call reg%add( trim(section), 'skip-bc-metrics-refresh', obj_mesh_motion%skip_bc_metrics_refresh, '.false.', &
                  'TEST ONLY: do not rebuild the ghost-cell metrics after the mesh moves', 'logical', .false. )

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
        call reg%add( trim(codename)//'-Multigrid', trim(option), obj_multigrid%iter_threshold(m), '0', 'Iterations for multigrid level '//trim(str(.true.,m))//' (re-read at runtime)', '>= 0', .false. )
      enddo
      call reg%add( trim(codename)//'-Multigrid', 'prolongated-diter', obj_multigrid%prolongated_diter, '0', &
        'Coarse levels: least level iterations between two writes of field/wall-prolongated, which are always written at the end of a level (0 = with every field write; re-read at runtime)', '>= 0', .false. )
      call reg%add( trim(codename)//'-Multigrid', 'cycle', obj_multigrid%cycle, 'none', &
        'FAS cycle on the levels below the running one: none = grid sequencing only, v, w (steady runs, integration-variables cons). '// &
        'w is the robust choice: v without IRS can limit-cycle on 4 or more levels', 'none, v, w', .false. )
      call reg%add( trim(codename)//'-Multigrid', 'pre-iter', obj_multigrid%pre_iter, '2', &
        'FAS: smoothing iterations on a level before its coarse-grid visit', '>= 0', .false. )
      call reg%add( trim(codename)//'-Multigrid', 'post-iter', obj_multigrid%post_iter, '1', &
        'FAS: smoothing iterations on a level after the coarse-grid correction', '>= 0', .false. )
      call reg%add( trim(codename)//'-Multigrid', 'coarse-iter', obj_multigrid%coarse_iter, '10', &
        'FAS: smoothing iterations on the coarsest level of a cycle', '> 0', .false. )
    else
      do m = 1, obj_multigrid%MGL
        write(option,'(A5,I0,A5)') 'level', m, '-iter'
        call reg%add( 'HYDRA-Multigrid', trim(option), obj_multigrid%iter_threshold(m), '0', 'Iterations for multigrid level '//trim(str(.true.,m))//' (re-read at runtime)', '>= 0', .false. )
      enddo
      call reg%add( 'HYDRA-Multigrid', 'prolongated-diter', obj_multigrid%prolongated_diter, '0', &
        'Coarse levels: least level iterations between two writes of field/wall-prolongated, which are always written at the end of a level (0 = with every field write; re-read at runtime)', '>= 0', .false. )
    endif

  end subroutine Register_Multigrid_Levels  

end module FUSS_Read_Numerics