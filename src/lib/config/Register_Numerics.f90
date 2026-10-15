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
    call reg%add( trim(section), 'dt-face-k', obj_time_scheme%dt_face_k, '.false.', &
                  'Cell time step from the larger face conductivity of each direction, evaluated as the fluxes '// &
                  'evaluate it (face-k), instead of the cell k: needed where a block holds a material jump, at which '// &
                  'the first cells of the less conductive material exchange heat through up to their neighbour''s k. '// &
                  'Changes no step on one material at constant k; a second pass over the cells with up to 9 property '// &
                  'lookups each (re-read at runtime)', 'logical', .false. )

    ! Time-accurate switch
    call reg%add( trim(section), 'time-accurate', obj_time_scheme%time_accurate, '.false.', 'Time accurate switch', 'logical', .true. )

    ! New state
    call reg%add( trim(section), 'integration-variables', obj_time_scheme%integration_variables, 'cons', 'Integration variables (cons/prim)', 'cons ,  prim', .false. )


    ! Implicit residual smoothing --------------------------
    obj_irs%description     = 'none'
    obj_irs%warning_message = 'none'
    obj_irs%error_message   = 'none'
    call reg%add( trim(section), 'irs', obj_irs%enabled, '.false.', 'Implicit Residual Smoothing (switched on only at setup)', 'logical', .false. )
    call reg%add( trim(section), 'irs-beta', obj_irs%beta, '0.0', 'IRS beta parameter', '>= 0', .false. )

    !! ------------------------------------------------------
    !! ------------------------------------------------------


    !! ------------------------------------------------------
    !! Space Scheme ------------------------------------------
    !! ------------------------------------------------------

    ! Conductivity of the interior diffusive faces ------------------
    call reg%add( trim(section), 'face-k', obj_space_scheme%face_k, 'mean-T', &
                  'Conductivity of the interior diffusive faces: mean-T interpolates the lower cell''s material '// &
                  'at the mean face temperature, once per face; cell-mean interpolates each cell once, at its own '// &
                  'temperature and material, and averages the two cells at the face (symmetric at material '// &
                  'interfaces, fewer lookups). Connection and chimera faces always use the two-sided interface '// &
                  'balance; walls always use the inner cell', 'mean-T, cell-mean', .false. )

    ! Multigrid levels --------------------------------------
    call Register_Multigrid_Levels(nmgl)

    !! ------------------------------------------------------
    !! ------------------------------------------------------

  end subroutine Register_Numerics


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