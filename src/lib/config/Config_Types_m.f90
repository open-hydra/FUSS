module FUSS_Config_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Parameters_m

  implicit none
  private

  !! ------------------------------------------------------
  !! Simulation Parameters --------------------------------
  !! ------------------------------------------------------
  type :: simulation_parameters_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    logical   :: newrun           ! Flag di restart
    real(R8)  :: res_threshold    ! Residuo min. per arresto esecuzione
    real(R8)  :: time_threshold   ! Tempo max per arresto esecuzione 
    integer   :: iter_threshold            ! Numero max iterate per arresto esecuzione
    ! Useful variables
    integer         :: iter_general     ! Number of iteration - including all MG levels
    integer         :: iter_from_call
    real(R8)        :: time_from_call
    real(R8), allocatable :: residuotot
    integer         :: nthreads         ! Number of threads for simulation
    real(R8)        :: cputime(2)       ! Simulation time duration
    integer         :: TODO             ! Decide solve and/or postprocess
    logical         :: HYDRA_time_accurate = .false.
    logical         :: HYDRA_postprocess   = .false.
    logical         :: HYDRA_MG = .false.
  end type simulation_parameters_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Input-Output -----------------------------------------
  !! ------------------------------------------------------
  type :: io_t
    character(len=llen)  :: warning_message
    character(len=llen)  :: error_message
    character(len=llen)  :: description
    ! USER-DEFINED INPUTS
    integer              :: sol_diter, res_diter   ! iteration interval to save the solution
    real(R8)             :: sol_dtime              ! time interval to save the solution
    logical              :: sol_overwrite      ! switch to overwrite the solution
    character(len=llen)  :: sol_format, ini_format   ! solution format (native,tecplot,vtk) and formatting (ascii,raw)
    integer              :: shell_diter  ! Shell update
    integer              :: ini_diter    ! input.ini update
    ! Useful variables
    character(len=llen)  :: nameinit    ! Initial file name
    character(len=llen)  :: namesource  ! Source term file name
    logical              :: write_thermo, write_transport, write_composition
    character(len=hlen)  :: Ovarnames, ORANSname
    integer              :: Onvar
    real(R8)             :: IOtime
    logical              :: write_wall
    character(len=hlen)  :: Owallnames
    integer              :: Onwall
    integer              :: unitRES             ! Residuals file unit
    character(len=llen)  :: unitRES_format      ! Format for residual_history_file
  end type io_t

  type :: io_probes_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    real(R8)            :: dtime
    integer             :: diter
    character(len=clen) :: file
    character(len=hlen) :: varnames
    integer             :: iloc(4)
    real(R8)            :: loc(3)
    ! Useful variables
    ! ...
  end type io_probes_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! NUMERICAL SCHEME !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Time Scheme ------------------------------------------
  !! ------------------------------------------------------
  type :: time_scheme_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    character(len=llen) :: solver_type      ! Solver type (explicit/implicit)
    real(R8)            :: vnn              ! Parametro di stabilita' diffusiva
    integer             :: rampa_vnn_iter   ! Numero di iterazioni per rampa di vnn
    logical             :: time_accurate    ! Flag per integrazione time-accurate
    integer             :: n_RK
    character(len=llen) :: integration_variables ! Integration variables (cons/prim)
    ! Useful variables
    ! ...
  end type time_scheme_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Implicit residual smoothing --------------------------
  !! ------------------------------------------------------
  type irs_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DRFINED INPUTS
    real(R8)  :: beta
    ! Useful variables
    logical   :: enabled
  end type irs_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Multigrid --------------------------------------------
  !! ------------------------------------------------------
  type :: multigrid_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    integer              :: MGL                ! Number of multigrid levels
    integer, allocatable :: iter_threshold(:)  ! Number of iterations for each level
    ! Useful variables
    integer, public :: MG_level    ! Current MG-level being solved
    logical, public :: change_MG   ! Logical variable to change MG-level
  end type multigrid_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! BC ---------------------------------------------------
  !! ------------------------------------------------------
  type :: io_bc_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    ! ...
    ! Useful variables
    logical, allocatable :: viscous_flag(:,:)
    logical, allocatable :: coupling_flag(:,:)
  end type io_bc_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!! OTHER TYPES !!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Properties Table -------------------------------------
  !! ------------------------------------------------------
  type :: table_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    ! ...
    ! Useful variables
    integer                               :: n
    character(len=clen), allocatable      :: name(:)
    real(R8), dimension(:,:), allocatable :: cp, rho, kappa, h
  end type table_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(simulation_parameters_t), public :: obj_sim_param
  type(io_t), public                    :: obj_io
  type(io_probes_t), allocatable, public:: obj_io_probes(:)
  type(io_bc_t), public                 :: obj_io_bc
  type(time_scheme_t), public           :: obj_time_scheme
  type(irs_t), public                   :: obj_irs
  type(multigrid_t), public             :: obj_multigrid
  type(table_t), public                 :: obj_table
  !! ------------------------------------------------------
  !! ------------------------------------------------------

end module FUSS_Config_Types_m