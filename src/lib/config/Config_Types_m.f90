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
    logical   :: newrun           ! New run (false = restart)
    real(R8)  :: res_threshold    ! Stop below this residual
    real(R8)  :: time_threshold   ! Stop at this time
    integer   :: iter_threshold   ! Stop after this many iterations
    ! Useful variables
    integer         :: iter_general     ! Number of iteration - including all MG levels
    integer         :: iter_from_call
    real(R8)        :: time_from_call
    real(R8), allocatable :: residuotot
    integer         :: nthreads         ! Number of threads for simulation
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
  !> What a solution file says about the run that wrote it, read back by a restart (Read_vtk_tec): the kind of the
  !> run, the global iteration count and the counters that decide the next iteration of the run (see
  !> restart_counter_names in IO_Solution)
  type :: restart_counters_t
    logical  :: time_found = .true.     ! the solution holds a time (vtk TIME, Tecplot SOLUTIONTIME)
    logical  :: steady     = .false.    ! the solution is that of a steady run (negative time, or Tecplot STRANDID = 0)
    logical  :: found      = .false.    ! the solution holds all four counters below
    logical  :: resumed    = .false.    ! the counters are applied (same kind of run): Wrap_Setup
    integer  :: cycle      = 0          ! global iteration count of the run that wrote it (CYCLE)
    integer  :: iteration  = 0          ! FUSS_ITERATION: iteration count of the written level (domain%iter)
    integer  :: output_number = 0       ! FUSS_OUTPUT_NUMBER: number of the last numbered output
    integer  :: grid_level = 1          ! FUSS_GRID_LEVEL: level being computed
    real(R8) :: last_output_time = 0.d0 ! FUSS_LAST_OUTPUT_TIME: time of the last time-based output
  end type restart_counters_t

  type :: io_t
    character(len=llen)  :: warning_message
    character(len=llen)  :: error_message
    character(len=llen)  :: description
    ! USER-DEFINED INPUTS
    integer              :: sol_diter, res_diter   ! iteration interval to save the solution
    real(R8)             :: sol_dtime              ! time interval to save the solution
    logical              :: sol_overwrite      ! switch to overwrite the solution
    character(len=llen)  :: sol_format, ic_format   ! solution format (native,tecplot,vtk) and formatting (ascii,raw)
    integer              :: shell_diter  ! Shell update
    integer              :: ini_diter    ! input.ini update
    character(len=llen)  :: ini_error    ! input.ini read again with an error: warn (the reading not taken) or stop
    integer              :: timer_diter  ! Wall-clock timing report (0 = off)
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
    integer              :: output_number = 0   ! Number of the last numbered output (field<n>); restored by a restart
    type(restart_counters_t) :: restart         ! What the solution that a restart reads says about its run
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
    real(R8)            :: vnn              ! Diffusion number (von Neumann)
    integer             :: rampa_vnn_iter   ! Iterations of the vnn ramp
    logical             :: time_accurate    ! Time-accurate (global dt) run
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
    logical   :: at_setup = .false.   ! IRS on at setup: only then are its arrays allocated
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
    integer              :: prolongated_diter = 0  ! Level iterations between two prolongated writes (0: with every field write)
    character(len=llen)  :: cycle = 'none'     ! FAS cycle: none (grid sequencing), v, w
    integer              :: pre_iter  = 2      ! Smoothing iterations before the coarse-grid visit
    integer              :: post_iter = 1      ! Smoothing iterations after the correction
    integer              :: coarse_iter = 10   ! Iterations on the coarsest level of a cycle
    ! Useful variables
    integer, public :: MG_level    ! Current MG-level being solved
    logical, public :: change_MG   ! Logical variable to change MG-level
    integer, public :: prolongated_last = 0  ! Level iteration of the last prolongated write on the current level
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
    integer                               :: Tmin, Tmax   ! First and last temperature of the tables [K]
    character(len=clen), allocatable      :: name(:)
    real(R8), dimension(:,:), allocatable :: cp, rho, kappa, h
    logical, dimension(:), allocatable    :: h_monotone   ! Energy column non-decreasing in T (co_T bisects it)
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