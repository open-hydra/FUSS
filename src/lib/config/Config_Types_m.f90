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
    character(len=llen)  :: sol_format, ic_format   ! solution format (native,tecplot,vtk) and formatting (ascii,raw)
    integer              :: shell_diter  ! Shell update
    integer              :: ini_diter    ! input.ini update
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
  !! Mesh motion ------------------------------------------
  !! ------------------------------------------------------
  !> Configuration for the moving-mesh (ALE) path. The mesh update itself lives
  !> in the MORPH component; this only selects and parameterises the motion law.
  !>
  !> `law = static` is the default and is exactly inert: MORPH forces swept
  !> volumes to zero bitwise for cells that do not move, so a run with
  !> mesh-motion compiled in reproduces a non-moving run bit-for-bit.
  !>
  !> `prescribed` is the analytic test law used to verify the ALE/GCL machinery
  !> (plan 05 section 5.1). It is a function of the REFERENCE node position and
  !> time, so every MPI rank computes identical node positions without
  !> communication, and a restart lands on exactly the same geometry.
  type :: mesh_motion_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    character(len=llen) :: law          ! static | prescribed | translation
    real(R8)            :: vel(3)       ! rigid translation velocity [m/s]
    real(R8)            :: amp(3)       ! displacement amplitude per coordinate [m]
    real(R8)            :: kx(3)        ! spatial wavenumbers [1/m]
    real(R8)            :: ky(3)
    real(R8)            :: kz(3)
    real(R8)            :: omega        ! temporal angular frequency [1/s]
    real(R8)            :: gcl_tol      ! max acceptable relative GCL residual
    logical             :: taper        ! taper motion to zero on the block boundary
    ! Useful variables
    logical             :: enabled      ! .true. unless law == 'static'
  end type mesh_motion_t
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

  !! ------------------------------------------------------
  !! Pyrolysis (plan 06) ----------------------------------
  !! ------------------------------------------------------
  !> One reaction of the flattened reaction list (material-major, phase-major,
  !> reaction-minor): d chi/dt = A T^b exp(-E/(R T)) (1 - chi)^m above T_thr.
  type, public :: pyro_reaction_t
    integer  :: phase = 0          ! solid phase i of its material
    real(R8) :: F = 0.0d0          ! mass fraction of the phase converted by this reaction
    real(R8) :: A = 0.0d0          ! pre-exponential factor [1/s]
    real(R8) :: E = 0.0d0          ! activation energy [J/mol]
    real(R8) :: m = 0.0d0          ! reaction order
    real(R8) :: b = 0.0d0          ! temperature exponent (PATO's n)
    real(R8) :: T_thr = 0.0d0      ! threshold temperature [K]
    real(R8) :: h_reac = 0.0d0     ! PATO's h[i][j]: parsed, not used (second closure only)
    real(R8) :: eps_rho_F = 0.0d0  ! eps_i0 rho_i0 F [kg/m3]: sub-phase mass per unit chi
    real(R8) :: w = 0.0d0          ! eps_rho_F / (rho_v - rho_c): weight of (1 - chi) in tau
  end type pyro_reaction_t

  !> Per-material data: the constants of pyrolysis.dat and the 1 K tables of
  !> virgin.dat, char.dat and gas.dat (array bounds = integer temperature).
  type, public :: pyro_material_t
    character(len=clen) :: name = ''
    integer  :: nphase = 0
    real(R8), allocatable :: rho_I(:), eps_I(:)      ! intrinsic density and volume fraction per phase
    integer,  allocatable :: nreac_phase(:)
    real(R8) :: eps_g_v = 0.0d0, eps_g_c = 0.0d0     ! virgin and char gas porosity
    real(R8) :: rho_v = 0.0d0, rho_c = 0.0d0         ! sum eps_I rho_I, and that minus sum eps_I rho_I F
    integer  :: nsub = 0                             ! reactions of this material
    logical  :: hbar_guard = .false.                 ! rho_v = rho_c: tau == 1, the weights are undefined
    integer  :: Tlo_s = 0, Thi_s = 0                 ! virgin/char table rows [K]
    real(R8), allocatable :: cp_v(:), h_v(:), k_v(:), em_v(:), ab_v(:)
    real(R8), allocatable :: cp_c(:), h_c(:), k_c(:), em_c(:), ab_c(:)
    integer  :: Np = 0, Tlo_g = 0, Thi_g = 0         ! gas table: p-levels and rows [K]
    real(R8), allocatable :: p_lev(:)                ! ascending
    real(R8), allocatable :: W_g(:,:), h_g(:,:), mu_g(:,:)   ! (T, p-level)
    integer  :: T_brk_lo = 0, T_brk_hi = 0           ! Newton bracket [K]
  end type pyro_material_t

  type, public :: pyrolysis_t
    character(len=llen) :: warning_message = 'none'
    character(len=llen) :: error_message = 'none'
    character(len=llen) :: description = 'none'
    ! USER-DEFINED INPUTS ([FUSS-Pyrolysis], plan 06 section 3.5)
    character(len=llen) :: gas_closure = 'type1'
    character(len=llen) :: outflow_direction = 'none'
    real(R8)            :: p_gas = 101325.0d0
    logical             :: hbar_guard_report = .false.
    ! Useful variables
    logical  :: enabled = .false.                    ! <=> the section is present (Scan_Ini); docgen never sets it
    integer  :: nmat = 0, nsub = 0                   ! materials; reactions over all materials
    integer, allocatable :: mat_nsub(:)              ! reactions per material
    integer, allocatable :: mat_sub_offset(:)        ! (0:nmat): reactions of material m are offset(m-1)+1 .. offset(m)
    type(pyro_reaction_t), allocatable :: reac(:)    ! (1:nsub) flattened list
    type(pyro_material_t), allocatable :: mat(:)     ! (1:nmat)
    integer  :: out_dir = 0, out_sign = 0            ! outflow direction (1..3) and sign, from outflow-direction
  end type pyrolysis_t
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
  type(mesh_motion_t), public           :: obj_mesh_motion
  type(table_t), public                 :: obj_table
  type(pyrolysis_t), public             :: obj_pyrolysis
  !! ------------------------------------------------------
  !! ------------------------------------------------------

end module FUSS_Config_Types_m