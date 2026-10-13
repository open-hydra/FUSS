module FUSS_Advanced_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Base_Types_m
  use FUSS_Parameters_m
  use FUSS_Series_Data_m
  use Lib_ORION_Data

  implicit none

  !! ------------------------------------------------------
  !! ------------------------------------------------------
  !! FUNDAMENTAL TYPES
  type :: block_type
    integer                                :: dim(3)           ! Number of cells in i-j-k (ghost not included)
    real(R8), allocatable                  :: vol(:,:,:)       ! Cell volume
    ! Moving-mesh (ALE) bookkeeping, allocated only when mesh motion is active.
    ! vol_old is the volume at the previous step; dV_swept(f,i,j,k) is the volume
    ! swept by face f of that cell over the step, signed positive outward.
    ! The discrete GCL guarantees vol - vol_old == sum_f dV_swept per cell.
    real(R8), allocatable                  :: vol_old(:,:,:)   ! Cell volume at step n
    real(R8), allocatable                  :: dV_swept(:,:,:,:)! (6,i,j,k) swept volumes
    type(FUSS_vector_3D_type), allocatable :: node(:,:,:)      ! Mesh grid points (including ghost)
    type(FUSS_vector_3D_type), allocatable :: node_old(:,:,:)  ! Mesh grid points at the previous step (moving mesh only)
    type(FUSS_vector_3D_type), allocatable :: node_disp(:,:,:) ! law = nodeshift: this step's displacement of the receding
                                                               ! face's nodes (Compute_Node_Recession; zero elsewhere)
    type(FUSS_tensor_3D_type), allocatable :: M(:,:,:)         ! Metric transformation tensor
    type(FUSS_vector_3D_type), allocatable :: dl(:,:,:)        ! Average cell length (in i/j/k direction). eg: dl%c(1) is sqrt(dx**2+dy**2+dz**2) of the cell in the i direction
    type(FUSS_d_metrics_type)              :: dir(3)           ! Direction object. Contains: i-faces, j-faces, k-faces; eg: dir(1)%face(i,j,k)%n
  end type block_type

  type :: bc_type
    integer                    :: i, j, k, b, f                             ! ijk coordinates, block and face in which the boundary element is located
    integer                    :: type                                      ! BC type (1,9)
    integer                    :: bs, is, js, ks, fs, d11, d12, d21, d22    ! BC 1 (connection) specifications
    type(FUSS_tensor_3D_type)  :: Mg(2)                                     ! Ghost cell metric tensor
    type(FUSS_vector_3D_type)  :: dlg(2)                                    ! Ghost cell average cell length
    real(R8)                   :: volg(2)                                   ! Ghost cell volume
    real(R8)                   :: Tg(6), matIDg                             ! Ghost cell primitive stencil
    real(R8), allocatable      :: chig(:)                                   ! Ghost chi stencil (nsub): allocated and zeroed with pyrolysis on, filled in plan 08
    real(R8)                   :: pg(6) = 0d0                               ! Ghost gas-pressure stencil (plan 07): filled by Fill_Ghost_p as Tg is
    integer                    :: ni(2)                                     ! BC chimera
    integer, allocatable       :: donorID(:,:)                              ! BC chimera
    real(R8), allocatable      :: volume_fraction(:)                        ! BC chimera
    real(R8)                   :: ext_flux                                  ! Multi-Solver Coupling              
  end type bc_type
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! ------------------------------------------------------
  ! SOLID EXTENSIONS
  type, extends(block_type) :: FUSS_block_type
    real(R8), dimension(:,:,:), allocatable :: T, TO           ! Primitive variables at time n and n-1
    real(R8), dimension(:,:,:), allocatable :: R               ! Residuals
    real(R8), dimension(:,:,:), allocatable :: RS1, RS2        ! Implicit smoothing residuals (temporary storage)
    real(R8), dimension(:,:,:), allocatable :: dtlocal         ! Local time step
    real(R8), dimension(:,:,:), allocatable :: qvol            ! Volumetric source term
    real(R8), dimension(:,:,:), allocatable :: matID           ! ID of material for each cell
    ! Pyrolysis state and caches (plan 06 section 3.1), allocated through the
    ! field list only when [FUSS-Pyrolysis] is present: chi(1:nsub, cells) is
    ! the component-first reaction advancement (ghost layers allocated and
    ! mirrored at wall faces, unread in phase 2), chiO its start-of-step copy,
    ! dchi the stage increment; rho_s/k_s/cp_s are the solid caches refreshed
    ! after every (T, chi) update; mdot_pg is the gas mass flux through the
    ! cell's downstream face along the outflow direction.
    real(R8), dimension(:,:,:,:), allocatable :: chi, chiO, dchi
    real(R8), dimension(:,:,:),   allocatable :: rho_s, k_s, cp_s, mdot_pg
    ! Darcy closure (plan 07 section 3.1), field-list rows allocated only
    ! under gas-closure = type2: p_pg the gas pressure (ghost layers, the
    ! restart/output field), p_pgO its start-of-step copy, p_new the sweep
    ! scratch, eta_pg / eta_prev the lagged equation-of-state storage
    ! coefficient eps_pg W/(R T) of this and the previous step, eta_st /
    ! eta_stO the coefficient of the STORED gas at the end / start of the step
    ! (G = eta_st p_pg, a state of the scheme), gam_pg the lagged
    ! face-coefficient tensor Gamma = rho_g K / mu (9 components), pi_pg the
    ! gas production used by the solve, u_pg the cell gas velocity written for
    ! output, T_prev the temperature at step n-1 from which the solve
    ! extrapolates its production rate to mid-step. eta_prev, eta_st and
    ! T_prev are restarted.
    real(R8), dimension(:,:,:),   allocatable :: p_pg, p_pgO, p_new, eta_pg, eta_prev, eta_st, eta_stO, pi_pg, T_prev
    real(R8), dimension(:,:,:,:), allocatable :: u_pg, gam_pg
  end type FUSS_block_type

  type, extends(bc_type) :: FUSS_bc_type
    ! Default-initialised: each wall BC type reads only its own
    ! subset of these, and a dispatch that names the wrong one must read 0, not
    ! whatever the heap held. -finit-derived does not reach allocatable arrays
    ! of this type, so the init sweep cannot catch that class without this.
    real(R8) :: qw = 0d0, Tw = 0d0, hconv = 0d0, qrad = 0d0, Tref = 0d0, eps = 0d0  ! BC viscous wall specifications
    ! Darcy closure (plan 07 sections 3.1, 3.4): pw the Dirichlet gas pressure
    ! of an open face (the optional 4th real of a 303/305 line, else p-wall),
    ! mdot_pg / edot_pg the gas mass [kg/m2/s] and enthalpy [W/m2] fluxes
    ! leaving through the face (outward positive), set by the solver.
    real(R8) :: pw = 0d0, mdot_pg = 0d0, edot_pg = 0d0
    ! Tw_face: the wall temperature the face's BC routine computed in the
    ! current RK stage (BC_Fluxes, 301-305), read by the type-2 gas-enthalpy
    ! flux of an open face; -1 where no wall BC has set it (the cell value is
    ! used there).
    real(R8) :: Tw_face = -1d0
    type(time_series_type) :: BCtime
    ! Ablating wall 306 (plan 08 section 2.1). Line 2 of bc.txt: lambda, T_inf
    ! and the two INPUT/ files (the B' table, registry index ibprime, 0 for
    ! none; the env table). The current env values are set once per step by
    ! Ablation_Update_Env (the env p_w goes into pw above); then the state of
    ! the last converged surface energy balance: T_w (also in Tw_face, which
    ! the type-2 open-face gas enthalpy reads), the conduction flux, the
    ! recession rate, the char mass flux, B' = B'_pg + B'_ch, the recession,
    ! the blowing correction beta (the Picard start of the next solve) and the
    ! B'_pg of the last chem = 1 solve (the fractional-chem rule); seb(1:8) =
    ! conduction, convection, pyrolysis-gas advection, char advection,
    ! alpha q_rad, eps sigma T_inf^4, -eps sigma T_w^4, residual F (seb.txt).
    real(R8) :: lambda = 0d0, Tinf = 0d0, rhoeUeCH0 = 0d0, hr = 0d0, qrad_env = 0d0, hconv_env = 0d0, &
                Tedge = 0d0, chem = 0d0
    real(R8) :: Tw_prev = 0d0, qw_last = 0d0, sdot = 0d0, mdot_ch = 0d0, Bprime = 0d0, s_total = 0d0
    real(R8) :: seb(8) = 0d0, beta = 1d0, hw = 0d0, Bpg = 0d0, Bch = 0d0, Bpg_lag = 0d0
    integer  :: seb_iter = 0, ibprime = 0, nclamp(4) = 0
    type(time_table_t) :: env
    ! Coupled ablating wall 307 (plan 08 section 2.5 ii): line 2 is the 103
    ! connectivity record, T_inf and the B' table (lambda = 0: the CFD's
    ! coefficient is already blown); env_now = p_w, rhoeUeCH, h_r, q_rad of
    ! the exchange, written by hydra before every FUSS step (standalone: by
    ! the TEST ONLY [FUSS-Ablation] coupled-env-file); p_w = 0 = never written
    real(R8) :: env_now(4) = 0d0
    ! TEST ONLY dump of a 101 face (plan 08 section 2.3 R5): the conduction
    ! flux BC_Connection put into the residual, and (type 2) the Darcy mass
    ! and enthalpy fluxes of the face, written to interface-flux.txt
    real(R8) :: q_if = 0d0, qm_if = 0d0, qh_if = 0d0
  end type FUSS_bc_type
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! ------------------------------------------------------
  ! COMPUOND TYPES
  type :: FUSS_domain_type
    real(R8)                                         :: time             ! Solution time
    real(R8)                                         :: dtglobal         ! Global dt for time accurate simulation
    integer                                          :: iter, itermax    ! Iteration number
    integer                                          :: nb, nbound       ! Number of blocks, number of boundary faces
    integer, dimension(:,:), allocatable             :: n_bf             ! Number of bc elements per faces per block
    type(FUSS_block_type), dimension(:), allocatable :: blk              ! Allocatable block type
    type(FUSS_bc_type),    dimension(:), allocatable :: bc               ! Allocatable bc object
    ! MPI local BC indices (built by build_local_bc_index)
    integer                                          :: n_local_bc = 0
    integer, dimension(:), allocatable               :: local_bc_idx
    integer                                          :: n_local_bs = 0
    integer, dimension(:), allocatable               :: local_bs_idx
  end type FUSS_domain_type

  type :: FUSS_simulation_type
    type(FUSS_domain_type), dimension(:), allocatable  :: domain
    type(ORION_data), dimension(:), allocatable        :: IOfield
  end type FUSS_simulation_type
  !! ------------------------------------------------------
  !! ------------------------------------------------------

end module FUSS_Advanced_Types_m