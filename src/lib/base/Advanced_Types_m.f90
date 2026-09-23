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
  end type FUSS_block_type

  type, extends(bc_type) :: FUSS_bc_type
    ! Default-initialised: each wall BC type reads only its own
    ! subset of these, and a dispatch that names the wrong one must read 0, not
    ! whatever the heap held. -finit-derived does not reach allocatable arrays
    ! of this type, so the init sweep cannot catch that class without this.
    real(R8) :: qw = 0d0, Tw = 0d0, hconv = 0d0, qrad = 0d0, Tref = 0d0, eps = 0d0  ! BC viscous wall specifications
    type(time_series_type) :: BCtime
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