!>@brief Types owned by MORPH (Mesh mOtion, Remeshing and Positioning Handler).
!>
!> MORPH is a geometry component: given node positions it produces every
!> geometric quantity a finite-volume solver needs, and advances those nodes
!> under a motion law. It knows nothing about temperature, enthalpy or fluxes.
!>
!> THE ONE RULE: this component must never `use` a FUSS module. The dependency
!> arrow points one way, always, so that extracting MORPH into a standalone
!> hydra solver is a directory move rather than a rewrite. The rule is checked
!> mechanically by test/run_morph_guards.sh -- not left to discipline.
!>
!> Note on the types below: they are deliberately structurally identical to
!> FUSS's geometry types, but Fortran treats structurally identical derived
!> types as DISTINCT types, so this does not buy a copy-free pass-through
!> (as originally hoped). It buys the one-way dependency, which is the point.
!> The adapter on the FUSS side performs one explicit pack/unpack; see
!> Adapter_Morph.f90 for the cost note and the aliasing escape hatch.
module Morph_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  public

  !> Status codes returned by MORPH entry points.
  !> MORPH NEVER calls `stop`: a library that aborts cannot be embedded in
  !> another solver. The caller inspects the status and decides.
  integer(I4), parameter :: MORPH_OK             = 0
  integer(I4), parameter :: MORPH_ERR_SINGULAR   = 1  !> metric tensor determinant is zero
  integer(I4), parameter :: MORPH_ERR_NEGVOL     = 2  !> non-positive cell volume
  integer(I4), parameter :: MORPH_ERR_JACOBIAN   = 3  !> Jacobian sign flip (tangled cell)
  integer(I4), parameter :: MORPH_ERR_SHAPE      = 4  !> inconsistent array extents
  integer(I4), parameter :: MORPH_ERR_GCL        = 5  !> discrete GCL identity violated

  integer(I4), parameter :: MORPH_MSGLEN = 256

  !! ------------------------------------------------------
  !! Geometric primitives
  !! ------------------------------------------------------
  type :: morph_vec3_t
    real(R8) :: c(3) = 0.0_R8
  end type morph_vec3_t

  type :: morph_tens3_t
    real(R8) :: c(3,3) = 0.0_R8
  end type morph_tens3_t

  !> One interface: outward unit normal and area.
  type :: morph_face_t
    real(R8) :: n(3) = 0.0_R8
    real(R8) :: A    = 0.0_R8
  end type morph_face_t

  !> All interfaces normal to one computational direction.
  type :: morph_dir_t
    type(morph_face_t), allocatable :: f(:,:,:)
  end type morph_dir_t

  !! ------------------------------------------------------
  !! Status
  !! ------------------------------------------------------
  type :: morph_status_t
    integer(I4)                  :: code = MORPH_OK
    character(len=MORPH_MSGLEN)  :: message = ''
    integer(I4)                  :: i = 0, j = 0, k = 0   !> offending cell, 0 if not cell-local
  contains
    procedure :: ok    => morph_status_ok
    procedure :: fail  => morph_status_fail
    procedure :: clear => morph_status_clear
  end type morph_status_t

  !! ------------------------------------------------------
  !! The geometry bundle
  !! ------------------------------------------------------
  !> Everything the ALE state update needs, produced by one Morph_Update call.
  !>
  !> `vol_old` and `dV_swept` exist so the caller can honour the discrete
  !> Geometric Conservation Law. `dV_swept(f,i,j,k)` is the volume swept by
  !> face `f` of cell (i,j,k) over the step, signed positive when the face
  !> moves outward. Face ordering matches FUSS's BC face convention:
  !>   1 = i-low, 2 = i-high, 3 = j-low, 4 = j-high, 5 = k-low, 6 = k-high.
  !>
  !> The GCL requires, per cell and to round-off:
  !>   vol - vol_old  ==  sum_f dV_swept(f,...)
  !> Morph_GCL_Residual returns exactly that mismatch so the caller can assert.
  !> `vol` is the SIGNED divergence-theorem volume throughout: it is what the
  !> solver consumes, what telescopes against dV_swept for the Geometric
  !> Conservation Law, and what makes `vol <= 0` a genuine inversion test.
  !>
  !> There was briefly a second field holding FUSS's legacy unsigned
  !> five-tetrahedron volume, to keep the static path bit-identical with
  !> Phase 0. That legacy formula was wrong for warped cells, so it was fixed
  !> globally instead and the two fields collapsed back into one -- one volume,
  !> one definition, no way to mix them.
  type :: morph_geom_t
    integer(I4)                       :: dim(3) = 0
    real(R8),            allocatable  :: vol(:,:,:)          !> signed cell volume
    real(R8),            allocatable  :: vol_old(:,:,:)      !> signed cell volume, previous step
    real(R8),            allocatable  :: dV_swept(:,:,:,:)   !> (6, i, j, k)
    type(morph_tens3_t), allocatable  :: M(:,:,:)            !> metric transformation tensor
    type(morph_vec3_t),  allocatable  :: dl(:,:,:)           !> average cell length per direction
    type(morph_dir_t)                 :: dir(3)              !> face normals and areas
    logical                           :: allocated_ = .false.
  end type morph_geom_t

contains

  pure logical function morph_status_ok ( self ) result ( res )
    class(morph_status_t), intent(in) :: self
    res = ( self%code == MORPH_OK )
  end function morph_status_ok


  pure subroutine morph_status_clear ( self )
    class(morph_status_t), intent(inout) :: self
    self%code = MORPH_OK
    self%message = ''
    self%i = 0; self%j = 0; self%k = 0
  end subroutine morph_status_clear


  !> Record a failure. First failure wins, so a loop can call this without
  !> masking the earliest (and usually most diagnostic) problem.
  pure subroutine morph_status_fail ( self, code, message, i, j, k )
    class(morph_status_t), intent(inout) :: self
    integer(I4),           intent(in)    :: code
    character(len=*),      intent(in)    :: message
    integer(I4), optional, intent(in)    :: i, j, k

    if ( self%code /= MORPH_OK ) return   ! keep the first failure

    self%code = code
    self%message = message
    if ( present(i) ) self%i = i
    if ( present(j) ) self%j = j
    if ( present(k) ) self%k = k
  end subroutine morph_status_fail

end module Morph_Types_m
