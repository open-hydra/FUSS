!>@brief An analytic prescribed motion law, used to exercise the ALE/GCL
!>       machinery without any physics (plan 05 section 5.1).
!>
!> Design points that matter:
!>
!> * The displacement is a function of the REFERENCE node position and time,
!>   not of the previous position. Two consequences, both deliberate:
!>     - it is exactly rank-independent, so every MPI rank computes identical
!>       node positions with no communication. That is what lets Phase 1 avoid
!>       touching the per-field-hardcoded Mod_GhostExchange (plan 05 4.4a).
!>     - it cannot drift WITHIN A RUN: the mesh at time t is the same however
!>       many steps were taken to reach it, and however uneven those steps were.
!>
!>   It does NOT survive a restart, and an earlier version of this comment
!>   wrongly claimed it did. The reference is captured from whatever geometry
!>   the law is first handed; after a restart that is the moved mesh out of the
!>   solution file, so the displacement is applied on top of itself. There is no
!>   symptom -- the mesh stays self-consistent and the GCL still holds -- so the
!>   combination is refused outright in Check_Mesh_Motion_Compatibility rather
!>   than left to be discovered. Use morph_motion_translation_t, which displaces
!>   incrementally and therefore carries no state across a restart, when a
!>   restartable moving mesh is what is wanted.
!>
!> * Each coordinate uses a DIFFERENT wavenumber and phase, so the motion is
!>   genuinely three-dimensional, non-uniform, and not aligned with any single
!>   computational direction. A GCL bug that happens to cancel for
!>   axis-aligned stretching will not survive this.
!>
!> * The displacement is smooth and bounded, so cells deform without tangling
!>   provided `amp` is a modest fraction of the local cell size.
module Morph_Motion_Prescribed
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Motion_m

  implicit none
  private

  public :: morph_motion_prescribed_t

  type, extends(morph_motion_t) :: morph_motion_prescribed_t
    !> Displacement amplitude per coordinate [m]. Zero disables that component.
    real(R8) :: amp(3)   = [0.0_R8, 0.0_R8, 0.0_R8]
    !> Spatial wavenumbers [1/m], one row per displaced coordinate.
    real(R8) :: kx(3)    = [1.0_R8, 0.0_R8, 0.0_R8]
    real(R8) :: ky(3)    = [0.0_R8, 1.0_R8, 0.0_R8]
    real(R8) :: kz(3)    = [0.0_R8, 0.0_R8, 1.0_R8]
    !> Temporal angular frequency [1/s].
    real(R8) :: omega    = 1.0_R8
    !> Phase offsets so the three components are not in lockstep.
    real(R8) :: phase(3) = [0.0_R8, 2.0943951023931953_R8, 4.1887902047863905_R8]

    !> Taper the displacement to zero on the block boundary, so only INTERIOR
    !> nodes move. Used by the conservation audit: with the domain boundary
    !> fixed and the walls adiabatic, total energy must be exactly conserved no
    !> matter how violently the interior mesh moves, because every ALE flux is
    !> then internal and telescopes. Without the taper the boundary sweeps and
    !> material genuinely enters or leaves, so energy is not expected to balance.
    logical :: taper_to_boundary = .false.

    !> Reference mesh, captured on the first call.
    type(morph_vec3_t), allocatable, private :: node_ref(:,:,:)
    logical,                         private :: have_ref = .false.
  contains
    procedure :: apply => prescribed_apply
  end type morph_motion_prescribed_t

contains

  subroutine prescribed_apply ( self, node, node_old, dim, t, dt, status )
    class(morph_motion_prescribed_t), intent(inout) :: self
    type(morph_vec3_t),               intent(inout) :: node(0:,0:,0:)
    type(morph_vec3_t),               intent(in)    :: node_old(0:,0:,0:)
    integer(I4),                      intent(in)    :: dim(3)
    real(R8),                         intent(in)    :: t, dt
    type(morph_status_t),             intent(out)   :: status
    ! Local
    integer(I4) :: i, j, k, c
    real(R8)    :: xr(3), tnew, arg, w
    real(R8), parameter :: PI = 3.14159265358979323846_R8

    call status%clear()

    ! Capture the reference mesh once. node_old on the first call is the mesh as
    ! read from the initial-condition file.
    if ( .not. self%have_ref ) then
      allocate( self%node_ref, source = node_old )
      self%have_ref = .true.
    endif

    if ( any( shape(self%node_ref) /= shape(node) ) ) then
      call status%fail ( MORPH_ERR_SHAPE, &
        'prescribed motion: node array shape changed since the reference was captured' )
      return
    endif

    tnew = t + dt

    !$omp parallel do collapse(3) private(i,j,k,c,xr,arg,w)
    do k = lbound(node,3), ubound(node,3)
    do j = lbound(node,2), ubound(node,2)
    do i = lbound(node,1), ubound(node,1)
      xr = self%node_ref(i,j,k)%c

      ! Boundary taper. sin(pi*xi) is exactly zero at xi = 0; at xi = 1 the
      ! computed value is not (sin(PI) = 1.2e-16; PI*30/30 rounds low, 5.7e-16),
      ! so a node there moves by about 1e-18 m, which rounds away against
      ! coordinates of order 0.1 m (measured margin in conservation_4blk: 4.2x
      ! below half an ulp, plan 11 stage 2). Two blocks that share such a node
      ! (one at xi = 1, the other at xi = 0) therefore agree bitwise by that
      ! margin, not by construction.
      ! Degenerate directions (dim == 1, e.g. the k direction of a 2-D mesh) are
      ! skipped rather than tapered to zero everywhere.
      w = 1.0_R8
      if ( self%taper_to_boundary ) then
        if ( dim(1) > 1 ) w = w * sin( PI * real(i,R8) / real(dim(1),R8) )
        if ( dim(2) > 1 ) w = w * sin( PI * real(j,R8) / real(dim(2),R8) )
        if ( dim(3) > 1 ) w = w * sin( PI * real(k,R8) / real(dim(3),R8) )
      endif

      do c = 1, 3
        arg = self%kx(c)*xr(1) + self%ky(c)*xr(2) + self%kz(c)*xr(3) &
              + self%omega*tnew + self%phase(c)
        node(i,j,k)%c(c) = xr(c) + w * self%amp(c) * sin(arg)
      enddo
    enddo; enddo; enddo
    !$omp end parallel do

  end subroutine prescribed_apply

end module Morph_Motion_Prescribed
