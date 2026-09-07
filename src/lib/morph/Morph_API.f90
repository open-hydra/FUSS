!>@brief The MORPH public façade. A solver calls this module and nothing else.
!>
!> Keeping the surface this narrow is what makes plan 09's extraction checklist
!> cheap: when MORPH becomes a standalone hydra solver, everything below stays
!> and only the caller changes.
!>
!> Typical use from a solver:
!>
!>   call Morph_Init   ( node, dim, geom, status )              ! once, at setup
!>   ...
!>   call Morph_Update ( node, node_old, dim, t, dt, law, geom, status )   ! each step
!>   call Morph_GCL_Residual ( geom, dim, r_abs, r_rel )        ! assert r_rel ~ 1e-15
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_API
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Metrics,           only: Morph_Metrics_Block, Morph_Mesh_Type, Morph_Geom_Allocate
  use Morph_GCL,               only: Morph_Cell_Volume_Signed, Morph_Swept_Volumes_Cell, &
                                     Morph_GCL_Residual, Morph_GCL_Volume_Gap
  use Morph_Quality,           only: Morph_Quality_Check
  use Morph_Motion_m,          only: morph_motion_t
  use Morph_Motion_Static,     only: morph_motion_static_t
  use Morph_Motion_Prescribed, only: morph_motion_prescribed_t

  implicit none
  private

  ! Types
  public :: morph_vec3_t, morph_tens3_t, morph_face_t, morph_dir_t
  public :: morph_geom_t, morph_status_t
  public :: morph_motion_t, morph_motion_static_t, morph_motion_prescribed_t

  ! Status codes
  public :: MORPH_OK, MORPH_ERR_SINGULAR, MORPH_ERR_NEGVOL, MORPH_ERR_JACOBIAN
  public :: MORPH_ERR_SHAPE, MORPH_ERR_GCL

  ! Entry points
  public :: Morph_Init
  public :: Morph_Update
  public :: Morph_Mesh_Type
  public :: Morph_GCL_Residual
  public :: Morph_GCL_Volume_Gap
  public :: Morph_Quality_Check
  public :: Morph_Geom_Allocate

contains

  !> Build the geometry for a mesh that is not moving yet. Call once at setup.
  !> Leaves vol_old == vol and dV_swept == 0, so a GCL residual taken straight
  !> after Morph_Init is exactly zero.
  subroutine Morph_Init ( node, dim, geom, status )
    type(morph_vec3_t),   intent(in)    :: node(0:,0:,0:)
    integer(I4),          intent(in)    :: dim(3)
    type(morph_geom_t),   intent(inout) :: geom
    type(morph_status_t), intent(out)   :: status

    call status%clear()

    call Morph_Metrics_Block ( node, dim, geom, status )

    geom%vol_old  = geom%vol
    geom%dV_swept = 0.0_R8

  end subroutine Morph_Init


  !> Advance the mesh one step under `law` and rebuild every geometric quantity
  !> the ALE update needs, including swept volumes consistent with the volume
  !> change so that the discrete GCL holds by construction.
  !>
  !> On entry `node` holds the positions at time t. On exit `node` holds the
  !> positions at t+dt and `node_old` holds those at t.
  subroutine Morph_Update ( node, node_old, dim, t, dt, law, geom, status )
    type(morph_vec3_t),    intent(inout) :: node(0:,0:,0:)
    type(morph_vec3_t),    intent(inout) :: node_old(0:,0:,0:)
    integer(I4),           intent(in)    :: dim(3)
    real(R8),              intent(in)    :: t, dt
    class(morph_motion_t), intent(inout) :: law
    type(morph_geom_t),    intent(inout) :: geom
    type(morph_status_t),  intent(out)   :: status
    ! Local
    type(morph_status_t) :: law_status
    real(R8)             :: min_vol, min_jac

    call status%clear()

    if ( .not. geom%allocated_ ) then
      call status%fail ( MORPH_ERR_SHAPE, &
        'Morph_Update called before Morph_Init' )
      return
    endif

    ! Roll the time levels BEFORE the law runs.
    node_old     = node
    geom%vol_old = geom%vol

    ! Move the nodes.
    call law%apply ( node, node_old, dim, t, dt, law_status )
    if ( .not. law_status%ok() ) then
      status = law_status
      return
    endif

    ! Rebuild geometry at the new positions.
    call Morph_Metrics_Block ( node, dim, geom, status )
    call fill_swept_volumes  ( node_old, node, dim, geom )

    ! Validity. Reported, never fatal here -- the caller decides.
    if ( status%ok() ) then
      call Morph_Quality_Check ( node, dim, min_vol, min_jac, status )
    endif

  end subroutine Morph_Update


  ! ---------------------------------------------------------------------------

  subroutine fill_swept_volumes ( node_old, node_new, dim, geom )
    type(morph_vec3_t), intent(in)    :: node_old(0:,0:,0:), node_new(0:,0:,0:)
    integer(I4),        intent(in)    :: dim(3)
    type(morph_geom_t), intent(inout) :: geom
    ! Local
    integer(I4) :: i, j, k
    real(R8)    :: xo(3,8), xn(3,8), dV(6)

    !$omp parallel do collapse(3) private(i,j,k,xo,xn,dV)
    do k = 1, dim(3)
    do j = 1, dim(2)
    do i = 1, dim(1)
      call gather_cell ( node_old, i, j, k, xo )
      call gather_cell ( node_new, i, j, k, xn )

      ! EXACT inertness for cells that did not move.
      !
      ! Evaluating the swept-volume formula on identical old/new nodes does not
      ! give exactly zero: the degenerate lateral quadrilateral (a,b,b,a) has
      ! centroid 0.25*(a+b+b+a), which does not round to exactly (a+b)/2, so the
      ! cross products leave a round-off residue (~1e-21 absolute, ~1e-16
      ! relative). Small, but the zero-motion regression (plan 05 section 5.2)
      ! demands BIT-identity with the static baseline, and a nonzero ALE flux
      ! however tiny would break it.
      !
      ! The test is per cell, not per block, so a partially moving mesh still
      ! gets exact zeros everywhere the nodes are unchanged.
      if ( all( xo == xn ) ) then
        geom%dV_swept(1:6,i,j,k) = 0.0_R8
      else
        call Morph_Swept_Volumes_Cell ( xo, xn, dV )
        geom%dV_swept(1:6,i,j,k) = dV
      endif
    enddo; enddo; enddo
    !$omp end parallel do

  end subroutine fill_swept_volumes


  pure subroutine gather_cell ( node, i, j, k, x )
    type(morph_vec3_t), intent(in)  :: node(0:,0:,0:)
    integer(I4),        intent(in)  :: i, j, k
    real(R8),           intent(out) :: x(3,8)

    x(:,1) = node(i-1,j-1,k-1)%c
    x(:,2) = node(i-1,j-1,k  )%c
    x(:,3) = node(i-1,j  ,k-1)%c
    x(:,4) = node(i-1,j  ,k  )%c
    x(:,5) = node(i  ,j-1,k-1)%c
    x(:,6) = node(i  ,j-1,k  )%c
    x(:,7) = node(i  ,j  ,k-1)%c
    x(:,8) = node(i  ,j  ,k  )%c

  end subroutine gather_cell

end module Morph_API
