!>@brief Rigid translation of the whole mesh at a constant velocity.
!>
!> WHY A LAW THIS TRIVIAL EARNS ITS OWN FILE
!> -----------------------------------------
!> It is the only motion for which the ALE remap's truncation error can be
!> written down in closed form BEFORE running anything, which is what plan 05
!> section 5.3 asks for ("state the tolerance derived from the discretisation
!> order, not chosen to pass").
!>
!> Under rigid translation by delta = v*dt every cell keeps its volume, so
!> sum_f dV_swept = 0 exactly, while the two faces normal to v sweep -A*delta
!> and +A*delta. The remap
!>
!>     h_new = [ h*V + sum_f h_upwind*dV_f ] / V
!>
!> then collapses to first-order upwind advection of h at Courant number
!> C = |v|*dt/dx, whose modified equation carries the numerical diffusivity
!>
!>     alpha_num = 0.5 * |v| * dx * ( 1 - C ).
!>
!> Every symbol on the right is known a priori, so the measured error in the
!> receding-surface conduction test can be COMPARED WITH A PREDICTION rather
!> than merely observed to be small. It also predicts first-order convergence
!> in dx -- which must not be misread as a defect in the GCL machinery: gate
!> 5.1 shows the geometry is exact to round-off, and the order here is the
!> upwind remap alone.
!>
!> DESIGN
!> ------
!> The displacement is INCREMENTAL (from node_old), not anchored to a cached
!> reference mesh the way Morph_Motion_Prescribed is. Two consequences:
!>
!>   * the law is stateless, so a restart needs nothing beyond the mesh that
!>     the solution file already carries. The oscillatory prescribed law cannot
!>     say the same -- see its header.
!>   * it stays rank-independent under MPI as long as every rank uses the same
!>     dt, which Set_Global_dt guarantees for a time-accurate run. Mesh motion
!>     in a non-time-accurate run is refused at input validation.
!>
!> Node-shifting -- surface nodes moving while the opposite boundary stays put,
!> with the interior redistributed -- is a DIFFERENT law and belongs with plan
!> 08, where the recession rate comes from the ablating wall boundary condition
!> instead of being prescribed. Nothing here needs to change for it to arrive.
module Morph_Motion_Translation
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Motion_m

  implicit none
  private

  public :: morph_motion_translation_t

  type, extends(morph_motion_t) :: morph_motion_translation_t
    !> Translation velocity [m/s]. Every node moves at exactly this velocity,
    !> boundary nodes included -- so unlike the prescribed law, material really
    !> does cross the domain boundary and the conservation statement is a
    !> balance against boundary fluxes rather than an exact constant.
    real(R8) :: vel(3) = [0.0_R8, 0.0_R8, 0.0_R8]
  contains
    procedure :: apply => translation_apply
  end type morph_motion_translation_t

contains

  subroutine translation_apply ( self, node, node_old, dim, t, dt, status )
    class(morph_motion_translation_t), intent(inout) :: self
    type(morph_vec3_t),                intent(inout) :: node(0:,0:,0:)
    type(morph_vec3_t),                intent(in)    :: node_old(0:,0:,0:)
    integer(I4),                       intent(in)    :: dim(3)
    real(R8),                          intent(in)    :: t, dt
    type(morph_status_t),              intent(out)   :: status
    ! Local
    integer(I4) :: i, j, k
    real(R8)    :: d(3)

    call status%clear()

    d = self%vel * dt

    ! A zero velocity must leave the nodes bit-identical, not merely unchanged
    ! to round-off: Morph_API's swept-volume fill tests old and new node
    ! positions for exact equality to decide whether a cell moved at all.
    ! Adding an exact zero would in fact be harmless here, but the early return
    ! also keeps a vel = 0 run free of any per-node work.
    if ( all( d == 0.0_R8 ) ) return

    !$omp parallel do collapse(3) private(i,j,k)
    do k = lbound(node,3), ubound(node,3)
    do j = lbound(node,2), ubound(node,2)
    do i = lbound(node,1), ubound(node,1)
      node(i,j,k)%c = node_old(i,j,k)%c + d
    enddo; enddo; enddo
    !$omp end parallel do

  end subroutine translation_apply

end module Morph_Motion_Translation
