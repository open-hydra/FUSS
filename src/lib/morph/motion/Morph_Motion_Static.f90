!>@brief The null motion law: nodes do not move.
!>
!> Required, not a convenience. It is what the zero-motion regression gate
!> (plan 05 section 5.2) runs against: with this law selected, the whole
!> moving-mesh machinery is compiled in and exercised, yet every result must
!> still reproduce the Phase 0 baseline bit-identically. A moving-mesh
!> implementation that cannot reproduce the static answer exactly is wrong, and
!> this law is how that claim gets tested rather than asserted.
module Morph_Motion_Static
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Motion_m

  implicit none
  private

  public :: morph_motion_static_t

  type, extends(morph_motion_t) :: morph_motion_static_t
  contains
    procedure :: apply => static_apply
  end type morph_motion_static_t

contains

  subroutine static_apply ( self, node, node_old, dim, t, dt, status )
    class(morph_motion_static_t), intent(inout) :: self
    type(morph_vec3_t),           intent(inout) :: node(0:,0:,0:)
    type(morph_vec3_t),           intent(in)    :: node_old(0:,0:,0:)
    integer(I4),                  intent(in)    :: dim(3)
    real(R8),                     intent(in)    :: t, dt
    type(morph_status_t),         intent(out)   :: status

    call status%clear()

    ! Assign rather than leave untouched: the caller is entitled to assume that
    ! `node` holds the new positions after apply(), whatever the law. Copying
    ! node_old makes that contract hold for the null law too, and it is exact,
    ! so bit-identity with the static path is preserved.
    node = node_old

  end subroutine static_apply

end module Morph_Motion_Static
