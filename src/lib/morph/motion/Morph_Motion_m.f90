!>@brief The mesh-motion law abstraction.
!>
!> A motion law is the ONLY thing that decides where nodes go. Everything else
!> in MORPH -- metrics, swept volumes, the GCL identity, quality checks -- is
!> written against node positions and is therefore independent of the law.
!>
!> That separation is what makes the later phases cheap:
!>   * plan 08 drives surface recession by supplying a displacement to
!>     morph_motion_nodeshift_t; it adds no ALE or GCL code.
!>   * remeshing with conservative interpolation arrives as another subclass,
!>     again touching no ALE or GCL code.
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_Motion_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m

  implicit none
  private

  public :: morph_motion_t, morph_apply_i

  !> Abstract base. `apply` must write the NEW node positions into `node`,
  !> having been given the previous positions in `node_old`.
  type, abstract :: morph_motion_t
    character(len=32) :: name = 'unnamed'
  contains
    procedure(morph_apply_i), deferred :: apply
  end type morph_motion_t

  abstract interface
    subroutine morph_apply_i ( self, node, node_old, dim, t, dt, status )
      import :: morph_motion_t, morph_vec3_t, morph_status_t, I4, R8
      class(morph_motion_t), intent(inout) :: self
      type(morph_vec3_t),    intent(inout) :: node(0:,0:,0:)
      type(morph_vec3_t),    intent(in)    :: node_old(0:,0:,0:)
      integer(I4),           intent(in)    :: dim(3)
      real(R8),              intent(in)    :: t     !> time at the START of the step
      real(R8),              intent(in)    :: dt
      type(morph_status_t),  intent(out)   :: status
    end subroutine morph_apply_i
  end interface

end module Morph_Motion_m
