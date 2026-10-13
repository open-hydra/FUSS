!>@brief Surface recession: the nodes of one block face move by a given
!> displacement, the opposite face stays put, and every grid line between them
!> is compressed uniformly (plan 08 section 2.3).
!>
!> The block face that moves is named by `dir` (the grid direction of the lines
!> that end on it: 1 = i, 2 = j, 3 = k) and `side` (1 = the low-index face,
!> 2 = the high-index face). `disp(a, b)` is the displacement of the surface
!> node of the line at transverse node indices (a, b): (j, k) for dir 1,
!> (i, k) for dir 2, (i, j) for dir 3. The caller (FUSS's recession, from the
!> ablating wall's recession rate) supplies it every step; the law has no state
!> of its own and, like the translation law, displaces INCREMENTALLY from
!> node_old, so a restart needs nothing beyond the mesh in the solution file.
!>
!> A node at arc length s from the pinned end of its line, on a line of length
!> L, moves by disp * s/L. Every cell of the line therefore shrinks by the same
!> fraction |disp|/L and a graded mesh keeps its grading; a cell cannot invert
!> before the whole line is consumed. s and L are summed from node_old along
!> the line, the same sum for both, so the surface node moves by exactly disp
!> (s/L = 1 bitwise) and the pinned node by exactly nothing.
!>
!> A line whose length would drop to |disp| or less returns MORPH_ERR_NEGVOL:
!> a status, never a stop (MORPH guard (c)); the caller decides. A zero
!> displacement everywhere leaves the nodes bit-identical, which the swept-
!> volume fill needs to give exact zeros.
module Morph_Motion_Nodeshift
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Motion_m

  implicit none
  private

  public :: morph_motion_nodeshift_t

  type, extends(morph_motion_t) :: morph_motion_nodeshift_t
    integer(I4) :: dir  = 0                          !> 1..3, the direction of the receding lines
    integer(I4) :: side = 0                          !> 1 = low-index face recedes, 2 = high-index face
    type(morph_vec3_t), allocatable :: disp(:,:)     !> surface-node displacement over the transverse node indices
  contains
    procedure :: apply => nodeshift_apply
  end type morph_motion_nodeshift_t

contains

  subroutine nodeshift_apply ( self, node, node_old, dim, t, dt, status )
    class(morph_motion_nodeshift_t), intent(inout) :: self
    type(morph_vec3_t),              intent(inout) :: node(0:,0:,0:)
    type(morph_vec3_t),              intent(in)    :: node_old(0:,0:,0:)
    integer(I4),                     intent(in)    :: dim(3)
    real(R8),                        intent(in)    :: t, dt
    type(morph_status_t),            intent(out)   :: status
    ! Local
    integer(I4) :: t1, t2, a, b, m, n, bad_a, bad_b
    integer(I4) :: idx(3)
    real(R8)    :: d(3), Lline, s, seg(3)
    real(R8)    :: arc(0:maxval(dim))
    logical     :: neg, moving

    call status%clear()

    if ( self%dir < 1 .or. self%dir > 3 .or. self%side < 1 .or. self%side > 2 ) then
      call status%fail ( MORPH_ERR_SHAPE, 'nodeshift: dir must be 1..3 and side 1 or 2' )
      return
    endif
    select case ( self%dir )
      case (1) ; t1 = 2 ; t2 = 3
      case (2) ; t1 = 1 ; t2 = 3
      case default ; t1 = 1 ; t2 = 2
    end select
    if ( .not. allocated(self%disp) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'nodeshift: no displacement' )
      return
    endif
    if ( lbound(self%disp,1) /= 0 .or. ubound(self%disp,1) /= dim(t1) .or. &
         lbound(self%disp,2) /= 0 .or. ubound(self%disp,2) /= dim(t2) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'nodeshift: the displacement array does not match the receding face' )
      return
    endif

    ! node already holds node_old (Morph_Update copies before calling apply).
    ! No displacement anywhere: leave every node bit-identical.
    moving = .false.
    do b = 0, dim(t2)
      do a = 0, dim(t1)
        if ( any( self%disp(a,b)%c /= 0.0_R8 ) ) moving = .true.
      enddo
    enddo
    if ( .not. moving ) return

    n = dim(self%dir)
    neg = .false. ; bad_a = 0 ; bad_b = 0

    !$omp parallel do collapse(2) private(a,b,m,idx,d,Lline,s,seg,arc) reduction(.or.:neg)
    do b = 0, dim(t2)
    do a = 0, dim(t1)
      d = self%disp(a,b)%c
      if ( all( d == 0.0_R8 ) ) cycle
      ! arc(m) = arc length from the pinned end to the node m steps away
      ! from it, summed over node_old
      arc(0) = 0.0_R8
      do m = 1, n
        seg = point( node_old, line_idx( m, a, b ) ) - point( node_old, line_idx( m - 1, a, b ) )
        arc(m) = arc(m-1) + sqrt( dot_product( seg, seg ) )
      enddo
      Lline = arc(n)
      if ( Lline - sqrt( dot_product( d, d ) ) <= 0.0_R8 ) then
        neg = .true.
      else
        do m = 1, n
          s = arc(m)
          idx = line_idx( m, a, b )
          node(idx(1),idx(2),idx(3))%c = node_old(idx(1),idx(2),idx(3))%c + d * ( s / Lline )
        enddo
      endif
    enddo
    enddo
    !$omp end parallel do

    if ( neg ) then
      ! locate the first offending line (serial, only on failure)
      outer: do b = 0, dim(t2)
        do a = 0, dim(t1)
          d = self%disp(a,b)%c
          if ( all( d == 0.0_R8 ) ) cycle
          Lline = 0.0_R8
          do m = 1, n
            seg = point( node_old, line_idx( m, a, b ) ) - point( node_old, line_idx( m - 1, a, b ) )
            Lline = Lline + sqrt( dot_product( seg, seg ) )
          enddo
          if ( Lline - sqrt( dot_product( d, d ) ) <= 0.0_R8 ) then
            bad_a = a ; bad_b = b
            exit outer
          endif
        enddo
      enddo outer
      ! the location: the line's two transverse node indices and the direction
      call status%fail ( MORPH_ERR_NEGVOL, 'nodeshift: the displacement consumes a whole grid line', &
                         bad_a, bad_b, self%dir )
    endif

  contains

    !> (i, j, k) of the node m steps from the pinned end on line (a, b)
    pure function line_idx ( mm, aa, bb ) result ( ijk )
      integer(I4), intent(in) :: mm, aa, bb
      integer(I4) :: ijk(3)
      integer(I4) :: l
      l = mm
      if ( self%side == 1 ) l = n - mm     ! the high-index face is pinned
      ijk(self%dir) = l
      ijk(t1) = aa
      ijk(t2) = bb
    end function line_idx

    pure function point ( nd, ijk ) result ( x )
      type(morph_vec3_t), intent(in) :: nd(0:,0:,0:)
      integer(I4),        intent(in) :: ijk(3)
      real(R8) :: x(3)
      x = nd(ijk(1),ijk(2),ijk(3))%c
    end function point

  end subroutine nodeshift_apply

end module Morph_Motion_Nodeshift
