!>@brief Discrete Geometric Conservation Law: signed cell volumes and
!>       face-swept volumes that satisfy
!>
!>           V^{n+1} - V^n  ==  sum_f dV_swept(f)                        (GCL)
!>
!>       exactly, to floating-point round-off, for ARBITRARY node motion.
!>
!> Morph_Cell_Volume_Signed below is THE cell volume for the whole code:
!> Morph_Metric_Tensor calls it, so the volume the solver uses and the volume
!> that telescopes against the swept volumes are the same number by construction.
!>
!> It replaced FUSS's original five-tetrahedron abs() sum, which was exact only
!> for planar-faced cells, could not detect a tangled cell, and could not
!> telescope. Morph_GCL_Volume_Gap is retained as a diagnostic for comparing
!> against that legacy formula.
!>
!> HOW THE IDENTITY IS MADE EXACT
!> ------------------------------
!> Every volume here is computed with the divergence-theorem form
!>
!>     V = (1/3) * sum_triangles ( x_centroid . n ) * area
!>
!> over a closed triangulated surface. The space-time region between the cell at
!> t^n and t^{n+1} is bounded by: the 6 faces at t^n (inward), the 6 faces at
!> t^{n+1} (outward), and, for each cell edge, one lateral quadrilateral swept by
!> that edge. Summing the six per-face swept volumes therefore gives
!> V^{n+1} - V^n PROVIDED:
!>
!>   (a) a face is triangulated identically when it bounds the cell and when it
!>       bounds the swept region, and
!>   (b) each lateral quadrilateral is triangulated identically from both of the
!>       two faces that share its edge, so the two contributions cancel.
!>
!> (b) is the subtle one: the diagonal must be chosen PER EDGE, not per face, or
!> the cancellation fails and the GCL is violated by O(mesh velocity). We fix the
!> diagonal by a deterministic rule on the edge's node indices -- see
!> swept_quad_volume.
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_GCL
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m

  implicit none
  private

  public :: Morph_Cell_Volume_Signed
  public :: Morph_Swept_Volumes_Cell
  public :: Morph_GCL_Residual
  public :: Morph_GCL_Volume_Gap
  public :: Morph_Swept_Outflow_Ratio

  !> Local node numbering used throughout this module, matching the ordering
  !> that Morph_Metric_Tensor expects:
  !>   1=(i-1,j-1,k-1) 2=(i-1,j-1,k) 3=(i-1,j,k-1) 4=(i-1,j,k)
  !>   5=(i  ,j-1,k-1) 6=(i  ,j-1,k) 7=(i  ,j,k-1) 8=(i  ,j,k)
  !>
  !> Faces, in the FUSS face convention, given as the 4 corner nodes.
  !> NB: this traversal is INWARD-right-handed; see ORIENT below.
  !>   1 = i-low  : 1,3,4,2        4 = j-high : 3,7,8,4
  !>   2 = i-high : 5,6,8,7        5 = k-low  : 1,5,7,3
  !>   3 = j-low  : 1,2,6,5        6 = k-high : 2,4,8,6
  integer(I4), parameter :: FACE_NODE(4,6) = reshape( [ &
       1, 3, 4, 2,   &   ! face 1, i-low
       5, 6, 8, 7,   &   ! face 2, i-high
       1, 2, 6, 5,   &   ! face 3, j-low
       3, 7, 8, 4,   &   ! face 4, j-high
       1, 5, 7, 3,   &   ! face 5, k-low
       2, 4, 8, 6 ], &   ! face 6, k-high
       [4,6] )

  !> The traversal order in FACE_NODE above is INWARD-right-handed (verified:
  !> it yields -1 for a unit cube). Rather than hand-reverse all six rows, which
  !> is easy to get subtly wrong, the orientation is corrected once here. Both
  !> the cell volume and the swept volumes apply it, so the GCL telescoping is
  !> unaffected -- it is a global sign, not a per-face one.
  real(R8), parameter :: ORIENT = -1.0_R8

contains

  !> Signed volume of a hexahedral cell from its 8 corner nodes.
  !> Positive for a right-handed (untangled) cell; negative if inverted.
  pure function Morph_Cell_Volume_Signed ( x ) result ( vol )
    real(R8), intent(in) :: x(3,8)
    real(R8)             :: vol
    ! Local
    integer(I4) :: f

    vol = 0.0_R8
    do f = 1, 6
      vol = vol + quad_flux( x(:,FACE_NODE(1,f)), x(:,FACE_NODE(2,f)), &
                             x(:,FACE_NODE(3,f)), x(:,FACE_NODE(4,f)) )
    enddo
    vol = ORIENT * vol / 3.0_R8

  end function Morph_Cell_Volume_Signed


  !> The six face-swept volumes for one cell, given its nodes at the old and new
  !> time levels. Signed: positive when the face sweeps outward, i.e. when it
  !> increases the cell volume.
  !>
  !> By construction sum(dV) == V_new - V_old to round-off; assert it with
  !> Morph_GCL_Residual rather than trusting this comment.
  pure subroutine Morph_Swept_Volumes_Cell ( x_old, x_new, dV )
    real(R8), intent(in)  :: x_old(3,8), x_new(3,8)
    real(R8), intent(out) :: dV(6)
    ! Local
    integer(I4) :: f, e, n1, n2
    real(R8)    :: acc

    do f = 1, 6
      acc = 0.0_R8

      ! New face, outward.
      acc = acc + quad_flux( x_new(:,FACE_NODE(1,f)), x_new(:,FACE_NODE(2,f)), &
                             x_new(:,FACE_NODE(3,f)), x_new(:,FACE_NODE(4,f)) )
      ! Old face, reversed orientation (it bounds the swept region inward).
      acc = acc - quad_flux( x_old(:,FACE_NODE(1,f)), x_old(:,FACE_NODE(2,f)), &
                             x_old(:,FACE_NODE(3,f)), x_old(:,FACE_NODE(4,f)) )

      ! Lateral quadrilaterals, one per face edge. Each is shared with the
      ! neighbouring face of the same cell in the opposite orientation, so these
      ! cancel exactly in sum(dV) -- provided the diagonal is chosen per EDGE.
      do e = 1, 4
        n1 = FACE_NODE(e,f)
        n2 = FACE_NODE(merge(1, e+1, e == 4), f)
        acc = acc + swept_quad_volume( x_old(:,n1), x_old(:,n2), &
                                       x_new(:,n2), x_new(:,n1), n1, n2 )
      enddo

      dV(f) = ORIENT * acc / 3.0_R8
    enddo

  end subroutine Morph_Swept_Volumes_Cell


  !> Per-cell GCL residual: | (V_new - V_old) - sum_f dV_swept |.
  !> Should be at round-off (relative to the cell volume) if the implementation
  !> is correct. Exposed so the solver can assert on it directly.
  subroutine Morph_GCL_Residual ( geom, dim, resid_max, resid_rel_max )
    type(morph_geom_t), intent(in)  :: geom
    integer(I4),        intent(in)  :: dim(3)
    real(R8),           intent(out) :: resid_max      !> absolute, m^3
    real(R8),           intent(out) :: resid_rel_max  !> normalised by cell volume
    ! Local
    integer(I4) :: i, j, k
    real(R8)    :: r, vscale

    resid_max     = 0.0_R8
    resid_rel_max = 0.0_R8

    !$omp parallel do collapse(3) private(i,j,k,r,vscale) &
    !$omp   reduction(max: resid_max, resid_rel_max)
    do k = 1, dim(3)
    do j = 1, dim(2)
    do i = 1, dim(1)
      r = abs( ( geom%vol(i,j,k) - geom%vol_old(i,j,k) ) &
               - sum( geom%dV_swept(1:6,i,j,k) ) )
      vscale = max( abs(geom%vol(i,j,k)), abs(geom%vol_old(i,j,k)) )
      resid_max = max( resid_max, r )
      if ( vscale > 0.0_R8 ) resid_rel_max = max( resid_rel_max, r / vscale )
    enddo; enddo; enddo
    !$omp end parallel do

  end subroutine Morph_GCL_Residual


  !> Largest fraction of its own volume that any cell sweeps OUT in this step.
  !>
  !> This is the moving-mesh analogue of a Courant number, and it is the exact
  !> positivity condition for the Lagrange+remap update rather than a heuristic.
  !> The remap gives cell c the coefficient
  !>
  !>     V_old + sum over faces with dV < 0 of dV
  !>
  !> on its own enthalpy. Once the volume swept out exceeds V_old that
  !> coefficient turns negative, the update starts extrapolating instead of
  !> averaging, and the scheme is unconditionally unstable however exact the GCL
  !> is. The GCL residual cannot see this: a mesh that moves much too far in one
  !> step is still perfectly self-consistent.
  !>
  !> Returns >= 0; the caller should refuse to continue at >= 1.
  subroutine Morph_Swept_Outflow_Ratio ( geom, dim, ratio, ci, cj, ck )
    type(morph_geom_t), intent(in)  :: geom
    integer(I4),        intent(in)  :: dim(3)
    real(R8),           intent(out) :: ratio
    integer(I4),        intent(out) :: ci, cj, ck   !> where the maximum was found
    ! Local
    integer(I4) :: i, j, k, f
    real(R8)    :: out_vol, r, vold

    ratio = 0.0_R8
    ci = 0; cj = 0; ck = 0

    ! Not parallelised: it is a reduction that has to carry the LOCATION of the
    ! maximum, the loop is only run when the mesh moves, and it costs a single
    ! pass over cells against the several the metric rebuild already does.
    do k = 1, dim(3)
    do j = 1, dim(2)
    do i = 1, dim(1)
      out_vol = 0.0_R8
      do f = 1, 6
        if ( geom%dV_swept(f,i,j,k) < 0.0_R8 ) out_vol = out_vol - geom%dV_swept(f,i,j,k)
      enddo

      vold = abs( geom%vol_old(i,j,k) )
      if ( vold <= 0.0_R8 ) cycle

      r = out_vol / vold
      if ( r > ratio ) then
        ratio = r
        ci = i; cj = j; ck = k
      endif
    enddo; enddo; enddo

  end subroutine Morph_Swept_Outflow_Ratio


  !> Largest relative disagreement between the signed volume used by the GCL and
  !> the legacy absolute-value volume used by the static path. Diagnostic only:
  !> it should be at round-off for untangled cells, and large if a cell inverted.
  pure function Morph_GCL_Volume_Gap ( x, vol_legacy ) result ( gap )
    real(R8), intent(in) :: x(3,8)
    real(R8), intent(in) :: vol_legacy
    real(R8)             :: gap
    ! Local
    real(R8) :: vs

    vs = Morph_Cell_Volume_Signed( x )
    if ( abs(vol_legacy) > 0.0_R8 ) then
      gap = abs( vs - vol_legacy ) / abs(vol_legacy)
    else
      gap = abs( vs - vol_legacy )
    endif

  end function Morph_GCL_Volume_Gap


  ! ---------------------------------------------------------------------------
  ! Internals
  ! ---------------------------------------------------------------------------

  !> 3 * (volume contribution) of a planar-or-warped quadrilateral face
  !> (a,b,c,d) taken in right-handed outward order, via the divergence theorem.
  !>
  !> The quad is split into 4 triangles about its centroid rather than 2 about a
  !> diagonal. That keeps the result independent of any diagonal choice, which
  !> removes an entire class of face/swept-region inconsistency, and is exact for
  !> a bilinear (warped) face.
  pure function quad_flux ( a, b, c, d ) result ( f3 )
    real(R8), intent(in) :: a(3), b(3), c(3), d(3)
    real(R8)             :: f3
    ! Local
    real(R8) :: m(3)

    m = 0.25_R8 * ( a + b + c + d )
    f3 =   tri_flux( a, b, m ) &
         + tri_flux( b, c, m ) &
         + tri_flux( c, d, m ) &
         + tri_flux( d, a, m )

  end function quad_flux


  !> 3 * (volume contribution) of triangle (p,q,r), right-handed outward:
  !>   ( centroid . normal ) * area  ==  (1/6) * (p+q+r) . ((q-p) x (r-p))
  !> The 1/2 from the cross-product area and the 1/3 from the centroid are folded
  !> in here; the caller divides the total by 3 once. Net factor 1/6.
  pure function tri_flux ( p, q, r ) result ( f3 )
    real(R8), intent(in) :: p(3), q(3), r(3)
    real(R8)             :: f3
    ! Local
    real(R8) :: u(3), v(3), n(3)

    u = q - p
    v = r - p
    n(1) = u(2)*v(3) - u(3)*v(2)
    n(2) = u(3)*v(1) - u(1)*v(3)
    n(3) = u(1)*v(2) - u(2)*v(1)

    f3 = 0.5_R8 * dot_product( (p + q + r) / 3.0_R8, n )

  end function tri_flux


  !> 3 * (volume contribution) of a lateral quadrilateral swept by one cell edge.
  !>
  !> The quad (a,b,c,d) is a^n -> b^n -> b^{n+1} -> a^{n+1}. Unlike quad_flux we
  !> must guarantee that the two faces sharing this edge produce EXACTLY equal
  !> and opposite contributions. Splitting about the centroid already achieves
  !> that (the centroid of the 4 nodes is independent of traversal direction and
  !> tri_flux is exactly antisymmetric under reversing a triangle), so the
  !> per-edge node indices are accepted only to make the deterministic ordering
  !> explicit and to allow a future diagonal-based variant to stay consistent.
  pure function swept_quad_volume ( a, b, c, d, n1, n2 ) result ( f3 )
    real(R8),    intent(in) :: a(3), b(3), c(3), d(3)
    integer(I4), intent(in) :: n1, n2
    real(R8)                :: f3

    ! Traverse from the lower-numbered node so that the two faces sharing this
    ! edge generate the identical triangulation in opposite orientation.
    if ( n1 < n2 ) then
      f3 =  quad_flux( a, b, c, d )
    else
      f3 = -quad_flux( d, c, b, a )
    endif

  end function swept_quad_volume

end module Morph_GCL
