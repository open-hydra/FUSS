!>@brief Geometry kernels: cell volume, metric tensor, cell lengths, face
!>       normals and areas, computed from node positions alone.
!>
!> These routines were MOVED here from FUSS's Lib_Metrics.f90 (Compute_Metric_Tensor,
!> Compute_Norm_Area, Check_Mesh_Type). Three deliberate differences from the
!> originals:
!>
!>   1. Failures return a morph_status_t instead of calling `stop`. MORPH is a
!>      library; the caller decides what to do. (The originals aborted on a
!>      singular metric tensor and on non-positive volume.)
!>   2. In Morph_Norm_Area the `A == 0` test now happens BEFORE the division by
!>      A rather than after. The original divided first and overwrote the result,
!>      which is numerically identical but raises a spurious IEEE divide-by-zero
!>      and would trap under -ffpe-trap. Degenerate faces occur routinely on the
!>      collapsed k-faces of 2-D meshes, so this is not a hypothetical path.
!>   3. THE CELL VOLUME IS NOW SIGNED AND EXACT. The original five-tetrahedron
!>      abs() sum was wrong for warped cells (see Morph_Metric_Tensor). This
!>      CHANGES RESULTS on any non-planar-faced mesh, so the Phase 0 regression
!>      baseline was regenerated when it landed -- deliberately, as a correctness
!>      fix rather than an inert refactor.
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_GCL, only: Morph_Cell_Volume_Signed

  implicit none
  private

  public :: Morph_Mesh_Type
  public :: Morph_Metric_Tensor
  public :: Morph_Norm_Area
  public :: Morph_Metrics_Block
  public :: Morph_Geom_Allocate

contains

  !> Allocate (or re-use) the geometry bundle for a block of `dim` cells.
  !> Idempotent: a second call with the same dim is a no-op, so this can sit on
  !> the per-step path without churning the allocator.
  subroutine Morph_Geom_Allocate ( geom, dim )
    type(morph_geom_t), intent(inout) :: geom
    integer(I4),        intent(in)    :: dim(3)
    ! Local
    integer(I4) :: d, im, jm, km

    if ( geom%allocated_ .and. all(geom%dim == dim) ) return

    call Morph_Geom_Free ( geom )

    im = dim(1); jm = dim(2); km = dim(3)
    geom%dim = dim

    ! Cell-centred quantities carry one ghost layer, matching FUSS's convention.
    allocate( geom%vol     (0:im+1, 0:jm+1, 0:km+1) )
    allocate( geom%vol_old (0:im+1, 0:jm+1, 0:km+1) )
    allocate( geom%M       (0:im+1, 0:jm+1, 0:km+1) )
    allocate( geom%dl      (0:im+1, 0:jm+1, 0:km+1) )
    allocate( geom%dV_swept(6, 1:im, 1:jm, 1:km) )

    geom%vol      = 0.0_R8
    geom%vol_old  = 0.0_R8
    geom%dV_swept = 0.0_R8

    ! Face arrays: interfaces normal to direction d have one extra plane in d.
    do d = 1, 3
      select case (d)
      case (1); allocate( geom%dir(1)%f(0:im, 1:jm, 1:km) )
      case (2); allocate( geom%dir(2)%f(1:im, 0:jm, 1:km) )
      case (3); allocate( geom%dir(3)%f(1:im, 1:jm, 0:km) )
      end select
    enddo

    geom%allocated_ = .true.

  end subroutine Morph_Geom_Allocate


  subroutine Morph_Geom_Free ( geom )
    type(morph_geom_t), intent(inout) :: geom
    integer(I4) :: d

    if (allocated(geom%vol))      deallocate(geom%vol)
    if (allocated(geom%vol_old))  deallocate(geom%vol_old)
    if (allocated(geom%M))        deallocate(geom%M)
    if (allocated(geom%dl))       deallocate(geom%dl)
    if (allocated(geom%dV_swept)) deallocate(geom%dV_swept)
    do d = 1, 3
      if (allocated(geom%dir(d)%f)) deallocate(geom%dir(d)%f)
    enddo
    geom%allocated_ = .false.
    geom%dim = 0

  end subroutine Morph_Geom_Free


  !> Classify the mesh as 1-D / 2-D / 2-D-axisymmetric / 3-D and return the
  !> axisymmetric wedge angle. Moved verbatim from FUSS Check_Mesh_Type, except
  !> that ndir and delthe are RETURNED rather than written into solver globals.
  subroutine Morph_Mesh_Type ( node, dim, ndir, delthe )
    type(morph_vec3_t), intent(in)  :: node(0:,0:,0:)
    integer(I4),        intent(in)  :: dim(3)
    integer(I4),        intent(out) :: ndir
    real(R8),           intent(out) :: delthe
    ! Local
    real(R8) :: theta1, theta2, theta(2)
    real(R8) :: r0, r1   ! radius sqrt(y^2+z^2) at the k=0 and k=1 probe nodes
    integer(I4) :: jm, km

    jm = dim(2); km = dim(3)
    delthe = 0.0_R8

    if ( dim(3) > 1 ) then
      ! 3D
      ndir = 3
      theta1 = atan2( node(0,1,0)%c(3), node(0,1,0)%c(2) )
      theta2 = atan2( node(0,1,km)%c(3), node(0,1,km)%c(2) )
      theta(1) = theta2 - theta1
      theta1 = atan2( node(0,jm,0)%c(3), node(0,jm,0)%c(2) )
      theta2 = atan2( node(0,jm,km)%c(3), node(0,jm,km)%c(2) )
      theta(2) = theta2 - theta1
      if ( (theta(1)-theta(2)) < 1.d-5 ) then
        delthe = theta(1)
      else
        delthe = 0.0_R8
      endif

    elseif ( dim(3) == 1 .and. dim(2) == 1 ) then
      ! 1D
      ndir = 1

    else
      ! 2D
      ndir = 2
      theta1 = atan2( node(0,1,0)%c(3), node(0,1,0)%c(2) )
      theta2 = atan2( node(0,1,1)%c(3), node(0,1,1)%c(2) )
      theta(1) = theta2 - theta1
      theta1 = atan2( node(0,jm,0)%c(3), node(0,jm,0)%c(2) )
      theta2 = atan2( node(0,jm,1)%c(3), node(0,jm,1)%c(2) )
      theta(2) = theta2 - theta1
      ! For true axisymmetry the radius r=sqrt(y^2+z^2) must be conserved
      ! between k-planes. A flat-plate mesh with z(k=0)=0 at the wall (y=0)
      ! gives r=0, which would falsely satisfy the theta-difference check.
      r0 = sqrt( node(0,1,0)%c(2)**2 + node(0,1,0)%c(3)**2 )
      r1 = sqrt( node(0,1,1)%c(2)**2 + node(0,1,1)%c(3)**2 )
      if ( (theta(1)-theta(2)) < 1.d-5 .and. &
           r0 > 1.d-10 .and. abs(r0-r1) < 1.d-5*r0 ) then
        delthe = theta(1)      ! 2D axisymmetric
      else
        delthe = 0.0_R8        ! 2D planar
      endif
    endif

  end subroutine Morph_Mesh_Type


  !> Metric transformation tensor, average cell lengths and cell volume for one
  !> hexahedral cell given its 8 corner nodes.
  !>
  !> Node ordering follows the FUSS convention:
  !>   N1=(i-1,j-1,k-1) N2=(i-1,j-1,k) N3=(i-1,j,k-1) N4=(i-1,j,k)
  !>   N5=(i  ,j-1,k-1) N6=(i  ,j-1,k) N7=(i  ,j,k-1) N8=(i  ,j,k)
  pure subroutine Morph_Metric_Tensor ( N1, N2, N3, N4, N5, N6, N7, N8, M, dl, vol, status )
    type(morph_vec3_t),   intent(in)    :: N1, N2, N3, N4, N5, N6, N7, N8
    type(morph_tens3_t),  intent(out)   :: M
    type(morph_vec3_t),   intent(out)   :: dl
    real(R8),             intent(out)   :: vol
    type(morph_status_t), intent(inout) :: status
    ! Local
    integer(I4) :: h
    real(R8)    :: det, A(3,3), cofactor(3,3)
    real(R8)    :: x(3,8)

    ! A is M^-1, the inverse metric tensor:
    !   A = [ xcs, ycs, zcs ; xet, yet, zet ; xzi, yzi, zzi ]
    A(1,:) = N5%c - N1%c + N6%c - N2%c + N7%c - N3%c + N8%c - N4%c
    A(2,:) = N3%c - N1%c + N4%c - N2%c + N7%c - N5%c + N8%c - N6%c
    A(3,:) = N2%c - N1%c + N4%c - N3%c + N6%c - N5%c + N8%c - N7%c
    A = 0.25_R8 * A

    det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

    if ( abs(det) == 0.0_R8 ) then
      ! The original aborted here. MORPH reports and falls back to the identity
      ! so the caller can decide; it must not silently continue unnoticed.
      call status%fail ( MORPH_ERR_SINGULAR, 'metric tensor determinant is zero' )
      M%c = 0.0_R8
      do h = 1, 3
        M%c(h,h) = 1.0_R8
      enddo
    else
      cofactor(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))
      cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      cofactor(1,3) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))
      cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      cofactor(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))
      cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      cofactor(3,1) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))
      cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      cofactor(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))
      M%c = cofactor / det
    endif

    do h = 1, 3
      dl%c(h) = sqrt( A(h,1)**2 + A(h,2)**2 + A(h,3)**2 )
    enddo

    ! Cell volume, SIGNED, by the divergence theorem over the six bilinear
    ! faces (see Morph_GCL).
    !
    ! This replaces FUSS's original formula, which summed five tetrahedra taking
    ! the ABSOLUTE value of each. That formula was exact only for planar-faced
    ! cells: for a unit cell with one corner displaced by d it returned
    ! 1 + d/3 instead of the exact trilinear value 1 + d/4 -- a 0.8% error at
    ! d = 0.1, systematic wherever the mesh is warped. Being unsigned it also
    ! reported a positive volume for a tangled cell, so it could not be used as
    ! a validity test, and it could not telescope against signed swept volumes,
    ! so it was unusable for the Geometric Conservation Law.
    !
    ! With the signed formula all three problems go away at once: the volume is
    ! exact for trilinear hexahedra, `vol <= 0` genuinely means inverted, and
    ! the same number satisfies the GCL.
    x(:,1) = N1%c;  x(:,2) = N2%c;  x(:,3) = N3%c;  x(:,4) = N4%c
    x(:,5) = N5%c;  x(:,6) = N6%c;  x(:,7) = N7%c;  x(:,8) = N8%c

    vol = Morph_Cell_Volume_Signed( x )

    if ( vol <= 0.0_R8 ) then
      call status%fail ( MORPH_ERR_NEGVOL, 'non-positive signed cell volume (inverted cell)' )
    endif

  end subroutine Morph_Metric_Tensor


  !> Face outward unit normals and areas for all three directions.
  !>
  !> Opens its OWN OpenMP parallel region, exactly as the original did, so it
  !> must NOT be called from inside an existing parallel region.
  subroutine Morph_Norm_Area ( node, dim, geom )
    type(morph_vec3_t), intent(in)    :: node(0:,0:,0:)
    integer(I4),        intent(in)    :: dim(3)
    type(morph_geom_t), intent(inout) :: geom
    ! Local
    real(R8)    :: Ai, snix, sniy, sniz, Aj, snjx, snjy, snjz, Ak, snkx, snky, snkz
    integer(I4) :: i, j, k, im, jm, km
    real(R8)    :: d1(3), d2(3), d3(3)
    real(R8)    :: snixx, sniyy, snizz, snjxx, snjyy, snjzz, snkxx, snkyy, snkzz
    real(R8)    :: scal, signi, signj, signk

    im = dim(1); jm = dim(2); km = dim(3)

    ! ---- Orientation: establish the sign that makes normals point outward.
    ! i direction
    d1 = node(1,1,1)%c - node(1,0,0)%c
    d2 = node(1,0,1)%c - node(1,1,0)%c
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
    snixx=d3(1); sniyy=d3(2); snizz=d3(3)
    d1 = 0.25d0*( node(1,1,1)%c + node(1,0,1)%c + node(1,0,0)%c + node(1,1,0)%c )
    d2 = 0.25d0*( node(0,1,1)%c + node(0,0,1)%c + node(0,0,0)%c + node(0,1,0)%c )
    d3 = d1 - d2
    scal=d3(1)*snixx+d3(2)*sniyy+d3(3)*snizz
    signi=sign(1.d0,scal)

    ! j direction
    d1 = node(0,1,1)%c - node(1,1,0)%c
    d2 = node(1,1,1)%c - node(0,1,0)%c
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
    snjxx=d3(1); snjyy=d3(2); snjzz=d3(3)
    d1 = 0.25d0*( node(1,1,1)%c + node(0,1,1)%c + node(0,1,0)%c + node(1,1,0)%c )
    d2 = 0.25d0*( node(1,0,1)%c + node(0,0,1)%c + node(0,0,0)%c + node(1,0,0)%c )
    d3 = d1 - d2
    scal=d3(1)*snjxx+d3(2)*snjyy+d3(3)*snjzz
    signj=sign(1.d0,scal)

    ! k direction
    d1 = node(0,1,1)%c - node(1,0,1)%c
    d2 = node(0,0,1)%c - node(1,1,1)%c
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
    snkxx=d3(1); snkyy=d3(2); snkzz=d3(3)
    d1 = 0.25d0*( node(1,1,1)%c + node(0,1,1)%c + node(0,0,1)%c + node(1,0,1)%c )
    d2 = 0.25d0*( node(1,1,0)%c + node(0,1,0)%c + node(0,0,0)%c + node(1,0,0)%c )
    d3 = d1 - d2
    scal=d3(1)*snkxx+d3(2)*snkyy+d3(3)*snkzz
    signk=sign(1.d0,scal)

    ! ---- Normals and areas
    !$omp parallel private (d1,d2,d3,i,j,k,snix,sniy,sniz,Ai,Aj,snjx,snjy,snjz,Ak,snkx,snky,snkz)

    !$omp do collapse(3)
    do k = 1, km ; do j = 1, jm ; do i = 0, im
      d1 = node(i,j,k)%c   - node(i,j-1,k-1)%c
      d2 = node(i,j-1,k)%c - node(i,j,k-1)%c
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
      Ai = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
      ! Guard BEFORE dividing -- see the module header note (2).
      if (Ai == 0d0) then
        snix = 0d0; sniy = 0d0; sniz = 0d0
      else
        snix = d3(1)/Ai*signi
        sniy = d3(2)/Ai*signi
        sniz = d3(3)/Ai*signi
      end if
      geom%dir(1)%f(i,j,k)%A = Ai
      geom%dir(1)%f(i,j,k)%n = [ snix, sniy, sniz ]
    enddo; enddo; enddo

    !$omp do collapse(3)
    do k = 1, km ; do j = 0, jm ; do i = 1, im
      d1 = node(i-1,j,k)%c - node(i,j,k-1)%c
      d2 = node(i,j,k)%c   - node(i-1,j,k-1)%c
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
      Aj = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
      if (Aj == 0d0) then
        snjx = 0d0; snjy = 0d0; snjz = 0d0
      else
        snjx = d3(1)/Aj*signj
        snjy = d3(2)/Aj*signj
        snjz = d3(3)/Aj*signj
      end if
      geom%dir(2)%f(i,j,k)%A = Aj
      geom%dir(2)%f(i,j,k)%n = [ snjx, snjy, snjz ]
    enddo; enddo; enddo

    !$omp do collapse(3)
    do k = 0, km ; do j = 1, jm ; do i = 1, im
      d1 = node(i-1,j,k)%c   - node(i,j-1,k)%c
      d2 = node(i-1,j-1,k)%c - node(i,j,k)%c
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
      Ak = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
      if (Ak == 0d0) then
        snkx = 0d0; snky = 0d0; snkz = 0d0
      else
        snkx = d3(1)/Ak*signk
        snky = d3(2)/Ak*signk
        snkz = d3(3)/Ak*signk
      end if
      geom%dir(3)%f(i,j,k)%A = Ak
      geom%dir(3)%f(i,j,k)%n = [ snkx, snky, snkz ]
    enddo; enddo; enddo
    !$omp end parallel

  end subroutine Morph_Norm_Area


  !> Full geometry for one block: volume, metric tensor, cell lengths (per cell)
  !> followed by face normals and areas.
  !>
  !> Opens OpenMP regions internally; do not call from inside a parallel region.
  subroutine Morph_Metrics_Block ( node, dim, geom, status )
    type(morph_vec3_t),   intent(in)    :: node(0:,0:,0:)
    integer(I4),          intent(in)    :: dim(3)
    type(morph_geom_t),   intent(inout) :: geom
    type(morph_status_t), intent(inout) :: status
    ! Local
    integer(I4) :: i, j, k
    type(morph_status_t) :: local_status

    call Morph_Geom_Allocate ( geom, dim )

    !$omp parallel private(i,j,k) firstprivate(local_status)
    call local_status%clear()
    !$omp do collapse(3)
    do k = 1, dim(3)
    do j = 1, dim(2)
    do i = 1, dim(1)
      call Morph_Metric_Tensor (       &
             node(i-1,j-1,k-1),        &
             node(i-1,j-1,k  ),        &
             node(i-1,j  ,k-1),        &
             node(i-1,j  ,k  ),        &
             node(i  ,j-1,k-1),        &
             node(i  ,j-1,k  ),        &
             node(i  ,j  ,k-1),        &
             node(i  ,j  ,k  ),        &
             geom%M(i,j,k), geom%dl(i,j,k), geom%vol(i,j,k), local_status )
      if ( .not. local_status%ok() ) then
        local_status%i = i; local_status%j = j; local_status%k = k
      endif
    enddo; enddo; enddo
    !$omp end do
    ! Merge the per-thread status into the caller's. Serialised so the reported
    ! cell is well defined rather than whichever thread wrote last.
    !$omp critical (morph_status_merge)
    if ( .not. local_status%ok() .and. status%ok() ) then
      call status%fail ( local_status%code, local_status%message, &
                         local_status%i, local_status%j, local_status%k )
    endif
    !$omp end critical (morph_status_merge)
    !$omp end parallel

    call Morph_Norm_Area ( node, dim, geom )

  end subroutine Morph_Metrics_Block

end module Morph_Metrics
