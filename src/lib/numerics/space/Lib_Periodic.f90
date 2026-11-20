!> Periodic faces (BC 201).
!>
!> A 201 record has the payload of a block connection 101 (bs is js ks fs d11 d12 d21 d22): face f of cell (i,j,k) of
!> block b takes as ghost cells the cells across face fs of cell (is,js,ks) of block bs, with the same index map as a
!> connection, as AFUS's periodic BC (type 6 of bcdef_file) and ATLAS/BCB's `type = periodic` write it: the two faces
!> are the two ends of a sector (rotational periodicity, as the side faces of one channel of a cooling jacket) or of a
!> period (translational periodicity), and the face of one end is a boundary of the other and vice versa.
!>
!> The values are those of a connection (the ghost cells copy the temperature and the material of the cells of the
!> other face, one flux per face computed by its owner and given to both cells, the exchanges of the connections
!> between ranks): IO_BC reads a 201 as a 101 and marks it with per /= 0. The positions are not: the cells of face fs
!> lie one period away from the ghost cells of face f. Every position of a cell of the other face that the geometry of
!> a connection takes is moved by the rigid motion x -> R x + t of the pair of faces, which takes face fs onto face f:
!> the ghost nodes of BC_Connect_Metrics (so the ghost metric Mg, its lengths and volume, of diffusive-metric = mean
!> as of the across-face column of the stencil metric), the ghost centre and the centres of the tangential neighbours
!> of the stencil face metric (Lib_Stencil_Metrics). With the moved positions a periodic face is, for both metrics, the
!> interior face that it is in the full (unrolled) mesh: a linear field along the faces, any field periodic with the
!> period of the pair, gives the fluxes of the full mesh.
!>
!> The motion of each pair of faces (b, f) <- (bs, fs) is fitted here, on every multigrid level, to points of the
!> records: the centres of the two faces of each record, and the points half a cell across them (outward of face f,
!> inward of face fs: the ghost cell is the image of the cell inside face fs). The fit is the least-squares rigid
!> motion of B. K. P. Horn (Closed-form solution of absolute orientation using unit quaternions, J. Opt. Soc. Am. A 4,
!> 1987): the rotation is the unit quaternion of the largest eigenvalue of a symmetric 4x4 matrix built from the
!> cross-covariance of the two sets (Jacobi rotations here), without reflections. No axis, angle or period is given in
!> input: they come out of the mesh, and are printed at set-up. The run stops when
!>  - the faces are not images of each other by a rigid motion with the index map of the records: a fitted point is
!>    farther than MISFIT_MAX of the half cell from the image of its partner (wrong pair of faces, or cells listed in
!>    another order on the two faces);
!>  - the points do not fix the motion (a pair of one-cell faces: two points on a line);
!>  - a face is periodic one way only, or the two directions of a pair do not invert each other: a point of the face
!>    moved there and back lands farther than MISFIT_MAX of the half cell from itself.
!> The records of all blocks and the nodes of all blocks are on every rank during the setup of the metrics: every rank
!> fits the same motions, and the stops are taken by every rank (mpi_stop_all).
module FUSS_Lib_Periodic
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_Periodic, Periodic_Move

  real(R8), parameter :: MISFIT_MAX = 1d-3    ! largest misfit of a fitted point, and largest distance of a point of the
                                              ! face from its image there and back, in units of the smallest half cell

contains

  !> Position x of a cell of the other face of a periodic pair moved onto the side of the face of the record.
  pure function Periodic_Move ( pair, x ) result ( y )
    use FUSS_Advanced_Types_m, only: periodic_pair_type
    implicit none
    type(periodic_pair_type), intent(in) :: pair
    real(R8), intent(in) :: x(3)
    real(R8) :: y(3)

    y = matmul ( pair % R, x ) + pair % t

  end function Periodic_Move


  !> Pairs of periodic faces of the domain of one multigrid level (domain%per) and the pair of each periodic record
  !> (bc%per): the motion of each pair fitted to the records (see the module). Called by Setup_Metrics before the ghost
  !> metrics, on every rank.
  subroutine Setup_Periodic ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Mod_MPI, only: mpi_stop_all, mpi_is_root
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer, allocatable  :: key(:,:), np(:), rev(:)
    real(R8), allocatable :: pc(:,:), qc(:,:), S(:,:,:), back(:)
    integer  :: r, g, n, ng, gi
    real(R8) :: p(3,2), q(3,2), h, Rot(3,3), dev
    character(len=400) :: msg

    if ( allocated(domain % per) ) deallocate ( domain % per )
    if ( .not. any ( domain % bc(:) % per /= 0 ) ) then
      allocate ( domain % per(0) )
      return
    end if

    ! 1. The pairs of faces: one per (b, f, bs, fs) of the periodic records, in the order of their first record
    allocate ( key(4, 6 * domain % nb * 6) )
    ng = 0
    do r = 1, domain % nbound
      if ( domain % bc(r) % per == 0 ) cycle
      associate ( bc => domain % bc(r) )
        gi = 0
        do g = 1, ng
          if ( all ( key(:,g) == [ bc % b, bc % f, bc % bs, bc % fs ] ) ) then
            gi = g
            exit
          end if
        end do
        if ( gi == 0 ) then
          if ( ng == size(key,2) ) call mpi_stop_all ( '[ERROR] periodic faces (201): too many pairs of faces' )
          ng = ng + 1
          key(:,ng) = [ bc % b, bc % f, bc % bs, bc % fs ]
          gi = ng
        end if
        bc % per = gi
      end associate
    end do
    allocate ( domain % per(ng) )
    do g = 1, ng
      domain % per(g) % b  = key(1,g) ; domain % per(g) % f  = key(2,g)
      domain % per(g) % bs = key(3,g) ; domain % per(g) % fs = key(4,g)
      domain % per(g) % hmin = huge(1d0)
    end do

    ! 2. Centroids of the two sets of points of each pair (p: face fs, the source; q: face f, the target)
    allocate ( np(ng), pc(3,ng), qc(3,ng), S(3,3,ng) )
    np = 0 ; pc = 0d0 ; qc = 0d0 ; S = 0d0
    do r = 1, domain % nbound
      g = domain % bc(r) % per
      if ( g == 0 ) cycle
      call Record_Points ( domain, r, p, q, h )
      np(g) = np(g) + 2
      pc(:,g) = pc(:,g) + p(:,1) + p(:,2)
      qc(:,g) = qc(:,g) + q(:,1) + q(:,2)
      domain % per(g) % hmin = min ( domain % per(g) % hmin, h )
    end do
    do g = 1, ng
      pc(:,g) = pc(:,g) / np(g)
      qc(:,g) = qc(:,g) / np(g)
    end do

    ! 3. Cross-covariance of the two sets: S(a,c) = sum (p_a - pc_a) (q_c - qc_c)
    do r = 1, domain % nbound
      g = domain % bc(r) % per
      if ( g == 0 ) cycle
      call Record_Points ( domain, r, p, q, h )
      do n = 1, 2
        S(:,:,g) = S(:,:,g) + spread ( p(:,n) - pc(:,g), 2, 3 ) * spread ( q(:,n) - qc(:,g), 1, 3 )
      end do
    end do

    ! 4. Rotation (Horn) and translation of each pair
    do g = 1, ng
      call Horn_Rotation ( S(:,:,g), Rot, dev )
      if ( dev < 0d0 ) then
        write(msg,'(A,4(I0,A))') '[ERROR] periodic faces (201): block ', key(1,g), ' face ', key(2,g), ' <- block ', &
          key(3,g), ' face ', key(4,g), ': the records do not fix the motion of the faces (one-cell faces)'
        call mpi_stop_all ( trim(msg) )
      end if
      domain % per(g) % R = Rot
      domain % per(g) % t = qc(:,g) - matmul ( Rot, pc(:,g) )
    end do

    ! 5. Misfit of the fitted points
    do r = 1, domain % nbound
      g = domain % bc(r) % per
      if ( g == 0 ) cycle
      call Record_Points ( domain, r, p, q, h )
      do n = 1, 2
        domain % per(g) % misfit = max ( domain % per(g) % misfit, &
                                         norm2 ( Periodic_Move ( domain % per(g), p(:,n) ) - q(:,n) ) )
      end do
    end do

    ! 6. The opposite direction of each pair, (bs, fs) <- (b, f), and the distance of the points of face f from their
    !    image there and back (the motion of the opposite direction, then the motion of the pair): zero to round-off when
    !    the two motions invert each other
    allocate ( rev(ng), back(ng) )
    rev = 0 ; back = 0d0
    do g = 1, ng
      do n = 1, ng
        if ( all ( key(:,n) == [ key(3,g), key(4,g), key(1,g), key(2,g) ] ) ) rev(g) = n
      end do
    end do
    do r = 1, domain % nbound
      g = domain % bc(r) % per
      if ( g == 0 ) cycle
      if ( rev(g) == 0 ) cycle
      call Record_Points ( domain, r, p, q, h )
      do n = 1, 2
        back(g) = max ( back(g), norm2 ( Periodic_Move ( domain % per(g), &
                                         Periodic_Move ( domain % per(rev(g)), q(:,n) ) ) - q(:,n) ) )
      end do
    end do

    ! 7. Summary (level 1, rank 0) and checks
    if ( mpi_is_root .and. domain % mg_level == 1 ) call Print_Pairs ( domain % per )
    do g = 1, ng
      associate ( pr => domain % per(g) )
        if ( pr % misfit > MISFIT_MAX * pr % hmin ) then
          write(msg,'(A,4(I0,A),ES10.3,A,ES10.3,A)') '[ERROR] periodic faces (201): block ', pr % b, ' face ', pr % f, &
            ' <- block ', pr % bs, ' face ', pr % fs, ': not images by a rigid motion (misfit ', pr % misfit, &
            ' m, half cell ', pr % hmin, ' m): wrong pair of faces, or cells in another order on the two faces'
          call mpi_stop_all ( trim(msg) )
        end if
        if ( rev(g) == 0 ) then
          write(msg,'(A,4(I0,A))') '[ERROR] periodic faces (201): block ', pr % b, ' face ', pr % f, ' <- block ', &
            pr % bs, ' face ', pr % fs, ': no periodic record of the opposite direction (a face is periodic one way only)'
          call mpi_stop_all ( trim(msg) )
        end if
        if ( back(g) > MISFIT_MAX * pr % hmin ) then
          write(msg,'(A,4(I0,A),ES10.3,A,ES10.3,A)') '[ERROR] periodic faces (201): block ', pr % b, ' face ', pr % f, &
            ' <- block ', pr % bs, ' face ', pr % fs, ': the motions of the two directions are not inverse (a point moved '// &
            'there and back lands ', back(g), ' m away, half cell ', pr % hmin, ' m)'
          call mpi_stop_all ( trim(msg) )
        end if
      end associate
    end do

  end subroutine Setup_Periodic


  !> The two pairs of points of periodic record r: p on face fs of the other block, q on face f of the block of the
  !> record. Point 1: the centres of the two faces. Point 2: half a cell across them, outward of face f (where the ghost
  !> cell is) and inward of face fs (where the cell of the other block is), h the distance of the centre of face f from
  !> the centre of its cell.
  subroutine Record_Points ( domain, r, p, q, h )
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    integer, intent(in)   :: r
    real(R8), intent(out) :: p(3,2), q(3,2), h
    ! Local
    real(R8) :: xf(3), nf(3), xs(3), ns(3), xc(3)

    associate ( bc => domain % bc(r) )
      call Face_Centre_Outward ( domain % blk(bc % b)  % node, bc % i,  bc % j,  bc % k,  bc % f,  xf, nf, xc )
      h = norm2 ( xf - xc )
      call Face_Centre_Outward ( domain % blk(bc % bs) % node, bc % is, bc % js, bc % ks, bc % fs, xs, ns, xc )
    end associate
    q(:,1) = xf
    p(:,1) = xs
    q(:,2) = xf + h * nf
    p(:,2) = xs - h * ns

  end subroutine Record_Points


  !> Centre (mean of the four nodes) and outward unit normal of face Fm of cell (i,j,k) of a block, and the centre of the
  !> cell (mean of its eight nodes). The normal of a face of zero area is zero.
  pure subroutine Face_Centre_Outward ( node, i, j, k, Fm, xf, nf, xc )
    use FUSS_Base_Types_m, only: FUSS_vector_3D_type
    implicit none
    type(FUSS_vector_3D_type), intent(in) :: node(0:,0:,0:)
    integer, intent(in)   :: i, j, k, Fm
    real(R8), intent(out) :: xf(3), nf(3), xc(3)
    ! Local
    real(R8) :: p1(3), p2(3), p3(3), p4(3), d1(3), d2(3), v(3)
    integer  :: a

    select case (Fm)
      case (1:2)
        a = i - mod(Fm,2)
        p1 = node(a,j-1,k-1) % c ; p2 = node(a,j,k-1) % c ; p3 = node(a,j,k) % c ; p4 = node(a,j-1,k) % c
      case (3:4)
        a = j - mod(Fm,2)
        p1 = node(i-1,a,k-1) % c ; p2 = node(i,a,k-1) % c ; p3 = node(i,a,k) % c ; p4 = node(i-1,a,k) % c
      case default
        a = k - mod(Fm,2)
        p1 = node(i-1,j-1,a) % c ; p2 = node(i,j-1,a) % c ; p3 = node(i,j,a) % c ; p4 = node(i-1,j,a) % c
    end select
    xf = 0.25d0 * ( p1 + p2 + p3 + p4 )
    xc = 0.125d0 * ( node(i-1,j-1,k-1) % c + node(i,j-1,k-1) % c + node(i-1,j,k-1) % c + node(i,j,k-1) % c &
                   + node(i-1,j-1,k  ) % c + node(i,j-1,k  ) % c + node(i-1,j,k  ) % c + node(i,j,k  ) % c )
    d1 = p3 - p1
    d2 = p4 - p2
    v  = [ d1(2)*d2(3) - d1(3)*d2(2), d1(3)*d2(1) - d1(1)*d2(3), d1(1)*d2(2) - d1(2)*d2(1) ]
    nf = 0d0
    if ( norm2(v) > 0d0 ) nf = v / norm2(v)
    ! The cross product follows the cyclic order of the face nodes, whose side depends on the face: outward is the side
    ! away from the centre of the cell
    if ( dot_product ( nf, xf - xc ) < 0d0 ) nf = -nf

  end subroutine Face_Centre_Outward


  !> Rotation R that best takes the set p onto the set q (q - qc ~ R (p - pc)), from their cross-covariance
  !> S(a,c) = sum (p_a - pc_a)(q_c - qc_c) (Horn 1987): the unit quaternion of the largest eigenvalue of the symmetric
  !> matrix N(S). dev is the gap of the two largest eigenvalues relative to the largest, negative when they coincide
  !> (to 1e-10): the sets do not fix the rotation.
  subroutine Horn_Rotation ( S, R, dev )
    implicit none
    real(R8), intent(in)  :: S(3,3)
    real(R8), intent(out) :: R(3,3), dev
    ! Local
    real(R8) :: N(4,4), V(4,4), lam(4), e(4)
    integer  :: m(1), i, l2

    N(1,1) = S(1,1) + S(2,2) + S(3,3)
    N(1,2) = S(2,3) - S(3,2) ; N(1,3) = S(3,1) - S(1,3) ; N(1,4) = S(1,2) - S(2,1)
    N(2,2) = S(1,1) - S(2,2) - S(3,3)
    N(2,3) = S(1,2) + S(2,1) ; N(2,4) = S(3,1) + S(1,3)
    N(3,3) = - S(1,1) + S(2,2) - S(3,3)
    N(3,4) = S(2,3) + S(3,2)
    N(4,4) = - S(1,1) - S(2,2) + S(3,3)
    N(2,1) = N(1,2) ; N(3,1) = N(1,3) ; N(4,1) = N(1,4) ; N(3,2) = N(2,3) ; N(4,2) = N(2,4) ; N(4,3) = N(3,4)

    call Jacobi_Eigen4 ( N, lam, V )
    m = maxloc ( lam )
    l2 = 0
    do i = 1, 4
      if ( i == m(1) ) cycle
      if ( l2 == 0 ) then
        l2 = i
      else if ( lam(i) > lam(l2) ) then
        l2 = i
      end if
    end do
    dev = ( lam(m(1)) - lam(l2) ) / max ( abs ( lam(m(1)) ), tiny(1d0) )
    if ( dev <= 1d-10 ) dev = -1d0

    e = V(:,m(1)) / norm2 ( V(:,m(1)) )     ! (w, x, y, z)
    R(1,1) = e(1)**2 + e(2)**2 - e(3)**2 - e(4)**2
    R(1,2) = 2d0 * ( e(2)*e(3) - e(1)*e(4) )
    R(1,3) = 2d0 * ( e(2)*e(4) + e(1)*e(3) )
    R(2,1) = 2d0 * ( e(2)*e(3) + e(1)*e(4) )
    R(2,2) = e(1)**2 - e(2)**2 + e(3)**2 - e(4)**2
    R(2,3) = 2d0 * ( e(3)*e(4) - e(1)*e(2) )
    R(3,1) = 2d0 * ( e(2)*e(4) - e(1)*e(3) )
    R(3,2) = 2d0 * ( e(3)*e(4) + e(1)*e(2) )
    R(3,3) = e(1)**2 - e(2)**2 - e(3)**2 + e(4)**2

  end subroutine Horn_Rotation


  !> Eigenvalues lam and eigenvectors (columns of V) of a symmetric 4x4 matrix A: cyclic Jacobi rotations until the
  !> off-diagonal part is below 1e-15 of the norm of A (W. H. Press et al., Numerical Recipes, sec. 11.1).
  subroutine Jacobi_Eigen4 ( A0, lam, V )
    implicit none
    real(R8), intent(in)  :: A0(4,4)
    real(R8), intent(out) :: lam(4), V(4,4)
    ! Local
    real(R8) :: A(4,4), off, nrm, theta, t, c, s, apq, app, aqq, akp, akq
    integer  :: sweep, p, q, k

    A = A0
    V = Identity4()
    nrm = sqrt ( sum ( A**2 ) )
    do sweep = 1, 100
      off = 0d0
      do p = 1, 3
        do q = p + 1, 4
          off = off + A(p,q)**2
        end do
      end do
      if ( sqrt(off) <= 1d-15 * nrm .or. nrm == 0d0 ) exit
      do p = 1, 3
        do q = p + 1, 4
          apq = A(p,q)
          if ( apq == 0d0 ) cycle
          app = A(p,p) ; aqq = A(q,q)
          theta = ( aqq - app ) / ( 2d0 * apq )
          t = sign ( 1d0, theta ) / ( abs(theta) + sqrt ( theta**2 + 1d0 ) )
          c = 1d0 / sqrt ( t**2 + 1d0 )
          s = t * c
          do k = 1, 4
            if ( k == p .or. k == q ) cycle
            akp = A(k,p) ; akq = A(k,q)
            A(k,p) = c * akp - s * akq ; A(p,k) = A(k,p)
            A(k,q) = s * akp + c * akq ; A(q,k) = A(k,q)
          end do
          A(p,p) = app - t * apq
          A(q,q) = aqq + t * apq
          A(p,q) = 0d0 ; A(q,p) = 0d0
          do k = 1, 4
            akp = V(k,p) ; akq = V(k,q)
            V(k,p) = c * akp - s * akq
            V(k,q) = s * akp + c * akq
          end do
        end do
      end do
    end do
    do k = 1, 4
      lam(k) = A(k,k)
    end do

  end subroutine Jacobi_Eigen4


  !> One line per pair of periodic faces: the motion (rotation angle and axis, translation) and the misfit.
  subroutine Print_Pairs ( per )
    use FUSS_Advanced_Types_m, only: periodic_pair_type
    implicit none
    type(periodic_pair_type), intent(in) :: per(:)
    ! Local
    integer  :: g
    real(R8) :: ang, w(3), ax(3)

    write(*,*)
    write(*,'(A)') ' Periodic faces (201): x(face f) = R x(face fs) + t'
    do g = 1, size(per)
      associate ( R => per(g) % R )
        ang = acos ( max ( -1d0, min ( 1d0, 0.5d0 * ( R(1,1) + R(2,2) + R(3,3) - 1d0 ) ) ) ) * 180d0 / acos(-1d0)
        w = [ R(3,2) - R(2,3), R(1,3) - R(3,1), R(2,1) - R(1,2) ]
        ax = 0d0
        if ( norm2(w) > 0d0 ) ax = w / norm2(w)
        write(*,'(A,I0,A,I0,A,I0,A,I0,A,F11.6,A,3F9.5,A,3ES11.3,A,ES9.2,A)') '   block ', per(g) % b, ' face ', per(g) % f, &
          ' <- block ', per(g) % bs, ' face ', per(g) % fs, ': rotation ', ang, ' deg about (', ax, '), t = (', &
          per(g) % t, ') m, misfit ', per(g) % misfit / per(g) % hmin, ' half cell'
      end associate
    end do

  end subroutine Print_Pairs


  pure function Identity4 () result ( I4m )
    real(R8) :: I4m(4,4)
    integer :: k
    I4m = 0d0
    do k = 1, 4
      I4m(k,k) = 1d0
    end do
  end function Identity4

end module FUSS_Lib_Periodic
