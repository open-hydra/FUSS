!> Face metrics of the diffusive fluxes with diffusive-metric = stencil.
!>
!> The face gradient takes, of the temperature, one difference across the face (the two cells) and one per tangential
!> direction (the mean of the central differences of the two cells), D T = (D_a, D_b, D_c), and maps them to Cartesian
!> components with a face metric: grad T = D T . M_f (Diffusive_Flux). Here M_f comes from the same differences of the
!> cell positions: J = [D_a x, D_b x, D_c x] (columns), M_f = J^-1. For a linear field T = T0 + g . x every difference
!> is D T = J^T g, so the face gradient is g on any mesh: skewed, stretched, with jumps of size, folded at the faces
!> between blocks. "In order to obtain the identity tensor numerically ... the two difference approximations must be the
!> same" (M. Vinokur, An analysis of finite-difference and finite-volume formulations of conservation laws, NASA
!> CR-177416, 1986, p. 26; J. Comput. Phys. 81, 1989). The metric of two cells, [1/2 (A_L + A_R)]^-1, has the right
!> column across the face but not along it: its tangential columns are the mean edges of the two cells, not the
!> differences of the centres that D_b and D_c take, so on skewed, non-smooth meshes the face gradient of a linear field
!> is wrong. (Port of MOSE's Lib_Stencil_Metrics, commit 000d28a, conduction subset: one scalar, no groups of
!> quantities.)
!>
!> Columns kept from the mean of the two cell metrics, so that its bits do not change where nothing has to: the column
!> across the face (1/2 (a_L + a_R) is x_R - x_L with the centre as the mean of the eight nodes), and the columns of the
!> index directions with one cell (a 1-D case, the k direction of a 2-D mesh), where the field does not vary.
!>
!> Positions: the centre of a cell is the mean of its eight nodes; the first layer of ghost cells takes the position at
!> which its value is right for a linear field: the centre of the source cell across a face between blocks (101);
!> behind a wall (301-305) and a chimera face (102) 2 x_f - x1 (x_f the centre of the face: Ghost_Wall's 2 Tw - T1
!> with the wall value at the face centre; ATLAS builds the receiver cell on linearly extrapolated ghost nodes, centred
!> at 2 x_f - x1); behind a symmetry face (300) and a null face (0, adiabatic, the same ghost copy) the mirror image of
!> x1 across the face (x1 itself on a face of zero area, an axis); behind 103 and 104, which carry a flux and copy the
!> cell, x1. The centres of every block are computed on every rank, from the nodes every rank holds until the end of
!> the setup of the metrics, so no message is needed; the face metrics are stored for the blocks of the rank only
!> (blk%fm(d)%Mf, interior faces) and in the boundary records (bc%Mf of the faces between blocks).
module FUSS_Lib_Stencil_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_Face_Metrics

  !> Cell centres of one block with one layer of ghost cells: x(1:3, 0:n1+1, 0:n2+1, 0:n3+1)
  type :: centres_t
    real(R8), allocatable :: x(:,:,:,:)
  end type centres_t

contains

  !> Centre of cell (i,j,k) of a block: the mean of its eight nodes.
  pure function Cell_Centre ( node, i, j, k ) result ( x )
    use FUSS_Base_Types_m, only: FUSS_vector_3D_type
    implicit none
    type(FUSS_vector_3D_type), intent(in) :: node(0:,0:,0:)
    integer, intent(in) :: i, j, k
    real(R8) :: x(3)

    x = 0.125d0 * ( node(i-1,j-1,k-1) % c + node(i,j-1,k-1) % c + node(i-1,j,k-1) % c + node(i,j,k-1) % c &
                  + node(i-1,j-1,k  ) % c + node(i,j-1,k  ) % c + node(i-1,j,k  ) % c + node(i,j,k  ) % c )

  end function Cell_Centre


  !> Centre (mean of the four nodes) and unit normal (direction of growing index) of face Fm of cell (i,j,k). The normal
  !> of a face of zero area (its nodes on an axis) is zero, as in Compute_Norm_Area: Mirror then gives the point itself.
  pure subroutine Face_Centre_Normal ( node, i, j, k, Fm, xf, nf )
    use FUSS_Base_Types_m, only: FUSS_vector_3D_type
    implicit none
    type(FUSS_vector_3D_type), intent(in) :: node(0:,0:,0:)
    integer, intent(in) :: i, j, k, Fm
    real(R8), intent(out) :: xf(3), nf(3)
    ! Local
    real(R8) :: p1(3), p2(3), p3(3), p4(3), d1(3), d2(3), v(3)
    integer :: a

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
    d1 = p3 - p1
    d2 = p4 - p2
    v = [ d1(2)*d2(3) - d1(3)*d2(2), d1(3)*d2(1) - d1(1)*d2(3), d1(1)*d2(2) - d1(2)*d2(1) ]
    if ( norm2(v) > 0d0 ) then
      nf = v / norm2(v)
    else
      nf = 0d0
    end if

  end subroutine Face_Centre_Normal


  !> Mirror image of x across the plane of a face (centre xf, unit normal nf).
  pure function Mirror ( x, xf, nf ) result ( y )
    implicit none
    real(R8), intent(in) :: x(3), xf(3), nf(3)
    real(R8) :: y(3)

    y = x - 2d0 * dot_product ( x - xf, nf ) * nf

  end function Mirror


  !> Positions of the cell centres (every block) and the face metrics of the blocks of this rank, for every multigrid
  !> level (Setup_Metrics calls it with the domain of the level, after the ghost metrics bc%Mg).
  subroutine Setup_Face_Metrics ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Lib_BC_Fluxes_Connection, only: Pair_Connections
    use FUSS_Config_Types_m, only: obj_space_scheme
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    type(centres_t), allocatable :: xc(:)
    integer :: b, i, j, k, n(3), d, r

    allocate ( xc(domain % nb) )

    ! 0. Partner and owner of every connection record (one flux per face, Connection_Fluxes)
    call Pair_Connections ( domain )

    ! diffusive-metric = mean needs no stencil metric: the mean of the two cell metrics at the interior faces, the
    ! two-sided balance at the connections and the cell metric at the walls, all from the cell metrics
    if ( .not. obj_space_scheme % stencil_metric ) return

    ! 1. Centres of the cells of every block
    do b = 1, domain % nb
      n = domain % blk(b) % dim
      allocate ( xc(b) % x(3, 0:n(1)+1, 0:n(2)+1, 0:n(3)+1) )
      xc(b) % x = 0d0
      do k = 1, n(3) ; do j = 1, n(2) ; do i = 1, n(1)
        xc(b) % x(:,i,j,k) = Cell_Centre ( domain % blk(b) % node, i, j, k )
      enddo ; enddo ; enddo
    enddo

    ! 2. Ghost cells of the faces that are not between blocks, every block
    do r = 1, domain % nbound
      if ( domain % bc(r) % type /= 101 ) call Ghost_Position ( domain, xc, r )
    enddo

    ! 3. Ghost cells of the faces between blocks: the centre of the source cell (computed in 1, not a ghost)
    do r = 1, domain % nbound
      if ( domain % bc(r) % type == 101 ) call Connection_Ghost_Position ( domain, xc, r )
    enddo

    ! 4. Face metrics of the blocks of this rank
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n = domain % blk(b) % dim
      if ( .not. allocated(domain % blk(b) % fm(1) % Mf) ) then
        allocate ( domain % blk(b) % fm(1) % Mf(max(n(1)-1,0), n(2), n(3)) )
        allocate ( domain % blk(b) % fm(2) % Mf(n(1), max(n(2)-1,0), n(3)) )
        allocate ( domain % blk(b) % fm(3) % Mf(n(1), n(2), max(n(3)-1,0)) )
      end if
      do d = 1, 3
        do k = 1, size ( domain % blk(b) % fm(d) % Mf, 3 )
        do j = 1, size ( domain % blk(b) % fm(d) % Mf, 2 )
        do i = 1, size ( domain % blk(b) % fm(d) % Mf, 1 )
          domain % blk(b) % fm(d) % Mf(i,j,k) % c = 0d0
        enddo ; enddo ; enddo
      enddo
      call Interior_Face_Metrics ( b, domain % blk(b), xc(b) )
    enddo
    do r = 1, domain % nbound
      if ( .not. is_local_block ( domain % bc(r) % b ) ) cycle
      select case ( domain % bc(r) % type )
        case (101)
          if ( domain % bc(r) % owner ) call Connection_Face_Metric ( domain, xc, r )   ! only the owner computes the face
        case (301:305)
          call Wall_Face_Metric ( domain, xc, r )
      end select
    enddo

    deallocate ( xc )

  end subroutine Setup_Face_Metrics


  !> Position of the first ghost cell of boundary record r (not a face between blocks).
  subroutine Ghost_Position ( domain, xc, r )
    use FUSS_Advanced_Types_m
    use FUSS_Parameters_m, only: guide
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    type(centres_t), intent(inout) :: xc(:)
    integer, intent(in) :: r
    ! Local
    integer :: b, Im, Jm, Km, Fm, Ig, Jg, Kg
    real(R8) :: x1(3), xf(3), nf(3)

    b  = domain % bc(r) % b
    Im = domain % bc(r) % i ; Jm = domain % bc(r) % j ; Km = domain % bc(r) % k
    Fm = domain % bc(r) % f
    Ig = Im - guide(Fm,1) ; Jg = Jm - guide(Fm,2) ; Kg = Km - guide(Fm,3)
    x1 = xc(b) % x(:,Im,Jm,Km)

    select case ( domain % bc(r) % type )

      case (0, 300)                             ! adiabatic (copy of the cell, mirror ghost): mirror image across the face
        call Face_Centre_Normal ( domain % blk(b) % node, Im, Jm, Km, Fm, xf, nf )
        xc(b) % x(:,Ig,Jg,Kg) = Mirror ( x1, xf, nf )

      case (301:305, 102)                       ! walls (2 Tw - T1) and chimera receivers (ATLAS): 2 x_f - x1
        call Face_Centre_Normal ( domain % blk(b) % node, Im, Jm, Km, Fm, xf, nf )
        xc(b) % x(:,Ig,Jg,Kg) = 2d0 * xf - x1

      case default                              ! 103, 104: the ghost value is a copy of the cell
        xc(b) % x(:,Ig,Jg,Kg) = x1

    end select

  end subroutine Ghost_Position


  !> Position of the first ghost cell of a face between blocks (101): the centre of the source cell.
  subroutine Connection_Ghost_Position ( domain, xc, r )
    use FUSS_Advanced_Types_m
    use FUSS_Parameters_m, only: guide
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    type(centres_t), intent(inout) :: xc(:)
    integer, intent(in) :: r
    ! Local
    integer :: b, bs, Ig, Jg, Kg

    b  = domain % bc(r) % b
    bs = domain % bc(r) % bs
    Ig = domain % bc(r) % i - guide(domain % bc(r) % f,1)
    Jg = domain % bc(r) % j - guide(domain % bc(r) % f,2)
    Kg = domain % bc(r) % k - guide(domain % bc(r) % f,3)
    xc(b) % x(:,Ig,Jg,Kg) = xc(bs) % x(:,domain % bc(r) % is, domain % bc(r) % js, domain % bc(r) % ks)

  end subroutine Connection_Ghost_Position


  !> Transposed face Jacobian of the mean of the metrics of two cells, (A_L + A_R)/2 with A = M^-1: its columns
  !> are the mean edges; the column across the face is x_R - x_L.
  function Mean_Edges ( M1, M2 ) result ( JT )
    use FUSS_Lib_Diffusive, only: Inverse3
    implicit none
    real(R8), intent(in) :: M1(3,3), M2(3,3)
    real(R8) :: JT(3,3)

    JT = 0.5d0 * ( Inverse3 ( M1 ) + Inverse3 ( M2 ) )

  end function Mean_Edges


  !> Metric J^-1 of a face of cell (i,j,k) of block b (the three differences of the stencil in the columns of J). The
  !> run stops when they do not span the space, |det J| <= 1e-12 |D_a x| |D_b x| |D_c x|: a degenerate face or cell,
  !> where the metric is not defined. A declared limit of diffusive-metric = stencil.
  function Face_Inverse ( JT, b, i, j, k, f, kind ) result ( Mf )
    use FUSS_Lib_Diffusive, only: Inverse3
    use FUSS_Mod_MPI, only: mpi_abort_all
    implicit none
    real(R8), intent(in) :: JT(3,3)
    integer, intent(in) :: b, i, j, k, f
    character(len=*), intent(in) :: kind
    real(R8) :: Mf(3,3)
    ! Local
    real(R8) :: det
    integer :: t

    det = JT(1,1) * ( JT(2,2)*JT(3,3) - JT(2,3)*JT(3,2) ) &
        - JT(1,2) * ( JT(2,1)*JT(3,3) - JT(2,3)*JT(3,1) ) &
        + JT(1,3) * ( JT(2,1)*JT(3,2) - JT(2,2)*JT(3,1) )
    if ( .not. ( abs(det) > 1d-12 * norm2(JT(:,1)) * norm2(JT(:,2)) * norm2(JT(:,3)) ) ) then
      write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] diffusive-metric = stencil: block ', b, ', cell', i, j, k, &
                                            ', face ', f, ' ('//kind//')'
      write(*,'(A)') ' [ERROR]   degenerate face or cell: the metric of the stencil is not defined there. The differences'
      write(*,'(A)') ' [ERROR]   of the cell centres it takes:'
      do t = 1, 3
        write(*,'(A,I0,A,3(1X,ES14.6))') ' [ERROR]   difference ', t, ':', JT(:,t)
      enddo
      call mpi_abort_all ( 'metric of diffusive-metric = stencil not defined (see above)' )
    end if
    Mf = Inverse3 ( JT )

  end function Face_Inverse


  !> Face metrics of the interior faces of block b: fm(d)%Mf(i,j,k) is the face between cell (i,j,k) and the next one
  !> along d.
  subroutine Interior_Face_Metrics ( b, blk, xc )
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: ndir
    implicit none
    integer, intent(in) :: b
    type(FUSS_block_type), intent(inout) :: blk
    type(centres_t), intent(in) :: xc
    ! Local
    integer :: i, j, k, n(3), d, t, i2, j2, k2
    integer :: e(3,3)
    real(R8) :: JT(3,3)

    n = blk % dim
    e = reshape ( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    do d = 1, 3
      do k = 1, n(3) - e(3,d) ; do j = 1, n(2) - e(2,d) ; do i = 1, n(1) - e(1,d)
        i2 = i + e(1,d) ; j2 = j + e(2,d) ; k2 = k + e(3,d)
        JT = Mean_Edges ( blk % M(i,j,k) % c, blk % M(i2,j2,k2) % c )
        do t = 1, 3
          if ( t == d .or. n(t) == 1 ) cycle
          if ( t == 3 .and. ndir == 2 ) cycle
          JT(:,t) = 0.25d0 * ( ( xc % x(:,i +e(1,t),j +e(2,t),k +e(3,t)) - xc % x(:,i -e(1,t),j -e(2,t),k -e(3,t)) ) &
                             + ( xc % x(:,i2+e(1,t),j2+e(2,t),k2+e(3,t)) - xc % x(:,i2-e(1,t),j2-e(2,t),k2-e(3,t)) ) )
        enddo
        blk % fm(d) % Mf(i,j,k) % c = Face_Inverse ( JT, b, i, j, k, 2*d, 'interior' )
      enddo ; enddo ; enddo
    enddo

  end subroutine Interior_Face_Metrics


  !> Face metric of the face of boundary record r between blocks, stored in the record (bc%Mf). The two cells in the
  !> order of growing index of block b (L, R) and their tangential neighbours: in block b, and for the ghost cell the
  !> cells of the source block that Fill_BC_Ghost_Connection takes for Tg(3:6).
  subroutine Connection_Face_Metric ( domain, xc, r )
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: ndir
    use FUSS_Parameters_m, only: guide
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    type(centres_t), intent(in) :: xc(:)
    integer, intent(in) :: r
    ! Local
    integer :: b, bs, Im, Jm, Km, Fm, dir, t, tb(2), q, n(3)
    integer :: e(3,3), sc(3,4)
    real(R8) :: JT(3,3), xg(3,4), dloc(3), dgh(3), Mc(3,3)

    b  = domain % bc(r) % b
    bs = domain % bc(r) % bs
    Im = domain % bc(r) % i ; Jm = domain % bc(r) % j ; Km = domain % bc(r) % k
    Fm = domain % bc(r) % f
    n  = domain % blk(b) % dim
    e  = reshape ( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    dir = (Fm + 1) / 2
    select case (dir)
      case (1) ; tb = [2, 3]
      case (2) ; tb = [1, 3]
      case default ; tb = [1, 2]
    end select

    ! Cells of the source block that hold Tg(3:6) (Fill_BC_Ghost_Connection), and their positions
    call Source_Neighbours ( domain % bc(r), sc )
    do q = 1, 4
      xg(:,q) = xc(bs) % x(:,sc(1,q),sc(2,q),sc(3,q))
    enddo

    ! Face Jacobian: the cells in the order L, R of block b (the ghost metric is that of the ghost cell, bc%Mg(1))
    Mc = domain % blk(b) % M(Im,Jm,Km) % c
    if ( mod(Fm,2) == 0 ) then
      JT = Mean_Edges ( Mc, domain % bc(r) % Mg(1) % c )
    else
      JT = Mean_Edges ( domain % bc(r) % Mg(1) % c, Mc )
    end if
    do q = 1, 2
      t = tb(q)
      if ( n(t) == 1 ) cycle
      if ( t == 3 .and. ndir == 2 ) cycle
      dloc = xc(b) % x(:,Im+e(1,t),Jm+e(2,t),Km+e(3,t)) - xc(b) % x(:,Im-e(1,t),Jm-e(2,t),Km-e(3,t))
      dgh  = xg(:,2*q) - xg(:,2*q-1)
      if ( mod(Fm,2) == 0 ) then
        JT(:,t) = 0.25d0 * ( dloc + dgh )
      else
        JT(:,t) = 0.25d0 * ( dgh + dloc )
      end if
    enddo
    domain % bc(r) % Mf = Face_Inverse ( JT, b, Im, Jm, Km, Fm, 'between blocks' )

  end subroutine Connection_Face_Metric


  !> Metric of the wall face of boundary record r (301-305): the gradient of Wall_Temperature takes the difference
  !> across the half cell, D_n = 2 s (Tw - T), and the differences of the cell values along the face, (T_hi - T_lo)/
  !> (hi - lo) for every tangential direction its loop uses; the metric takes the same differences of the positions,
  !> J(:,n) = 2 s (x_f - x_c) and J(:,d) = (x_hi - x_lo)/(hi - lo). Two variants as in the gradient: 1 with the
  !> ghost cells of the block (flux, ghosts = .true.), 2 one-sided inside the block (ghost fill, wall output,
  !> ghosts = .false.). The directions the loop skips keep the cell edge; a face with no tangential difference
  !> keeps the cell metric itself.
  subroutine Wall_Face_Metric ( domain, xc, r )
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: ndir
    use FUSS_Lib_Diffusive, only: Inverse3
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    type(centres_t), intent(in) :: xc(:)
    integer, intent(in) :: r
    ! Local
    integer :: b, Im, Jm, Km, Fm, dir, d, v, n(3), c(3), lo(3), hi(3), modfm2
    real(R8) :: JT(3,3), Ac(3,3), xf(3), nf(3)
    logical :: changed

    b  = domain % bc(r) % b
    Im = domain % bc(r) % i ; Jm = domain % bc(r) % j ; Km = domain % bc(r) % k
    Fm = domain % bc(r) % f
    c  = [Im, Jm, Km]
    n  = domain % blk(b) % dim
    dir = (Fm + 1) / 2
    modfm2 = 1 - 2 * mod(Fm,2)
    call Face_Centre_Normal ( domain % blk(b) % node, Im, Jm, Km, Fm, xf, nf )
    Ac = Inverse3 ( domain % blk(b) % M(Im,Jm,Km) % c )        ! columns: the edges of the cell
    do v = 1, 2
      JT = Ac
      JT(:,dir) = 2d0 * modfm2 * ( xf - xc(b) % x(:,Im,Jm,Km) )
      changed = .false.
      do d = 1, ndir
        if ( d == dir ) cycle
        lo = c ; hi = c
        lo(d) = lo(d) - 1
        hi(d) = hi(d) + 1
        if ( v == 2 ) then
          lo(d) = max ( lo(d), 1 )
          hi(d) = min ( hi(d), n(d) )
          if ( hi(d) == lo(d) ) cycle
        end if
        JT(:,d) = ( xc(b) % x(:,hi(1),hi(2),hi(3)) - xc(b) % x(:,lo(1),lo(2),lo(3)) ) / ( hi(d) - lo(d) )
        changed = .true.
      enddo
      if ( changed ) then
        domain % bc(r) % wM(:,:,v) = Face_Inverse ( JT, b, Im, Jm, Km, Fm, 'wall' )
      else
        domain % bc(r) % wM(:,:,v) = domain % blk(b) % M(Im,Jm,Km) % c
      end if
    enddo
    domain % bc(r) % wstencil = .true.

  end subroutine Wall_Face_Metric


  !> Indices in the source block of the cells of Tg(3:6) of a connection record, as in Fill_BC_Ghost_Connection.
  subroutine Source_Neighbours ( rec, sc )
    use FUSS_Advanced_Types_m, only: FUSS_bc_type
    implicit none
    type(FUSS_bc_type), intent(in) :: rec
    integer, intent(out) :: sc(3,4)
    ! Local
    integer :: Is, Js, Ks

    Is = rec % is ; Js = rec % js ; Ks = rec % ks
    select case ( rec % fs )
      case (1:2)
        sc(:,1) = [Is, Js - rec % d11, Ks - rec % d12]
        sc(:,2) = [Is, Js + rec % d11, Ks + rec % d12]
        sc(:,3) = [Is, Js - rec % d21, Ks - rec % d22]
        sc(:,4) = [Is, Js + rec % d21, Ks + rec % d22]
      case (3:4)
        sc(:,1) = [Is - rec % d11, Js, Ks - rec % d12]
        sc(:,2) = [Is + rec % d11, Js, Ks + rec % d12]
        sc(:,3) = [Is - rec % d21, Js, Ks - rec % d22]
        sc(:,4) = [Is + rec % d21, Js, Ks + rec % d22]
      case default
        sc(:,1) = [Is - rec % d11, Js - rec % d12, Ks]
        sc(:,2) = [Is + rec % d11, Js + rec % d12, Ks]
        sc(:,3) = [Is - rec % d21, Js - rec % d22, Ks]
        sc(:,4) = [Is + rec % d21, Js + rec % d22, Ks]
    end select

  end subroutine Source_Neighbours

end module FUSS_Lib_Stencil_Metrics
