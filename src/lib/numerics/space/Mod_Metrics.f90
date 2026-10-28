module FUSS_Mod_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none

contains

  subroutine Setup_Metrics ( domain )
    use FUSS_Base_Types_m
    use FUSS_Advanced_Types_m
    use FUSS_Global_m
    use FUSS_Lib_Metrics
    use FUSS_Lib_Stencil_Metrics, only: Setup_Face_Metrics
    use FUSS_Mod_MPI, only: is_local_block, mpi_allreduce_max_r8, mpi_stop_if_any
    use FUSS_Mod_dt, only: Skew_Factor
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    logical :: bad
    ! Local
    integer :: b, i, j, k, g
    real(R8) :: fmax
    integer :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(FUSS_vector_3D_type) :: N1, N2, N3, N4, N5, N6, N7, N8

    call Check_Mesh_Type ( domain )

    ! A cell folded at one of its corners, or turned against its block, stops the run here
    call Check_Folded_Cells ( domain, bad )
    call mpi_stop_if_any ( bad, 'degenerate or folded cell in the mesh (see above)', '' )
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      ! Compute metric tensor and cell dimension across i,j,k. => block % M, & block % dl      
      !$omp do collapse(3) private(i, j, k, N1, N2, N3, N4, N5, N6, N7, N8)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        N1 % c = domain % blk(b) % node(i-1,j-1,k-1) % c
        N2 % c = domain % blk(b) % node(i-1,j-1,k  ) % c
        N3 % c = domain % blk(b) % node(i-1,j  ,k-1) % c
        N4 % c = domain % blk(b) % node(i-1,j  ,k  ) % c
        N5 % c = domain % blk(b) % node(i  ,j-1,k-1) % c
        N6 % c = domain % blk(b) % node(i  ,j-1,k  ) % c
        N7 % c = domain % blk(b) % node(i  ,j  ,k-1) % c
        N8 % c = domain % blk(b) % node(i  ,j  ,k  ) % c
        call Compute_Metric_Tensor ( N1, N2, N3, N4, N5, N6, N7, N8, domain % blk(b) % M(i,j,k), domain % blk(b) % dl(i,j,k), domain % blk(b) % vol(i,j,k), &
                                     loc = [b, i, j, k, 0] )
      enddo; enddo; enddo
      !$omp end parallel
      ! A degenerate or inverted cell (Compute_Metric_Tensor wrote the first one): no further block here; every rank
      ! stops together below
      if ( bad_metric ) exit

      ! Compute Normal & Area
      call Compute_Norm_Area ( domain % blk(b) )
    enddo
    call Stop_On_Bad_Metric ( )

    ! A block too thin for the ghost cells and the face gradients that its boundary faces need: every rank stops together,
    ! before the ghost metrics are built (BC_Connect_Metrics would read node planes outside the block across)
    call Check_Thin_Blocks ( domain, bad )
    call mpi_stop_if_any ( bad, 'block too thin for its connections or its faces (see above)', '' )

    ! Largest time-step skew factor of the cells of the finest level (1 unless 3-D cells are skewed in all pairs of
    ! directions), for the set-up summary and for the tests: the factor is 1 on every cell of a 2-D mesh
    if ( domain % mg_level == 1 ) then
      fmax = 1d0
      if ( ndir == 3 ) then
        do b = 1, domain % nb
          if (.not. is_local_block(b)) cycle
          do k = 1, domain % blk(b) % dim(3)
          do j = 1, domain % blk(b) % dim(2)
          do i = 1, domain % blk(b) % dim(1)
            fmax = max ( fmax, Skew_Factor ( domain % blk(b) % M(i,j,k) % c ) )
          enddo; enddo; enddo
        enddo
      end if
      call mpi_allreduce_max_r8 ( fmax, skew_factor_max )
    end if

    !$omp parallel private(g)
    ! Create the nodes for gc layers of ghost cell (of the symmetry faces the first layer only: below)
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, domain % nbound
      Bm = domain % bc(i) % b
      ! Ghost metrics land in blk(Bm) and bc(i), read only by the rank owning Bm (BC_Fluxes skips the others), and the
      ! sign check reads the metric of the first cell, computed by that rank only
      if (.not. is_local_block(Bm)) cycle
      Im = domain % bc(i) % i 
      Jm = domain % bc(i) % j 
      Km = domain % bc(i) % k 
      Fm = domain % bc(i) % f
      select case ( domain % bc(i) % type )
        case(101) ! block connection
          Bs = domain % bc(i) % bs
          Is = domain % bc(i) % is 
          Js = domain % bc(i) % js 
          Ks = domain % bc(i) % ks 
          Fs = domain % bc(i) % fs
          d11s = domain % bc(i) % d11
          d12s = domain % bc(i) % d12
          d21s = domain % bc(i) % d21
          d22s = domain % bc(i) % d22
          call BC_Connect_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                    Is, Js, Ks, Fs, domain % blk(Bs), d11s, d12s, d21s, d22s, &
                                    domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg)
        case(300) ! symmetry
          call BC_Symmetry_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                     domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, 1 )
        case default  ! chimera included
          call BC_Extrapolate_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                        domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, Bm, &
                                        chimera = ( domain % bc(i) % type == 102 ) )
      end select
    enddo
    !$omp end do
    ! The further layers of the symmetry faces copy the cells g - 1 behind the face, which in a direction with fewer
    ! cells than g are ghost cells of the opposite face: one layer at a time, after every face has written the previous
    ! one (the barrier at the end of each loop), so that no face reads a ghost cell that another one is writing
    do g = 2, gc
      !$omp do schedule (dynamic) private(i, Bm)
      do i = 1, domain % nbound
        if ( domain % bc(i) % type /= 300 ) cycle
        Bm = domain % bc(i) % b
        if (.not. is_local_block(Bm)) cycle
        call BC_Symmetry_Metrics ( domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k, &
                                   domain % bc(i) % f, domain % blk(Bm), &
                                   domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, g )
      enddo
      !$omp end do
    enddo
    !$omp end parallel
    ! A degenerate or inverted ghost cell (Compute_Metric_Tensor wrote the first one): every rank stops together
    call Stop_On_Bad_Metric ( )

    ! Face metrics of the diffusive fluxes from the positions of the cell centres
    call Setup_Face_Metrics ( domain )

  end subroutine Setup_Metrics


  !> A cell is the image of the unit cube by the trilinear map of its eight nodes, and a finite volume needs that map
  !> one-to-one: its Jacobian keeps the sign of the block in the whole cell. At a corner the Jacobian is
  !> J = e_i . (e_j x e_k), with e_i, e_j and e_k the three edges from that corner along i, j and k. The run stops at
  !> setup when, at a corner of a cell of a local block, s J <= 1e-12 |e_i| |e_j| |e_k| (Hadamard's bound, the measure of
  !> the determinant check of Compute_Metric_Tensor), with s the sign of the determinant of the mean edges of the first
  !> cell of the block, the cell from which Compute_Norm_Area takes the orientation of the normals; a degenerate first
  !> cell, which gives the block no orientation, stops the run there, as that determinant check would. A corner with an
  !> edge of zero length (nodes that coincide on an axis) is not checked. On a 2-D mesh, plane or a wedge whose nodes are
  !> off the axis, the corners decide for the whole cell; in 3-D they are a necessary condition (docs/theory/numerics.md,
  !> Grid Metrics).
  !> The volume of a cell, the integral of that Jacobian, has the sign of the block when the cell is not turned over, and
  !> Compute_Metric_Tensor stops a cell whose volume has not the sign of the determinant at its centre. Setup_Metrics calls
  !> this on every multigrid level, whose nodes are every other node of the level above: a stop there names the level.
  subroutine Check_Folded_Cells ( domain, bad )
    use FUSS_Advanced_Types_m
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    logical, intent(out) :: bad   !< A degenerate or folded cell was found on this rank (its message is written)
    ! Local
    integer  :: b, i, j, k, ca, cb, cc, n(3)
    real(R8) :: s, det, hm, ei(3), ej(3), ek(3), h, Jv, Am(3,3)
    character(len=32) :: lev

    bad = .false.
    lev = ':'
    if ( domain % mg_level > 1 ) write(lev,'(A,I0,A)') ' of multigrid level ', domain % mg_level, ':'

    do b = 1, domain % nb
      if ( .not. is_local_block(b) ) cycle
      n = domain % blk(b) % dim
      associate ( node => domain % blk(b) % node )
        ! Orientation of the block: the determinant of the mean edges of its first cell, as Compute_Metric_Tensor; a
        ! degenerate first cell gives none, and stops the run here as it would there
        Am(1,:) = node(1,0,0)%c - node(0,0,0)%c + node(1,0,1)%c - node(0,0,1)%c &
                + node(1,1,0)%c - node(0,1,0)%c + node(1,1,1)%c - node(0,1,1)%c
        Am(2,:) = node(0,1,0)%c - node(0,0,0)%c + node(0,1,1)%c - node(0,0,1)%c &
                + node(1,1,0)%c - node(1,0,0)%c + node(1,1,1)%c - node(1,0,1)%c
        Am(3,:) = node(0,0,1)%c - node(0,0,0)%c + node(0,1,1)%c - node(0,1,0)%c &
                + node(1,0,1)%c - node(1,0,0)%c + node(1,1,1)%c - node(1,1,0)%c
        Am  = 0.25d0 * Am
        det = Triple ( Am(1,:), Am(2,:), Am(3,:) )
        hm  = norm2 ( Am(1,:) ) * norm2 ( Am(2,:) ) * norm2 ( Am(3,:) )
        if ( .not. ( abs(det) > 1d-12 * hm ) ) then
          write(*,'(A,I0,A,ES12.4,A,ES12.4)') ' [ERROR] Mesh'//trim(lev)//' block ', b, &
            ', cell 1 1 1: degenerate, so the block has no orientation; det = ', det, ', |A1| |A2| |A3| = ', hm
          bad = .true.
          return
        end if
        s = sign ( 1d0, det )
        do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)
          do cc = 0, 1; do cb = 0, 1; do ca = 0, 1
            ei = node(i,      j-1+cb, k-1+cc)%c - node(i-1,    j-1+cb, k-1+cc)%c
            ej = node(i-1+ca, j,      k-1+cc)%c - node(i-1+ca, j-1,    k-1+cc)%c
            ek = node(i-1+ca, j-1+cb, k     )%c - node(i-1+ca, j-1+cb, k-1   )%c
            h = norm2 ( ei ) * norm2 ( ej ) * norm2 ( ek )
            if ( h == 0d0 ) cycle          ! an edge of zero length: nodes that coincide on an axis
            Jv = Triple ( ei, ej, ek )
            if ( s * Jv > 1d-12 * h ) cycle
            write(*,'(A,I0,A,3(1X,I0),A,3(1X,I0),A)') ' [ERROR] Mesh'//trim(lev)//' block ', b, ', cell', i, j, k, &
                                                     ', corner node', i-1+ca, j-1+cb, k-1+cc, &
                                                     ': the cell is folded or turned against its block'
            write(*,'(A,ES12.4,A,ES12.4,A,I0)') ' [ERROR]   e_i . (e_j x e_k) = ', Jv, ', |e_i| |e_j| |e_k| = ', h, &
                                               ', sign of the block (from its first cell): ', nint(s)
            write(*,'(A)') ' [ERROR]   At every corner of a cell the edges along i, j and k must keep the orientation of the block.'
            write(*,'(A,3(1X,ES14.6))') ' [ERROR]   corner node:', node(i-1+ca, j-1+cb, k-1+cc)%c
            bad = .true.
            return
          enddo; enddo; enddo
        enddo; enddo; enddo
      end associate
    enddo

  contains

    pure function Triple ( u, v, w ) result ( t )
      real(R8), intent(in) :: u(3), v(3), w(3)
      real(R8)             :: t
      t = u(1) * ( v(2)*w(3) - v(3)*w(2) ) + u(2) * ( v(3)*w(1) - v(1)*w(3) ) + u(3) * ( v(1)*w(2) - v(2)*w(1) )
    end function Triple

  end subroutine Check_Folded_Cells

  !> Blocks too thin for the ghost cells and the face gradients of their boundary faces. Two defects, one stop each, on
  !> every multigrid level (Setup_Metrics runs on each; the domain is that of the level):
  !>  - A connection (101): the gc ghost cells of its face are the first gc cells of the block across it, and their
  !>    metrics come from the first gc node planes of that block (BC_Connect_Metrics). With fewer cells across its face
  !>    (the same rule as MOSE's Pair_Connections, MOSE 2e31d8a), the second ghost cell is a ghost cell of that block,
  !>    which another record fills, and the node plane lies outside the node array: a build with bounds checking stops in
  !>    BC_Connect_Metrics; a release build reads past the array and runs on with a wrong or zero metric (a wrong wall
  !>    flux, or a stop on a cell with zero corner nodes). The record names the block, its face and the cell count.
  !>  - A direction d with one cell, other than the k direction of a 1-D or 2-D mesh (the field does not vary along it,
  !>    no flux reads its ghost cells; j is checked also in a 1-D mesh, whose walls on the j faces read it): the face gradient takes no difference along d (the interior
  !>    faces keep the cell edge in that column, the walls skip the one-sided difference of the wall output and of the
  !>    ghost fill: hi = lo), which loses the part of the gradient that the face normal takes from d. That part is
  !>    M(d,:) . n_e, zero where the lines of d are orthogonal to those of every other direction e (Cartesian blocks:
  !>    thin3d, one cell in k), and the stop is for cells where cos(M(d,:), M(e,:)) exceeds COS_MAX. A sheared block
  !>    one cell across gave a wrong wall flux (qw error 0.16 to 0.30 of k |g|) with no message and no out-of-range read.
  subroutine Check_Thin_Blocks ( domain, bad )
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: ndir, gc
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    logical, intent(out) :: bad   !< A thin block was found on this rank (its message is written)
    ! Local
    character(len=3), parameter :: DIRS = 'ijk'
    real(R8), parameter :: COS_MAX = 1d-3   ! cosine of the angle between the lines of two directions: 0 for orthogonal cells
    integer  :: r, b, bs, f, fs, n, nb, d, e, i, j, k, ie(3), iw(3)
    real(R8) :: c, cmax, vd(3), ve(3)
    logical, allocatable :: seen(:,:)
    character(len=32) :: lev

    bad = .false.
    nb  = domain % nb
    lev = ''
    if ( domain % mg_level > 1 ) write(lev,'(A,I0)') ' of multigrid level ', domain % mg_level

    ! Connections: the block across has at least gc cells across its face (one message per block face)
    allocate ( seen(nb,6) )
    seen = .false.
    do r = 1, domain % nbound
      if ( domain % bc(r) % type /= 101 ) cycle
      b  = domain % bc(r) % b
      f  = domain % bc(r) % f
      bs = domain % bc(r) % bs
      fs = domain % bc(r) % fs
      if ( .not. is_local_block(b) ) cycle
      if ( bs < 1 .or. bs > nb .or. fs < 1 .or. fs > 6 ) cycle      ! Pair_Connections stops on these
      n = domain % blk(bs) % dim((fs + 1) / 2)
      if ( n >= gc .or. seen(b,f) ) cycle
      seen(b,f) = .true.
      write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') ' [ERROR] Mesh'//trim(lev)//': block ', b, ', face ', f, &
        ' (connection record ', r, '): block ', bs, ' across it has ', n, ' cell(s) across its face ', fs, &
        ', fewer than the ', gc, ' that the ghost cells of a connection need'
      bad = .true.
    enddo
    deallocate ( seen )

    ! One cell in an index direction of the mesh, with lines of another direction that are not orthogonal to it
    do b = 1, nb
      if ( .not. is_local_block(b) ) cycle
      do d = 1, 3
        if ( d == 3 .and. ndir < 3 ) cycle      ! the k direction of a 1-D or 2-D mesh: the field does not vary along it
        if ( domain % blk(b) % dim(d) /= 1 ) cycle
        do e = 1, 3
          if ( e == d ) cycle
          cmax = 0d0 ; iw = 1
          do k = 1, domain % blk(b) % dim(3)
          do j = 1, domain % blk(b) % dim(2)
          do i = 1, domain % blk(b) % dim(1)
            vd = domain % blk(b) % M(i,j,k) % c(d,:)
            ve = domain % blk(b) % M(i,j,k) % c(e,:)
            c  = norm2 ( vd ) * norm2 ( ve )
            if ( c > 0d0 ) c = abs ( dot_product ( vd, ve ) ) / c
            if ( c > cmax ) then
              cmax = c ; iw = [i, j, k]
            end if
          enddo ; enddo ; enddo
          if ( cmax <= COS_MAX ) cycle
          ie = iw
          write(*,'(A,I0,A,I0,A,I0,A,ES10.3,A)') ' [ERROR] Mesh'//trim(lev)//': block ', b, ' has one cell in direction ', d, &
            ' and its lines in direction ', e, ' are not orthogonal to it (cos up to ', cmax, ', cell '// &
            trim(Cell_Text(ie))//'):'
          write(*,'(A)') ' [ERROR]   the face gradient takes no difference along a direction with one cell, so the heat flux'// &
            ' of the faces of this block would miss that part'
          write(*,'(A,I0,A,I0,A)') ' [ERROR]   give block ', b, ' at least 2 cells in direction ', d, &
            ', or make its direction '//DIRS(d:d)//' orthogonal to the others'
          bad = .true.
          exit
        enddo
      enddo
    enddo

  contains

    function Cell_Text ( c3 ) result ( t )
      integer, intent(in) :: c3(3)
      character(len=40) :: t
      write(t,'(I0,1X,I0,1X,I0)') c3
    end function Cell_Text

  end subroutine Check_Thin_Blocks


end module FUSS_Mod_Metrics