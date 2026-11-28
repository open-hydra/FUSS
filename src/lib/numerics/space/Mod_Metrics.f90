module FUSS_Mod_Metrics
  use FUSS_Mod_MPI, only: fuss_abort, is_local_block
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

contains

  !> Geometry of every block and every ghost-cell interface, on every multigrid level (plan 14 D9).
  !>
  !> The interior cells (volume, metric tensor, cell lengths, face normals and areas) come from MORPH through Adapter_Morph
  !> (plan 09): a moving mesh refreshes them, vol_old and the swept volumes from the same kernel every step (Update_Mesh),
  !> so the discrete GCL holds by construction. Lib_Metrics gives the mesh type, the ghost-cell metrics and the set-up
  !> checks of the cells; Lib_Stencil_Metrics the face metrics of the diffusive fluxes, refreshed by Update_Mesh too.
  subroutine Setup_Metrics ( domain )
    use FUSS_Base_Types_m
    use FUSS_Advanced_Types_m
    use FUSS_Global_m
    use FUSS_Lib_Metrics, only: Check_Mesh_Type, Compute_Metric_Tensor, Stop_On_Bad_Metric, bad_metric, skew_factor_max
    use FUSS_Adapter_Morph, only: Adapter_Block_Metrics
    use FUSS_Lib_Stencil_Metrics, only: Setup_Face_Metrics
    use FUSS_Lib_Periodic, only: Setup_Periodic
    use FUSS_Mod_MPI, only: is_local_block, mpi_allreduce_max_r8, mpi_stop_if_any
    use FUSS_Mod_dt, only: Skew_Factor
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    logical  :: bad, ok
    integer  :: b, i, j, k, ci, cj, ck
    real(R8) :: fmax, volx
    character(len=256) :: message
    type(FUSS_vector_3D_type) :: N1, N2, N3, N4, N5, N6, N7, N8, dlx
    type(FUSS_tensor_3D_type) :: Mx

    ! Mesh type (ndir, delthe) from the nodes of every block (upstream a207291; plan 14 D-S7-2)
    call Check_Mesh_Type ( domain )

    ! A cell folded at one of its corners, or turned against its block, stops the run here
    call Check_Folded_Cells ( domain, bad )
    call mpi_stop_if_any ( bad, 'degenerate or folded cell in the mesh (see above)', '' )

    ! Interior cells: MORPH (plan 09). A block it rejects stops every rank together (plan 14 D22)
    bad = .false.
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call Adapter_Block_Metrics ( domain%blk(b), domain%nb, b, ok, message, ci, cj, ck )
      if ( .not. ok ) then
        write(*,'(A)')            ' [ERROR] mesh geometry is invalid'
        write(*,'(A,A)')          '         ', trim(message)
        write(*,'(A,4(I0,A))')    '         at block ', b, ', cell (', ci, ',', cj, ',', ck, ')'
        ! A whole block numbered left-handed: ours refuses it, upstream accepts it (plan 14 S9d, D-S9d-6, O-LEFT-HANDED)
        if ( Left_Handed_Block ( domain%blk(b) ) ) then
          write(*,'(A,I0,A)') '         Block ', b, ' is left-handed in (i, j, k): e_i . (e_j x e_k) < 0 at every corner of'
          write(*,'(A)')      '         every cell, a block numbered the other way, not an inverted cell. FUSS refuses a'
          write(*,'(A)')      '         left-handed block (its metrics need a positive volume): number the block right-handed'
          write(*,'(A)')      '         (reverse one index direction, or swap two), or mirror the mesh (one coordinate to its'
          write(*,'(A)')      '         negative) when the problem allows it.'
        endif
        bad = .true.
        exit
      endif
    enddo
    call mpi_stop_if_any ( bad, 'mesh geometry is invalid (see above)', '' )

    ! The size, sign and volume checks of Compute_Metric_Tensor on the interior cells, into scratch: MORPH's values stay
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp parallel do collapse(3) private(i, j, k, N1, N2, N3, N4, N5, N6, N7, N8, Mx, dlx, volx)
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
        call Compute_Metric_Tensor ( N1, N2, N3, N4, N5, N6, N7, N8, Mx, dlx, volx, loc = [b, i, j, k, 0] )
      enddo; enddo; enddo
      !$omp end parallel do
      if ( bad_metric ) exit
    enddo
    call Stop_On_Bad_Metric ( )

    ! A block too thin for the ghost cells and the face gradients that its boundary faces need: every rank stops together,
    ! before the ghost metrics are built (BC_Connect_Metrics would read node planes outside the block across)
    call Check_Thin_Blocks ( domain, bad )
    call mpi_stop_if_any ( bad, 'block too thin for its connections or its faces (see above)', '' )

    ! Periodic faces (201): the rigid motion of each pair of faces, before the ghost metrics that it moves
    call Setup_Periodic ( domain )

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

    call Setup_BC_Metrics ( domain )

    ! Face metrics of the diffusive fluxes from the positions of the cell centres (after the ghost metrics: the
    ! connections read bc%Mg(1))
    call Setup_Face_Metrics ( domain )

  end subroutine Setup_Metrics


  !> Ghost-cell metrics of every boundary record of the blocks of this rank, at set-up and on the per-step refresh of a
  !> moving mesh (Update_Mesh): they are functions of the block geometry, and on a moving mesh left at their t = 0 values
  !> the geometry the solver describes at the boundary is not the one it solves on.
  !>
  !> bc%Mg IS LOAD-BEARING. At a block connection (101) it is the neighbour's cell metric: the two-sided interface balance
  !> of diffusive-metric = mean (Connection_Face_Flux_Tint), the normal column of the stencil face metric and the width
  !> of the neighbour in the series conductivity (Lib_Stencil_Metrics, Connection_Face_Flux) all read it. Measured on
  !> test/numerics-features/conjugate_two_layer (dx2/dx1 = 4, exact solution piecewise linear) before plan 14: correct
  !> max |T - T_exact| = 4.5e-09 K, every bc%Mg halved 2.9e+01 K; on multimat_Plate (equal spacing at every interface)
  !> the same sabotage changes nothing. The Darcy flux of type 2 reads it per side (Mod_Darcy, Connection_Terms).
  !> bc%dlg and bc%volg are written here and read nowhere; they stay as intent(out) dummies so that nothing is left
  !> undefined.
  !>
  !> Every record of a block of another rank is skipped, at set-up and on the refresh alike: the ghost metrics of a block
  !> are read only by its rank, and the sign check of Compute_Metric_Tensor reads the metric of the first cell, which only
  !> that rank computed (MORPH's loop in Setup_Metrics skips the other blocks).
  subroutine Setup_BC_Metrics ( domain, refresh )
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc
    use FUSS_Lib_Metrics, only: BC_Connect_Metrics, BC_Symmetry_Metrics, BC_Extrapolate_Metrics, Stop_On_Bad_Metric
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    !> .true. on the per-step refresh of a moving mesh (Update_Mesh)
    logical, intent(in), optional :: refresh
    ! Local
    logical :: is_refresh
    integer :: i, g
    integer :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s

    is_refresh = .false.
    if ( present(refresh) ) is_refresh = refresh

    !$omp parallel private(g)
    ! The first ghost layer of every record (of the symmetry faces the first layer only: below)
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, domain % nbound
      Bm = domain % bc(i) % b
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
          if ( domain % bc(i) % per > 0 ) then
            ! periodic face (201): the nodes of the source block moved by the motion of the pair
            call BC_Connect_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                      Is, Js, Ks, Fs, domain % blk(Bs), d11s, d12s, d21s, d22s, &
                                      domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, &
                                      Rp = domain % per(domain % bc(i) % per) % R, &
                                      tp = domain % per(domain % bc(i) % per) % t )
          else
            call BC_Connect_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                      Is, Js, Ks, Fs, domain % blk(Bs), d11s, d12s, d21s, d22s, &
                                      domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg)
          end if
        case(300) ! symmetry
          call BC_Symmetry_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                     domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, 1 )
        case(103, 307) ! coupled faces (plan 08 decision 3): hydra overwrites
          ! this metric once, after FUSS's setup, with the MOSE neighbour's, and
          ! its interface temperature reads it; a refresh must keep hydra's (no
          ! receding block owns a coupled face, guard abl_receding_coupled)
          if ( .not. is_refresh ) &
            call BC_Extrapolate_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                          domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, Bm, &
                                          chimera = .false. )
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
    ! A degenerate or inverted ghost cell (Compute_Metric_Tensor wrote the first one): every rank stops together, on the
    ! refresh too (one reduction of a flag per step)
    call Stop_On_Bad_Metric ( )

  end subroutine Setup_BC_Metrics


  !> Advance the mesh one step and refresh all geometry. Called once per time
  !> step, NOT once per Runge-Kutta stage: the swept volumes describe the motion
  !> over the whole step, and re-applying the law per stage would move the mesh
  !> n_RK times as far.
  !>
  !> Must be called OUTSIDE any OpenMP parallel region -- the MORPH kernels open
  !> their own, exactly as Setup_Metrics does.
  !>
  !> `dt` here only decides WHERE the mesh is placed. The discrete GCL holds for
  !> whatever motion actually occurred, because the swept volumes are computed
  !> from the old and new node positions themselves, not from dt. So an imperfect
  !> dt estimate costs accuracy, never conservation.
  subroutine Update_Mesh ( domain, t, dt )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_mesh_motion
    use FUSS_Adapter_Morph,  only: Adapter_Update_Mesh
    use FUSS_Mod_Recession,  only: Compute_Node_Recession, Assert_Interface_Nodes
    use FUSS_Mod_GhostExchange, only: exchange_node_planes
    use FUSS_Mod_Periodic_Checks, only: Check_Periodic_Faces
    use FUSS_Lib_Stencil_Metrics, only: Setup_Face_Metrics
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    real(R8),               intent(in)    :: t, dt
    ! Local
    integer  :: b, ci, cj, ck, si, sj, sk
    logical  :: ok, nodeshift
    real(R8) :: gcl_rel, sweep_ratio
    character(len=256) :: message

    if ( .not. obj_mesh_motion%enabled ) return

    ! Surface recession (plan 08 section 2.3): one displacement per surface
    ! node, over every block that shares it, BEFORE any per-block law.
    nodeshift = ( trim(obj_mesh_motion%law) == 'nodeshift' )
    if ( nodeshift ) call Compute_Node_Recession ( domain, dt )

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle

      call Adapter_Update_Mesh ( domain%blk(b), domain%nb, b, t, dt, &
                                 gcl_rel, sweep_ratio, si, sj, sk, &
                                 ok, message, ci, cj, ck )

      if ( .not. ok ) then
        write(*,'(A)')         ' [ERROR] mesh update failed'
        write(*,'(A,A)')       '         ', trim(message)
        write(*,'(A,4(I0,A))') '         at block ', b, ', cell (', ci, ',', cj, ',', ck, ')'
        call fuss_abort ( 1 )
      endif

      ! Assert the discrete GCL every step. A violation means the geometry and
      ! the ALE state update disagree about how much volume was swept, which
      ! corrupts the solution silently -- so it is fatal, not a warning.
      if ( gcl_rel > obj_mesh_motion%gcl_tol ) then
        write(*,'(A)')          ' [ERROR] discrete Geometric Conservation Law violated'
        write(*,'(A,ES12.4,A,ES12.4)') '         relative residual ', gcl_rel, &
                                       ' exceeds tolerance ', obj_mesh_motion%gcl_tol
        write(*,'(A,I0)')       '         block ', b
        call fuss_abort ( 1 )
      endif

      ! Moving-mesh stability. The GCL says the geometry is self-consistent; it
      ! says nothing about the mesh having moved too far in one step. Once a
      ! cell sweeps out more than its own volume the remap's coefficient on that
      ! cell's own state goes negative and the update extrapolates. dt is chosen
      ! from conduction alone and carries no mesh-velocity limit, so nothing
      ! upstream prevents this.
      if ( sweep_ratio >= 1.0d0 ) then
        write(*,'(A)')          ' [ERROR] the mesh moved too far in one time step'
        write(*,'(A,ES12.4,A)') '         a cell swept out ', sweep_ratio, ' times its own volume;'
        write(*,'(A)')          '         the conservative remap is unstable at or above 1.'
        write(*,'(A,4(I0,A))')  '         worst at block ', b, ', cell (', si, ',', sj, ',', sk, ')'
        write(*,'(A)')          '         Reduce the mesh velocity or the time step (lower vnn).'
        call fuss_abort ( 1 )
      endif
    enddo

    ! A block connection between two ranks builds its ghost metric from the
    ! neighbour's node planes: bring them to this rank first (plan 11 stage 2).
    call exchange_node_planes ( domain )

    ! Periodic faces (201): the two faces of every pair must still be images by the pair's motion of the set-up, the
    ! motion the ghost metrics below are rebuilt with (plan 14 S9d, Mod_Periodic_Checks); after the exchange, which
    ! brings the other face of a pair across ranks
    call Check_Periodic_Faces ( domain )

    ! The shared nodes of every block interface must still coincide bitwise; after
    ! the exchange, so that a face shared with another rank's block is checked
    ! against that block's moved nodes (plan 11 stage 4). A periodic face (201)
    ! shares no node: Check_Periodic_Faces above checks it by images.
    if ( nodeshift ) call Assert_Interface_Nodes ( domain )

    ! The ghost-cell metrics and then the face metrics of the diffusive fluxes are functions of the geometry that has just
    ! moved: rebuilt in that order (Setup_Face_Metrics reads bc%Mg(1) of the connections and the node planes brought
    ! above). Computed once, the face metrics would keep the t = 0 geometry in every interior, wall and connection flux:
    ! the stale inverse-metric red control of the chimera merge (2026-10-17) moved ATC 2.2 / 2.3 by 54 / 265 K,
    ! invisible to the 1-vs-N gate. TEST ONLY skip-bc-metrics-refresh keeps both stale: the red control of the
    ! receding conjugate case.
    if ( .not. obj_mesh_motion%skip_bc_metrics_refresh ) then
      call Setup_BC_Metrics ( domain, refresh = .true. )
      call Setup_Face_Metrics ( domain )
    endif

  end subroutine Update_Mesh

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

  !> A left-handed block in (i, j, k): e_i . (e_j x e_k) < 0 at every corner of every cell, the edges from the corner as
  !> Check_Folded_Cells takes them (a corner with an edge of zero length, on an axis, is skipped). Ours' interior metrics
  !> (MORPH) need a positive volume and refuse such a block, which upstream accepts (plan 14 S9d, D-S9d-6, O-LEFT-HANDED);
  !> Setup_Metrics asks only after MORPH refused the block, to name the cause.
  function Left_Handed_Block ( blk ) result ( left )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    implicit none
    type(FUSS_block_type), intent(in) :: blk
    logical :: left
    ! Local
    integer  :: i, j, k, ca, cb, cc, ntest
    real(R8) :: ei(3), ej(3), ek(3)

    left = .false.
    ntest = 0
    associate ( node => blk % node )
      do k = 1, blk % dim(3); do j = 1, blk % dim(2); do i = 1, blk % dim(1)
        do cc = 0, 1; do cb = 0, 1; do ca = 0, 1
          ei = node(i,      j-1+cb, k-1+cc)%c - node(i-1,    j-1+cb, k-1+cc)%c
          ej = node(i-1+ca, j,      k-1+cc)%c - node(i-1+ca, j-1,    k-1+cc)%c
          ek = node(i-1+ca, j-1+cb, k     )%c - node(i-1+ca, j-1+cb, k-1   )%c
          if ( norm2 ( ei ) * norm2 ( ej ) * norm2 ( ek ) == 0d0 ) cycle
          ntest = ntest + 1
          if ( .not. ( ei(1) * ( ej(2)*ek(3) - ej(3)*ek(2) ) + ei(2) * ( ej(3)*ek(1) - ej(1)*ek(3) ) &
                     + ei(3) * ( ej(1)*ek(2) - ej(2)*ek(1) ) < 0d0 ) ) return
        enddo; enddo; enddo
      enddo; enddo; enddo
    end associate
    ! every tested corner is negative; a block with no corner tested (every one skipped) is not called left-handed
    left = ( ntest > 0 )

  end function Left_Handed_Block

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
