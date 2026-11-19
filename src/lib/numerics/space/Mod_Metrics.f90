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
          call BC_Connect_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                    Is, Js, Ks, Fs, domain % blk(Bs), d11s, d12s, d21s, d22s, &
                                    domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg)
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

    ! The shared nodes of every block interface must still coincide bitwise; after
    ! the exchange, so that a face shared with another rank's block is checked
    ! against that block's moved nodes (plan 11 stage 4).
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

end module FUSS_Mod_Metrics
