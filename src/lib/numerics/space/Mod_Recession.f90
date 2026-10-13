!> The ablating wall's per-step and per-output bookkeeping (plan 08 sections
!> 2.1, 2.3 R3, 3.2): the env values of every 306 face at the step's time, the
!> restart state of the faces, and the two files the faces write.
!>
!> ENV VALUES. The env table of a 306 face (t p_w rhoeUeCH0 h_r q_rad chem
!> hconv T_edge) is evaluated ONCE per step, after the clock has advanced to
!> t^(n+1) and before the type-2 pressure solve, which reads the face's p_w as
!> its Dirichlet value: a refresh inside the BC routine (stage time) would hand
!> the backward-Euler p^(n+1) the previous step's p_w, and p-wall on step 1.
!> Every RK stage runs at t^(n+1), so the stored vector is what a stage-time
!> evaluation would give. Also called once at setup, after the clock is set.
!>
!> RESTART STATE. Besides the node positions the solution file carries, a
!> face's balance has state: the recession so far, the wall temperature (the
!> Newton start), the recession rate (it moves the mesh one step later), the
!> blowing correction (the Picard start) and the B'_pg of the last chemistry-
!> on solve (the fractional-chem rule). OUTPUT/<prefix>bc306-state.dat holds
!> them, rewritten at every solution output with 17 significant digits (an
!> exact binary64 round trip: the converged T_w depends on its start at the
!> seb-tol level) and read back when newrun = false.
!>
!> SEB FILE. OUTPUT/<prefix>seb.txt: one row per 306 face per solution output
!> with every term of the last converged balance (plan 08 section 3.2).
!>
!> RECESSION (law = nodeshift, plan 08 section 2.3). A block that owns 306
!> faces recedes from that block face: recede_dir(b) is the grid direction of
!> its lines (1..3) and recede_side(b) the face (1 = low index, 2 = high),
!> derived from the 306 face id; 0 = the block does not move. Every step,
!> before the per-block motion law, Compute_Node_Recession turns the faces'
!> recession rates into ONE displacement per surface node: the area-weighted
!> mean rate of every 306 face adjacent to the node, in every block, so that a
!> node shared by two blocks moves identically in both (a per-block law fed
!> per-face rates would tear the mesh). The node moves along its grid line,
!> its normal component -sdot dt. Assert_Interface_Nodes then checks after the
!> move that the four nodes of every 101 face agree bitwise on both sides.
!> The surface-node map is built once, at setup, by exact coordinate equality
!> of the block-boundary nodes (they come from one MESH/ file).
module FUSS_Mod_Recession
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Ablation_Update_Env, Setup_Ablation, Write_Ablation_Outputs
  public :: Setup_Recession, Compute_Node_Recession, Assert_Interface_Nodes, Write_Interface_Flux
  public :: recede_dir, recede_side

  logical, save :: seb_file_started = .false.
  logical, save :: if_file_started  = .false.

  integer, allocatable, save :: recede_dir(:), recede_side(:)

  !> One distinct surface node: its coordinates at setup, every (b, i, j, k)
  !> incarnation, and the bc indices of the 306 faces that contain it.
  type :: surf_node_t
    real(R8) :: x(3) = 0d0
    integer  :: ninc = 0, nface = 0
    integer, allocatable :: inc(:,:)      ! (4, ninc): b, i, j, k
    integer, allocatable :: face(:)       ! bc indices
  end type surf_node_t
  type(surf_node_t), allocatable, save :: smap(:)

  real(R8), save :: min_thick = 0d0             ! resolved min-thickness [m]
  real(R8), save :: min_thick_restart = -1d0    ! the value bc306-state.dat carried, < 0 = none
  integer,  save :: nstep_mesh = 0              ! nodeshift steps so far (the TEST ONLY perturbation fires at 2)

contains

  !> Env values of every 306 face at domain%time (serial; outside the
  !> parallel region).
  subroutine Ablation_Update_Env ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_ablation
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: i
    real(R8) :: v(7)

    if ( .not. obj_ablation%enabled ) return
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 306 ) cycle
      v = domain%bc(i)%env%eval( domain%time )
      domain%bc(i)%pw        = v(1)
      domain%bc(i)%rhoeUeCH0 = v(2)
      domain%bc(i)%hr        = v(3)
      domain%bc(i)%qrad_env  = v(4)
      domain%bc(i)%chem      = v(5)
      domain%bc(i)%hconv_env = v(6)
      domain%bc(i)%Tedge     = v(7)
    enddo

  end subroutine Ablation_Update_Env


  !> Once at setup, after Setup_BC: on a restart, the faces' state from
  !> OUTPUT/<prefix>bc306-state.dat (every 306 face must have its row, and
  !> every row its face).
  subroutine Setup_Ablation ( domain, newrun )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_ablation
    use FUSS_Global_m,       only: FUSS_phase_prefix
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    logical,                intent(in)    :: newrun
    ! Local
    character(len=512) :: file, line
    integer  :: u, ios, ln, i, b, ii, jj, kk, f, nrows
    real(R8) :: v(5)
    logical, allocatable :: seen(:)

    if ( .not. obj_ablation%enabled ) return
    if ( newrun ) return

    file = 'OUTPUT/'//trim(FUSS_phase_prefix)//'bc306-state.dat'
    open( newunit=u, file=trim(file), status='old', action='read', iostat=ios )
    if ( ios /= 0 ) then
      write(*,'(A)') ' [ERROR] restart with ablating faces (306): '//trim(file)//' is missing'
      write(*,'(A)') '         it holds each face''s recession, wall temperature, recession rate and Picard start'
      error stop 1
    endif
    allocate( seen(domain%nbound) )
    seen = .false.
    ln = 0 ; nrows = 0
    do
      read( u, '(A)', iostat=ios ) line
      if ( ios /= 0 ) exit
      ln = ln + 1
      line = adjustl(line)
      if ( line(1:15) == '# min-thickness' ) then
        read( line(16:), *, iostat=ios ) min_thick_restart
        if ( ios /= 0 ) call Fail( 'expected "# min-thickness <value>"' )
        cycle
      endif
      if ( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
      read( line, *, iostat=ios ) b, ii, jj, kk, f, v
      if ( ios /= 0 ) call Fail( 'expected "b i j k f s_total Tw_prev sdot beta Bpg_lag"' )
      nrows = nrows + 1
      do i = 1, domain%nbound
        if ( domain%bc(i)%type /= 306 ) cycle
        if ( domain%bc(i)%b == b .and. domain%bc(i)%i == ii .and. domain%bc(i)%j == jj .and. &
             domain%bc(i)%k == kk .and. domain%bc(i)%f == f ) exit
      enddo
      if ( i > domain%nbound ) call Fail( 'no 306 face of bc.txt matches this row' )
      if ( seen(i) ) call Fail( 'a second row for the same face' )
      seen(i) = .true.
      domain%bc(i)%s_total = v(1)
      domain%bc(i)%Tw_prev = v(2)
      domain%bc(i)%sdot    = v(3)
      domain%bc(i)%beta    = v(4)
      domain%bc(i)%Bpg_lag = v(5)
      domain%bc(i)%Tw_face = v(2)
    enddo
    close( u )
    do i = 1, domain%nbound
      if ( domain%bc(i)%type == 306 .and. .not. seen(i) ) then
        write(*,'(A,I0,A)') ' [ERROR] '//trim(file)//': no row for the 306 face of bc.txt record ', i, &
          ' (written by a run with another bc.txt?)'
        error stop 1
      endif
    enddo

  contains

    subroutine Fail ( reason )
      character(len=*), intent(in) :: reason
      write(*,'(A,I0,A)') ' [ERROR] '//trim(file)//' line ', ln, ': '//reason
      error stop 1
    end subroutine Fail

  end subroutine Setup_Ablation


  !> At every solution output (root): the seb.txt rows and the restart state.
  subroutine Write_Ablation_Outputs ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_ablation, obj_sim_param
    use FUSS_Global_m,       only: FUSS_phase_prefix
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    character(len=512) :: file
    integer :: u, ios, i

    if ( .not. obj_ablation%enabled ) return

    ! seb.txt: appended, header on a new run's first write
    file = 'OUTPUT/'//trim(FUSS_phase_prefix)//'seb.txt'
    if ( .not. seb_file_started ) then
      if ( obj_sim_param%newrun ) then
        open( newunit=u, file=trim(file), status='replace', action='write', iostat=ios )
        if ( ios == 0 ) write(u,'(A)') '# time b i j k f T_w T_c beta Bpg Bch h_w rhoeUeCH0 h_r mdot_pg mdot_ch iter '// &
          'nclamp_p nclamp_Bg nclamp_T nclamp_rc  cond conv adv_pg adv_ch alpha_qrad eps_Tinf4 m_eps_Tw4 F   (SI)'
      else
        open( newunit=u, file=trim(file), status='unknown', position='append', action='write', iostat=ios )
      endif
      seb_file_started = .true.
    else
      open( newunit=u, file=trim(file), status='old', position='append', action='write', iostat=ios )
    endif
    if ( ios /= 0 ) then
      write(*,'(A)') ' [ERROR] ablating wall: cannot open '//trim(file)
      error stop 1
    endif
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 306 ) cycle
      associate ( bc => domain%bc(i) )
        write(u,'(ES24.16,5I8,10ES24.16,5I8,8ES24.16)') domain%time, bc%b, bc%i, bc%j, bc%k, bc%f, &
          bc%Tw_prev, domain%blk(bc%b)%T(bc%i,bc%j,bc%k), bc%beta, bc%Bpg, bc%Bch, bc%hw, bc%rhoeUeCH0, bc%hr, &
          bc%mdot_pg, bc%mdot_ch, bc%seb_iter, bc%nclamp, bc%seb
      end associate
    enddo
    close(u)

    ! bc306-state.dat: the restart state, rewritten at every output
    file = 'OUTPUT/'//trim(FUSS_phase_prefix)//'bc306-state.dat'
    open( newunit=u, file=trim(file), status='replace', action='write', iostat=ios )
    if ( ios /= 0 ) then
      write(*,'(A)') ' [ERROR] ablating wall: cannot open '//trim(file)
      error stop 1
    endif
    write(u,'(A)') '# b i j k f s_total Tw_prev sdot beta Bpg_lag'
    write(u,'(A,ES24.16)') '# min-thickness ', min_thick
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 306 ) cycle
      associate ( bc => domain%bc(i) )
        write(u,'(5I8,5ES24.16)') bc%b, bc%i, bc%j, bc%k, bc%f, bc%s_total, bc%Tw_prev, bc%sdot, bc%beta, bc%Bpg_lag
      end associate
    enddo
    close(u)

  end subroutine Write_Ablation_Outputs


  !> Once at setup, after Setup_BC and the bc-level checks, law = nodeshift
  !> only (serial: MPI is refused with 306). The receding face of every block,
  !> the surface-node map, its rules (plan 08 section 2.3) -- (a) the 306 patch
  !> is closed across 101 interfaces, (b) blocks sharing a receding line give
  !> it the same length -- and the resolved min-thickness.
  subroutine Setup_Recession ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_mesh_motion, obj_ablation
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, b, f, n, nb, nbn, nsf, q, r, lo, ii, jj, kk, nd(3,4), cnt
    integer, allocatable :: bid(:,:), perm(:), sfi(:), sperm(:)
    real(R8), allocatable :: bx(:,:), sx(:,:)
    real(R8) :: L1, L2, hmin
    logical :: bad, found

    if ( trim(obj_mesh_motion%law) /= 'nodeshift' ) return
    nb = domain%nb
    if ( allocated(recede_dir) ) deallocate( recede_dir, recede_side )
    allocate( recede_dir(nb), recede_side(nb) )
    recede_dir = 0 ; recede_side = 0
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 306 ) cycle
      b = domain%bc(i)%b ; f = domain%bc(i)%f
      recede_dir(b)  = ( f + 1 ) / 2
      recede_side(b) = 2 - mod( f, 2 )
    enddo

    ! every block-boundary node, sorted by its coordinates
    nbn = 0
    do b = 1, nb
      nbn = nbn + n_boundary_nodes( domain%blk(b)%dim )
    enddo
    allocate( bx(3,nbn), bid(4,nbn), perm(nbn) )
    n = 0
    do b = 1, nb
      do kk = 0, domain%blk(b)%dim(3) ; do jj = 0, domain%blk(b)%dim(2) ; do ii = 0, domain%blk(b)%dim(1)
        if ( ii == 0 .or. ii == domain%blk(b)%dim(1) .or. jj == 0 .or. jj == domain%blk(b)%dim(2) .or. &
             kk == 0 .or. kk == domain%blk(b)%dim(3) ) then
          n = n + 1
          bx(:,n)  = domain%blk(b)%node(ii,jj,kk)%c
          bid(:,n) = [ b, ii, jj, kk ]
        endif
      enddo ; enddo ; enddo
    enddo
    call sort_coords( bx, perm )

    ! the (306 face, node) pairs, sorted by the node's coordinates: equal
    ! coordinates are one surface node
    nsf = 4 * obj_ablation%n_faces
    allocate( sx(3,nsf), sfi(nsf), sperm(nsf) )
    n = 0
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 306 ) cycle
      call face_nodes( domain%bc(i)%f, domain%bc(i)%i, domain%bc(i)%j, domain%bc(i)%k, nd )
      do q = 1, 4
        n = n + 1
        sx(:,n) = domain%blk(domain%bc(i)%b)%node(nd(1,q),nd(2,q),nd(3,q))%c
        sfi(n) = i
      enddo
    enddo
    nsf = n
    call sort_coords( sx(:,1:nsf), sperm(1:nsf) )
    cnt = 0
    do n = 1, nsf
      if ( n == 1 ) then
        cnt = 1
      elseif ( any( sx(:,sperm(n)) /= sx(:,sperm(n-1)) ) ) then
        cnt = cnt + 1
      endif
    enddo
    if ( allocated(smap) ) deallocate( smap )
    allocate( smap(cnt) )
    cnt = 0
    do n = 1, nsf
      if ( n == 1 ) then
        cnt = 1
      elseif ( any( sx(:,sperm(n)) /= sx(:,sperm(n-1)) ) ) then
        cnt = cnt + 1
      endif
      smap(cnt)%x = sx(:,sperm(n))
      call push_face( smap(cnt), sfi(sperm(n)) )
    enddo

    ! incarnations: the boundary nodes with the same coordinates
    do n = 1, size(smap)
      lo = first_not_below( bx, perm, smap(n)%x )
      smap(n)%ninc = 0
      do r = lo, nbn
        if ( any( bx(:,perm(r)) /= smap(n)%x ) ) exit
        call push_inc( smap(n), bid(:,perm(r)) )
      enddo
    enddo

    ! rule (a): every incarnation sits on a 306 face of its own block
    bad = .false.
    do n = 1, size(smap)
      do r = 1, smap(n)%ninc
        found = .false.
        do q = 1, smap(n)%nface
          if ( domain%bc(smap(n)%face(q))%b == smap(n)%inc(1,r) ) found = .true.
        enddo
        if ( .not. found ) then
          write(*,'(A,I0,A,3ES14.6,A)') ' [ERROR] law = nodeshift: block ', smap(n)%inc(1,r), &
            ' shares the node (', smap(n)%x, ') of a 306 face but has no 306 face containing it:'
          write(*,'(A)') '         the 306 patch must be closed across 101 interfaces (plan 08 section 2.3 (a)).'
          bad = .true.
        endif
      enddo
      if ( bad ) exit
    enddo
    ! rule (b): a node shared by receding blocks ends lines of the same length
    if ( .not. bad ) then
      do n = 1, size(smap)
        if ( smap(n)%ninc < 2 ) cycle
        L1 = line_length( domain, smap(n)%inc(:,1) )
        do r = 2, smap(n)%ninc
          L2 = line_length( domain, smap(n)%inc(:,r) )
          if ( L2 /= L1 ) then
            write(*,'(A,I0,A,I0,A,3ES14.6,A)') ' [ERROR] law = nodeshift: blocks ', smap(n)%inc(1,1), ' and ', &
              smap(n)%inc(1,r), ' share the receding node (', smap(n)%x, ') on lines of different length'
            write(*,'(A,2ES24.16,A)') '         (', L1, L2, ' m): the interface nodes would move apart '// &
              '(plan 08 section 2.3 (b)).'
            bad = .true.
            exit
          endif
        enddo
        if ( bad ) exit
      enddo
    endif
    if ( bad ) error stop 1

    ! min-thickness: the key, else the value a restart carried, else 10 x the
    ! smallest cell on the receding lines of this (new) run's mesh, capped at
    ! half the shortest receding line (10 cells of a coarse line can be longer
    ! than the line: 4 cells of 5 mm would give 50 mm on a 20 mm line)
    if ( obj_ablation%min_thickness > 0.0d0 ) then
      min_thick = obj_ablation%min_thickness
    elseif ( min_thick_restart >= 0.0d0 ) then
      min_thick = min_thick_restart
    else
      hmin = huge(1.0d0) ; L1 = huge(1.0d0)
      do n = 1, size(smap)
        do r = 1, smap(n)%ninc
          hmin = min( hmin, line_min_cell( domain, smap(n)%inc(:,r) ) )
          L1   = min( L1, line_length( domain, smap(n)%inc(:,r) ) )
        enddo
      enddo
      min_thick = min( 10.0d0 * hmin, 0.5d0 * L1 )
    endif
    write(*,'(A,I0,A,I0,A,ES12.4,A)') '   Receding blocks (nodeshift)   ', count( recede_dir > 0 ), &
      ', surface nodes ', size(smap), ', min-thickness ', min_thick, ' m'

  end subroutine Setup_Recession


  !> Once per step, from Update_Mesh, before the per-block law: the recession
  !> so far, then every surface node's displacement into blk%node_disp of
  !> every incarnation, and the consumed-block stop. Serial (Update_Mesh runs
  !> outside the parallel region).
  subroutine Compute_Node_Recession ( domain, dt )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_ablation
    use FUSS_Lib_BC_Fluxes,  only: Face_Index
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    real(R8),               intent(in)    :: dt
    ! Local
    integer  :: i, b, n, q, r, dr, fi, fj, fk, inc(4), nb_idx(3)
    real(R8) :: A, sum_a, sum_as, sdn, nrm(3), e(3), xs(3), cosang, disp(3), L

    ! the recession rates (TEST ONLY override) and the recession so far
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 306 ) cycle
      if ( obj_ablation%sdot_override >= 0.0d0 ) domain%bc(i)%sdot = obj_ablation%sdot_override
      domain%bc(i)%s_total = domain%bc(i)%s_total + domain%bc(i)%sdot * dt
    enddo

    do b = 1, domain%nb
      if ( recede_dir(b) == 0 ) cycle
      do r = 0, domain%blk(b)%dim(3) ; do q = 0, domain%blk(b)%dim(2) ; do n = 0, domain%blk(b)%dim(1)
        domain%blk(b)%node_disp(n,q,r)%c = 0.0d0
      enddo ; enddo ; enddo
    enddo

    do n = 1, size(smap)
      sum_a = 0.0d0 ; sum_as = 0.0d0 ; nrm = 0.0d0
      do q = 1, smap(n)%nface
        associate ( bc => domain%bc(smap(n)%face(q)) )
          call Face_Index ( bc%f, dr, bc%i, bc%j, bc%k, fi, fj, fk )
          A = domain%blk(bc%b)%dir(dr)%f(fi,fj,fk)%A
          sum_a  = sum_a  + A
          sum_as = sum_as + A * bc%sdot
          nrm    = nrm + A * domain%blk(bc%b)%dir(dr)%f(fi,fj,fk)%N * real( 1 - 2 * mod(bc%f,2), R8 )
        end associate
      enddo
      sdn = sum_as / sum_a
      if ( sdn == 0.0d0 ) cycle               ! no motion: the nodes stay bitwise
      nrm = nrm / norm2( nrm )
      ! the line direction, from the node toward its neighbour on the line
      inc = smap(n)%inc(:,1)
      b = inc(1)
      xs = domain%blk(b)%node(inc(2),inc(3),inc(4))%c
      nb_idx = inc(2:4)
      if ( recede_side(b) == 2 ) then
        nb_idx(recede_dir(b)) = nb_idx(recede_dir(b)) - 1
      else
        nb_idx(recede_dir(b)) = nb_idx(recede_dir(b)) + 1
      endif
      e = domain%blk(b)%node(nb_idx(1),nb_idx(2),nb_idx(3))%c - xs
      e = e / norm2( e )
      cosang = abs( dot_product( nrm, e ) )
      disp = ( sdn * dt / cosang ) * e
      do r = 1, smap(n)%ninc
        inc = smap(n)%inc(:,r)
        domain%blk(inc(1))%node_disp(inc(2),inc(3),inc(4))%c = disp
        L = line_length( domain, inc )
        if ( L - norm2( disp ) <= min_thick ) then
          write(*,'(A,I0,A,3(I0,A))') ' [ERROR] ablating block consumed: block ', inc(1), ', surface node (', &
            inc(2), ',', inc(3), ',', inc(4), ')'
          write(*,'(A,ES12.4,A,ES12.4,A)') '         the receding line would keep ', L - norm2( disp ), &
            ' m, at or below min-thickness ', min_thick, ' m'
          error stop 2
        endif
      enddo
    enddo

  end subroutine Compute_Node_Recession


  !> After the move, before the ghost metrics are rebuilt: the four nodes of
  !> every 101 face agree bitwise with its partner's (as sorted coordinate
  !> triples), else error stop 2. TEST ONLY: perturb-interface-node moves one
  !> shared node of the lowest receding block by 1e-12 m in x at step 2.
  subroutine Assert_Interface_Nodes ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_ablation
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: i, q, br, nd(3,4), ns(3,4)
    real(R8) :: xa(3,4), xb(3,4)
    integer  :: pa(4), pb(4)

    nstep_mesh = nstep_mesh + 1
    if ( obj_ablation%perturb_interface_node .and. nstep_mesh == 2 ) then
      br = 0
      do i = 1, domain%nb
        if ( recede_dir(i) > 0 ) then
          br = i
          exit
        endif
      enddo
      do i = 1, domain%nbound
        if ( domain%bc(i)%type /= 101 .or. domain%bc(i)%b /= br ) cycle
        call face_nodes( domain%bc(i)%f, domain%bc(i)%i, domain%bc(i)%j, domain%bc(i)%k, nd )
        domain%blk(br)%node(nd(1,1),nd(2,1),nd(3,1))%c(1) = domain%blk(br)%node(nd(1,1),nd(2,1),nd(3,1))%c(1) + 1.0d-12
        write(*,'(A,I0,A,3(I0,A))') '   TEST ONLY perturb-interface-node: block ', br, ', node (', &
          nd(1,1), ',', nd(2,1), ',', nd(3,1), ') moved by 1e-12 m in x'
        exit
      enddo
    endif

    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 101 ) cycle
      associate ( bc => domain%bc(i) )
        call face_nodes( bc%f,  bc%i,  bc%j,  bc%k,  nd )
        call face_nodes( bc%fs, bc%is, bc%js, bc%ks, ns )
        do q = 1, 4
          xa(:,q) = domain%blk(bc%b )%node(nd(1,q),nd(2,q),nd(3,q))%c
          xb(:,q) = domain%blk(bc%bs)%node(ns(1,q),ns(2,q),ns(3,q))%c
        enddo
        call sort_coords( xa, pa )
        call sort_coords( xb, pb )
        do q = 1, 4
          if ( any( xa(:,pa(q)) /= xb(:,pb(q)) ) ) then
            write(*,'(A,I0,A,I0,A)') ' [ERROR] the mesh tore at a block interface: blocks ', bc%b, ' and ', bc%bs, &
              ' no longer share a node of their 101 face'
            write(*,'(A,3ES24.16)') '         block side: ', xa(:,pa(q))
            write(*,'(A,3ES24.16)') '         other side: ', xb(:,pb(q))
            write(*,'(A,4(I0,A))')  '         bc.txt record ', i, ', cell (', bc%i, ',', bc%j, ',', bc%k, ')'
            error stop 2
          endif
        enddo
      end associate
    enddo

  end subroutine Assert_Interface_Nodes


  !> TEST ONLY (dump-interface-flux): at every solution output, one row per
  !> 101 face with the flux BC_Connection put into the residual (and, type 2,
  !> the Darcy mass and enthalpy fluxes), for the pairing check of the
  !> receding conjugate case.
  subroutine Write_Interface_Flux ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_ablation, obj_sim_param
    use FUSS_Global_m,       only: FUSS_phase_prefix
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    character(len=512) :: file
    integer :: u, ios, i

    if ( .not. obj_ablation%dump_interface_flux ) return
    file = 'OUTPUT/'//trim(FUSS_phase_prefix)//'interface-flux.txt'
    if ( .not. if_file_started .and. obj_sim_param%newrun ) then
      open( newunit=u, file=trim(file), status='replace', action='write', iostat=ios )
      if ( ios == 0 ) write(u,'(A)') '# time b i j k f bs is js ks fs q_if qm_if qh_if   (W, kg/s, W through the face)'
    else
      open( newunit=u, file=trim(file), status='unknown', position='append', action='write', iostat=ios )
    endif
    if_file_started = .true.
    if ( ios /= 0 ) then
      write(*,'(A)') ' [ERROR] cannot open '//trim(file)
      error stop 1
    endif
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 101 ) cycle
      associate ( bc => domain%bc(i) )
        write(u,'(ES24.16,10I8,3ES24.16)') domain%time, bc%b, bc%i, bc%j, bc%k, bc%f, bc%bs, bc%is, bc%js, bc%ks, &
          bc%fs, bc%q_if, bc%qm_if, bc%qh_if
      end associate
    enddo
    close( u )

  end subroutine Write_Interface_Flux


  ! ---------------------------------------------------------------------------

  !> The 4 nodes (i, j, k) of face f of cell (ci, cj, ck).
  pure subroutine face_nodes ( f, ci, cj, ck, nd )
    integer, intent(in)  :: f, ci, cj, ck
    integer, intent(out) :: nd(3,4)
    integer :: a, b, q
    q = 0
    do b = 0, 1
      do a = 0, 1
        q = q + 1
        select case ( f )
          case (1:2) ; nd(:,q) = [ ci - mod(f,2), cj - 1 + a, ck - 1 + b ]
          case (3:4) ; nd(:,q) = [ ci - 1 + a, cj - mod(f,2), ck - 1 + b ]
          case default ; nd(:,q) = [ ci - 1 + a, cj - 1 + b, ck - mod(f,2) ]
        end select
      enddo
    enddo
  end subroutine face_nodes


  pure integer function n_boundary_nodes ( d ) result ( n )
    integer, intent(in) :: d(3)
    n = ( d(1) + 1 ) * ( d(2) + 1 ) * ( d(3) + 1 ) - max( d(1) - 1, 0 ) * max( d(2) - 1, 0 ) * max( d(3) - 1, 0 )
  end function n_boundary_nodes


  !> Length of the receding line through the surface node inc = (b, i, j, k):
  !> the sum of its segments from the pinned face to that node.
  function line_length ( domain, inc ) result ( L )
    use FUSS_Advanced_Types_m
    type(FUSS_domain_type), intent(in) :: domain
    integer,                intent(in) :: inc(4)
    real(R8) :: L
    integer  :: b, d, m, a(3), c(3)
    b = inc(1) ; d = recede_dir(b)
    L = 0.0d0
    a = inc(2:4)
    do m = 1, domain%blk(b)%dim(d)
      a(d) = m - 1 ; c = a ; c(d) = m
      L = L + norm2( domain%blk(b)%node(c(1),c(2),c(3))%c - domain%blk(b)%node(a(1),a(2),a(3))%c )
    enddo
  end function line_length


  !> The shortest cell along that line.
  function line_min_cell ( domain, inc ) result ( h )
    use FUSS_Advanced_Types_m
    type(FUSS_domain_type), intent(in) :: domain
    integer,                intent(in) :: inc(4)
    real(R8) :: h
    integer  :: b, d, m, a(3), c(3)
    b = inc(1) ; d = recede_dir(b)
    h = huge(1.0d0)
    a = inc(2:4)
    do m = 1, domain%blk(b)%dim(d)
      a(d) = m - 1 ; c = a ; c(d) = m
      h = min( h, norm2( domain%blk(b)%node(c(1),c(2),c(3))%c - domain%blk(b)%node(a(1),a(2),a(3))%c ) )
    enddo
  end function line_min_cell


  subroutine push_face ( sn, i )
    type(surf_node_t), intent(inout) :: sn
    integer,           intent(in)    :: i
    integer, allocatable :: tmp(:)
    if ( .not. allocated(sn%face) ) allocate( sn%face(4) )
    if ( sn%nface == size(sn%face) ) then
      allocate( tmp(2*size(sn%face)) ) ; tmp(1:sn%nface) = sn%face(1:sn%nface) ; call move_alloc( tmp, sn%face )
    endif
    sn%nface = sn%nface + 1
    sn%face(sn%nface) = i
  end subroutine push_face


  subroutine push_inc ( sn, id )
    type(surf_node_t), intent(inout) :: sn
    integer,           intent(in)    :: id(4)
    integer, allocatable :: tmp(:,:)
    if ( .not. allocated(sn%inc) ) allocate( sn%inc(4,4) )
    if ( sn%ninc == size(sn%inc,2) ) then
      allocate( tmp(4,2*size(sn%inc,2)) ) ; tmp(:,1:sn%ninc) = sn%inc(:,1:sn%ninc) ; call move_alloc( tmp, sn%inc )
    endif
    sn%ninc = sn%ninc + 1
    sn%inc(:,sn%ninc) = id
  end subroutine push_inc


  pure logical function lex_lt ( a, b )
    real(R8), intent(in) :: a(3), b(3)
    if ( a(1) /= b(1) ) then
      lex_lt = a(1) < b(1)
    elseif ( a(2) /= b(2) ) then
      lex_lt = a(2) < b(2)
    else
      lex_lt = a(3) < b(3)
    endif
  end function lex_lt


  !> perm = the columns of x in lexicographic order (a stable merge sort).
  subroutine sort_coords ( x, perm )
    real(R8), intent(in)  :: x(:,:)
    integer,  intent(out) :: perm(:)
    integer, allocatable :: tmp(:)
    integer :: n, w, lo, mid, hi, i, a, b
    n = size(x,2)
    perm = [( i, i = 1, n )]
    allocate( tmp(n) )
    w = 1
    do while ( w < n )
      lo = 1
      do while ( lo <= n )
        mid = min( lo + w - 1, n ) ; hi = min( lo + 2*w - 1, n )
        a = lo ; b = mid + 1 ; i = lo
        do while ( a <= mid .and. b <= hi )
          if ( lex_lt( x(:,perm(b)), x(:,perm(a)) ) ) then
            tmp(i) = perm(b) ; b = b + 1
          else
            tmp(i) = perm(a) ; a = a + 1
          endif
          i = i + 1
        enddo
        do while ( a <= mid ) ; tmp(i) = perm(a) ; a = a + 1 ; i = i + 1 ; enddo
        do while ( b <= hi )  ; tmp(i) = perm(b) ; b = b + 1 ; i = i + 1 ; enddo
        perm(lo:hi) = tmp(lo:hi)
        lo = lo + 2*w
      enddo
      w = 2*w
    enddo
  end subroutine sort_coords


  !> The first position r of the sorted perm with x(:,perm(r)) not below xq.
  integer function first_not_below ( x, perm, xq ) result ( lo )
    real(R8), intent(in) :: x(:,:), xq(3)
    integer,  intent(in) :: perm(:)
    integer :: hi, mid
    lo = 1 ; hi = size(perm) + 1
    do while ( lo < hi )
      mid = ( lo + hi ) / 2
      if ( lex_lt( x(:,perm(mid)), xq ) ) then
        lo = mid + 1
      else
        hi = mid
      endif
    enddo
  end function first_not_below

end module FUSS_Mod_Recession
