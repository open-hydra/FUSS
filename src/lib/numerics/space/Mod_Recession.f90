!> The ablating wall's per-step and per-output bookkeeping (plan 08 sections
!> 2.1, 2.3 R3, 3.2): the env values of every 306 face at the step's time, the
!> restart state of the faces, and the two files the faces write.
!>
!> ENV VALUES. The env table of a 306 face (t p_w rhoeUeCH0 h_r q_rad chem
!> hconv T_edge) -- for a coupled ablating wall 307, the exchange values
!> env_now (p_w, rhoeUeCH, h_r, q_rad; chem = 1, hconv = 0), which hydra
!> writes before every FUSS step -- is evaluated ONCE per step, after the clock has advanced to
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
!>
!> SMOOTHING (plan 08 section 2.6, thesis Eq. 5.1). Before the displacement,
!> the nodal rate is averaged along each surface grid line over
!> smoothing-width nodes (5: nu_k = (nu_k-2 + .. + nu_k+2)/5), in two passes,
!> one per surface direction; at a line end the window keeps the nodes that
!> exist (truncated). The lines are walked on the surface-node map, so a line
!> continues through a 101 seam whatever the index orientation of the next
!> block (from a node to the neighbour that shares no 306 face with the
!> previous one), and a shared node gets one value. Written as nu_k + mean(
!> nu_j - nu_k), so a uniform rate is returned bit for bit.
!> The surface-node map is built once, at setup, by exact coordinate equality
!> of the block-boundary nodes (they come from one MESH/ file).
module FUSS_Mod_Recession
  use FUSS_Mod_MPI, only: fuss_abort
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Series_Data_m, only: time_table_t

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
  !> smoothing windows: win(m, f, n) = the m-th node of the window of surface
  !> node n along its surface direction f (m = 1..h one side, h+1..2h the
  !> other, nearest first; 0 = beyond the line end), h = (smoothing-width-1)/2
  integer, allocatable, save :: win(:,:,:)
  integer, save :: hwin = 0
  !> TEST ONLY [FUSS-Ablation] coupled-env-file: the exchange values of the
  !> 307 faces of a standalone run (t p_w rhoeUeCH h_r q_rad), one table for
  !> every face (coupled_env), or one per face when the file holds `// face b i
  !> j k f` blocks (env_face, by bc index: hydra's env-dump replayed, plan 12);
  !> neither = hydra
  type(time_table_t), save :: coupled_env
  type(time_table_t), allocatable, save :: env_face(:)

contains

  !> Env values of every 306 and 307 face at domain%time (serial; outside the
  !> parallel region). A 307 face takes its exchange values env_now (from the
  !> coupled-env-file table when that TEST ONLY key is set): hydra writes them
  !> before every step, not before FUSS's setup, so at setup (setup = .true.)
  !> a face without them keeps p-wall; at a step it stops the run. A face
  !> block that opens with p_w = 0 (hydra's values at setup) has none before
  !> its second row.
  subroutine Ablation_Update_Env ( domain, setup )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_ablation
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    logical, intent(in), optional :: setup
    ! Local
    integer  :: i
    real(R8) :: v(7)
    logical  :: at_setup, early

    if ( .not. obj_ablation%enabled ) return
    at_setup = .false.
    if ( present(setup) ) at_setup = setup
    do i = 1, domain%nbound
      select case ( domain%bc(i)%type )
        case (306)
          v = domain%bc(i)%env%eval( domain%time )
          domain%bc(i)%pw        = v(1)
          domain%bc(i)%rhoeUeCH0 = v(2)
          domain%bc(i)%hr        = v(3)
          domain%bc(i)%qrad_env  = v(4)
          domain%bc(i)%chem      = v(5)
          domain%bc(i)%hconv_env = v(6)
          domain%bc(i)%Tedge     = v(7)
        case (307)
          early = .false.
          if ( coupled_env%nrow > 0 ) domain%bc(i)%env_now = coupled_env%eval( domain%time )
          if ( allocated(env_face) ) then
            domain%bc(i)%env_now = env_face(i)%eval( domain%time )
            early = env_face(i)%var(1,1) == 0.0d0 .and. domain%time < env_face(i)%time(2)
          endif
          if ( early .or. .not. ( domain%bc(i)%env_now(1) > 0.0d0 ) ) then
            if ( at_setup ) cycle
            write(*,'(A,I0,A)') ' [ERROR] coupled ablating wall (307), bc.txt record ', i, &
              ': no exchange values (p_w = 0) at this step.'
            if ( early ) then
              write(*,'(A,ES24.16,A,ES24.16,A)') '         t = ', domain%time, &
                ' s precedes the first exchange row of the face''s coupled-env-file block (t = ', env_face(i)%time(2), ' s).'
            else
              write(*,'(A)') '         hydra writes them before every FUSS step; a standalone run needs '// &
                '[FUSS-Ablation] coupled-env-file.'
            endif
            call fuss_abort ( 1 )
          endif
          domain%bc(i)%pw        = domain%bc(i)%env_now(1)
          domain%bc(i)%rhoeUeCH0 = domain%bc(i)%env_now(2)
          domain%bc(i)%hr        = domain%bc(i)%env_now(3)
          domain%bc(i)%qrad_env  = domain%bc(i)%env_now(4)
          domain%bc(i)%chem      = 1.0d0
          domain%bc(i)%hconv_env = 0.0d0
          domain%bc(i)%Tedge     = domain%bc(i)%Tinf
      end select
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
    logical  :: blocks
    logical, allocatable :: seen(:)

    ! TEST ONLY coupled-env-file: the exchange values of the 307 faces
    if ( trim(obj_ablation%coupled_env_file) /= 'none' ) then
      if ( obj_ablation%n_coupled == 0 ) then
        write(*,'(A)') ' [ERROR] [FUSS-Ablation] coupled-env-file is set, but bc.txt has no coupled ablating wall (307).'
        call fuss_abort ( 1 )
      endif
      call Load_Face_Blocks ( domain, 'INPUT/'//trim(obj_ablation%coupled_env_file), blocks )
      if ( .not. blocks ) then
        call coupled_env%load( 'INPUT/'//trim(obj_ablation%coupled_env_file), 4 )
        do i = 1, coupled_env%nrow
          v(1:4) = coupled_env%var(:,i)
          if ( .not. ( v(1) > 0.0d0 ) .or. v(2) < 0.0d0 .or. v(4) < 0.0d0 ) then
            write(*,'(A,I0,A)') ' [ERROR] coupled-env-file INPUT/'//trim(obj_ablation%coupled_env_file)//' row ', i, &
              ': needs p_w > 0, rhoeUeCH >= 0 and q_rad >= 0'
            call fuss_abort ( 1 )
          endif
        enddo
      endif
    endif

    if ( .not. obj_ablation%enabled ) return
    if ( newrun ) return

    file = 'OUTPUT/'//trim(FUSS_phase_prefix)//'bc306-state.dat'
    open( newunit=u, file=trim(file), status='old', action='read', iostat=ios )
    if ( ios /= 0 ) then
      write(*,'(A)') ' [ERROR] restart with ablating faces (306): '//trim(file)//' is missing'
      write(*,'(A)') '         it holds each face''s recession, wall temperature, recession rate and Picard start'
      call fuss_abort ( 1 )
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
        if ( domain%bc(i)%type /= 306 .and. domain%bc(i)%type /= 307 ) cycle
        if ( domain%bc(i)%b == b .and. domain%bc(i)%i == ii .and. domain%bc(i)%j == jj .and. &
             domain%bc(i)%k == kk .and. domain%bc(i)%f == f ) exit
      enddo
      if ( i > domain%nbound ) call Fail( 'no 306 or 307 face of bc.txt matches this row' )
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
      if ( ( domain%bc(i)%type == 306 .or. domain%bc(i)%type == 307 ) .and. .not. seen(i) ) then
        write(*,'(A,I0,A)') ' [ERROR] '//trim(file)//': no row for the 306 face of bc.txt record ', i, &
          ' (written by a run with another bc.txt?)'
        call fuss_abort ( 1 )
      endif
    enddo

  contains

    subroutine Fail ( reason )
      character(len=*), intent(in) :: reason
      write(*,'(A,I0,A)') ' [ERROR] '//trim(file)//' line ', ln, ': '//reason
      call fuss_abort ( 1 )
    end subroutine Fail

  end subroutine Setup_Ablation


  !> The `// face b i j k f` blocks of a coupled-env-file (plan 12: hydra's
  !> env-dump replayed in a standalone run). Such a line opens the block of
  !> one 307 face (its bc indices, as in bc306-state.dat and hydra's
  !> wall-file); the rows below it, `t p_w rhoeUeCH h_r q_rad`, are that
  !> face's table (linear in t, clamped, as the one-table file). Every 307 face
  !> has exactly one block. A block's first row may carry p_w = 0: no exchange
  !> values yet, hydra's state at FUSS's setup (Ablation_Update_Env). Other
  !> `//` lines are comments. A file without a `// face` line is the one-table
  !> file: found = .false., and nothing is checked here.
  subroutine Load_Face_Blocks ( domain, file, found )
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(in)  :: domain
    character(len=*),       intent(in)  :: file
    logical,                intent(out) :: found
    ! Local
    character(len=1024) :: line, t, s
    integer  :: u, ios, ln, nb, nr, k, i, ib, r, orphan
    integer,  allocatable :: key(:,:), hline(:), owner(:), rline(:)
    real(R8), allocatable :: rows(:,:)
    real(R8) :: v(5)

    found = .false.
    open( newunit=u, file=trim(file), status='old', action='read', iostat=ios )
    if ( ios /= 0 ) return                   ! the one-table loader names the file
    allocate( key(5,16), hline(16), rows(5,64), owner(64), rline(64) )
    nb = 0 ; nr = 0 ; ln = 0 ; orphan = 0
    do
      read( u, '(A)', iostat=ios ) line
      if ( ios /= 0 ) exit
      ln = ln + 1
      t = adjustl(line)
      if ( len_trim(t) == 0 ) cycle
      if ( t(1:2) == '//' ) then
        s = adjustl(t(3:))
        if ( s(1:4) /= 'face' .or. s(5:5) /= ' ' ) cycle      ! a comment
        if ( nb == size(key,2) ) call Grow_Blocks ()
        nb = nb + 1
        hline(nb) = ln
        read( s(5:), *, iostat=ios ) key(:,nb)
        if ( ios /= 0 ) call Fail( 'a block header is "// face b i j k f": five integers, the bc indices of a 307 face' )
        cycle
      endif
      if ( nb == 0 ) then
        if ( orphan == 0 ) orphan = ln
        cycle
      endif
      read( t, *, iostat=ios ) v
      if ( ios /= 0 ) call Fail( 'expected 5 numbers (t p_w rhoeUeCH h_r q_rad)' )
      if ( any( v /= v ) .or. any( abs(v) > huge(1.0d0) ) ) call Fail( 'value is not finite' )
      if ( nr > 0 ) then
        if ( owner(nr) == nb .and. v(1) <= rows(1,nr) ) call Fail( 't must increase strictly from row to row' )
      endif
      if ( nr == size(rows,2) ) call Grow_Rows ()
      nr = nr + 1
      rows(:,nr) = v
      owner(nr) = nb
      rline(nr) = ln
    enddo
    close( u )
    if ( nb == 0 ) return
    found = .true.
    if ( orphan > 0 ) then
      ln = orphan
      call Fail( 'a row before the first "// face" line: with face blocks, every row belongs to the block above it' )
    endif

    allocate( env_face(domain%nbound) )
    do k = 1, nb
      ln = hline(k)
      write(s,'(A,5(1X,I0),A)') '"// face', key(:,k), '"'
      if ( count( owner(1:nr) == k ) < 2 ) call Fail( 'the block '//trim(s)//' has fewer than two rows' )
      ib = 0
      do i = 1, domain%nbound
        if ( domain%bc(i)%type /= 307 ) cycle
        if ( domain%bc(i)%b == key(1,k) .and. domain%bc(i)%i == key(2,k) .and. domain%bc(i)%j == key(3,k) .and. &
             domain%bc(i)%k == key(4,k) .and. domain%bc(i)%f == key(5,k) ) then
          ib = i
          exit
        endif
      enddo
      if ( ib == 0 ) call Fail( 'the block '//trim(s)//' names no coupled ablating wall (307) of bc.txt' )
      if ( env_face(ib)%nrow > 0 ) call Fail( 'a second block for the 307 face '//trim(s) )
      env_face(ib)%n    = 4
      env_face(ib)%nrow = count( owner(1:nr) == k )
      allocate( env_face(ib)%time(env_face(ib)%nrow), env_face(ib)%var(4, env_face(ib)%nrow) )
      r = 0
      do i = 1, nr
        if ( owner(i) /= k ) cycle
        r = r + 1
        env_face(ib)%time(r)  = rows(1,i)
        env_face(ib)%var(:,r) = rows(2:5,i)
        if ( ( .not. ( rows(2,i) > 0.0d0 ) .and. .not. ( r == 1 .and. rows(2,i) == 0.0d0 ) ) .or. &
             rows(3,i) < 0.0d0 .or. rows(5,i) < 0.0d0 ) then
          ln = rline(i)
          call Fail( 'needs p_w > 0 (p_w = 0 only on the first row of a block: no exchange values yet), '// &
                     'rhoeUeCH >= 0 and q_rad >= 0' )
        endif
      enddo
    enddo
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 307 .or. env_face(i)%nrow > 0 ) cycle
      write(*,'(A,I0,A,5(1X,I0),A)') ' [ERROR] coupled-env-file '//trim(file)//': no "// face" block for the 307 face '// &
        'of bc.txt record ', i, ' ("// face', domain%bc(i)%b, domain%bc(i)%i, domain%bc(i)%j, domain%bc(i)%k, &
        domain%bc(i)%f, '"); with face blocks every 307 face needs exactly one'
      call fuss_abort ( 1 )
    enddo

  contains

    subroutine Fail ( reason )
      character(len=*), intent(in) :: reason
      write(*,'(A,I0,A)') ' [ERROR] coupled-env-file '//trim(file)//' line ', ln, ': '//reason
      call fuss_abort ( 1 )
    end subroutine Fail

    subroutine Grow_Blocks ()
      integer, allocatable :: ki(:,:), hi(:)
      allocate( ki(5, 2*size(key,2)), hi(2*size(key,2)) )
      ki(:,1:size(key,2)) = key
      hi(1:size(key,2))   = hline
      call move_alloc( ki, key )
      call move_alloc( hi, hline )
    end subroutine Grow_Blocks

    subroutine Grow_Rows ()
      real(R8), allocatable :: tr(:,:)
      integer,  allocatable :: to(:), tl(:)
      allocate( tr(5, 2*size(rows,2)), to(2*size(rows,2)), tl(2*size(rows,2)) )
      tr(:,1:size(rows,2)) = rows
      to(1:size(rows,2))   = owner
      tl(1:size(rows,2))   = rline
      call move_alloc( tr, rows )
      call move_alloc( to, owner )
      call move_alloc( tl, rline )
    end subroutine Grow_Rows

  end subroutine Load_Face_Blocks


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
      call fuss_abort ( 1 )
    endif
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 306 .and. domain%bc(i)%type /= 307 ) cycle
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
      call fuss_abort ( 1 )
    endif
    write(u,'(A)') '# b i j k f s_total Tw_prev sdot beta Bpg_lag'
    write(u,'(A,ES24.16)') '# min-thickness ', min_thick
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 306 .and. domain%bc(i)%type /= 307 ) cycle
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
    integer, allocatable :: bid(:,:), perm(:), sfi(:), sperm(:), pnode(:)
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
    allocate( smap(cnt), pnode(nsf) )
    cnt = 0
    do n = 1, nsf
      if ( n == 1 ) then
        cnt = 1
      elseif ( any( sx(:,sperm(n)) /= sx(:,sperm(n-1)) ) ) then
        cnt = cnt + 1
      endif
      smap(cnt)%x = sx(:,sperm(n))
      call push_face( smap(cnt), sfi(sperm(n)) )
      pnode(sperm(n)) = cnt
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
    if ( bad ) call fuss_abort ( 1 )

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

    ! the smoothing windows (pair 4 (o - 1) + q is corner q of the o-th 306 face)
    call Setup_Smoothing ( reshape( pnode(1:nsf), [4, nsf/4] ) )

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
    real(R8), allocatable :: sdn_n(:)

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

    ! the nodal rates: the area-weighted mean of the adjacent 306 faces, then
    ! the smoothing along the surface lines
    allocate( sdn_n(size(smap)) )
    do n = 1, size(smap)
      sum_a = 0.0d0 ; sum_as = 0.0d0
      do q = 1, smap(n)%nface
        associate ( bc => domain%bc(smap(n)%face(q)) )
          call Face_Index ( bc%f, dr, bc%i, bc%j, bc%k, fi, fj, fk )
          A = domain%blk(bc%b)%dir(dr)%f(fi,fj,fk)%A
          sum_a  = sum_a  + A
          sum_as = sum_as + A * bc%sdot
        end associate
      enddo
      sdn_n(n) = sum_as / sum_a
    enddo
    if ( hwin > 0 ) then
      call Smooth_Pass ( sdn_n, 1 )
      call Smooth_Pass ( sdn_n, 2 )
    endif

    do n = 1, size(smap)
      sdn = sdn_n(n)
      if ( sdn == 0.0d0 ) cycle               ! no motion: the nodes stay bitwise
      nrm = 0.0d0
      do q = 1, smap(n)%nface
        associate ( bc => domain%bc(smap(n)%face(q)) )
          call Face_Index ( bc%f, dr, bc%i, bc%j, bc%k, fi, fj, fk )
          A = domain%blk(bc%b)%dir(dr)%f(fi,fj,fk)%A
          nrm = nrm + A * domain%blk(bc%b)%dir(dr)%f(fi,fj,fk)%N * real( 1 - 2 * mod(bc%f,2), R8 )
        end associate
      enddo
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
          call fuss_abort ( 2 )
        endif
      enddo
    enddo

  end subroutine Compute_Node_Recession


  !> The smoothing windows of every surface node (setup, law = nodeshift):
  !> fc(q, o) = the surface node of corner q of the o-th 306 face, corners in
  !> face_nodes' order (0,0) (1,0) (0,1) (1,1) of the face's two transverse
  !> directions, so the face's edges along direction 1 are corners 1-2 and
  !> 3-4, along direction 2 corners 1-3 and 2-4. A node's window along its
  !> direction f starts from its edge neighbour along f in its first face and
  !> walks both ways: the next node after (prev, cur) is the neighbour of cur
  !> that shares no 306 face with prev (on a quad surface, the node straight
  !> ahead), until the line ends or h nodes are taken.
  subroutine Setup_Smoothing ( fc )
    use FUSS_Config_Types_m, only: obj_ablation
    implicit none
    integer, intent(in) :: fc(:,:)
    ! Local
    integer :: nn, o, n, f, m, q, a, prev, cur, nxt, side, first
    integer, allocatable :: nbr(:,:), nnb(:), firstf(:)
    integer, parameter :: MAXNB = 16
    integer, parameter :: EDGES(2,4) = reshape( [1,2, 3,4, 1,3, 2,4], [2,4] )   ! dir 1: 1-2, 3-4; dir 2: 1-3, 2-4

    hwin = ( obj_ablation%smoothing_width - 1 ) / 2
    if ( allocated(win) ) deallocate( win )
    if ( hwin <= 0 ) return
    nn = size(smap)
    allocate( nbr(MAXNB, nn), nnb(nn), win(2*hwin, 2, nn), firstf(nn) )
    nbr = 0 ; nnb = 0 ; win = 0 ; firstf = 0
    do o = 1, size(fc, 2)
      do q = 1, 4
        call add_nbr( fc(EDGES(1,q),o), fc(EDGES(2,q),o) )
        call add_nbr( fc(EDGES(2,q),o), fc(EDGES(1,q),o) )
        if ( firstf(fc(q,o)) == 0 ) firstf(fc(q,o)) = o
      enddo
    enddo
    do n = 1, nn
      first = firstf(n)
      do f = 1, 2
        ! the edge neighbour of n along direction f in the first face containing n
        a = 0
        do q = 1, 2
          if ( fc(EDGES(1,2*(f-1)+q),first) == n ) a = fc(EDGES(2,2*(f-1)+q),first)
          if ( fc(EDGES(2,2*(f-1)+q),first) == n ) a = fc(EDGES(1,2*(f-1)+q),first)
        enddo
        do side = 1, 2
          prev = n
          cur  = a
          if ( side == 2 ) cur = straight( a, n )          ! the other side: ahead of n, coming from a
          if ( side == 2 ) prev = n
          do m = 1, hwin
            if ( cur == 0 ) exit
            win( (side-1)*hwin + m, f, n ) = cur
            nxt  = straight( prev, cur )
            prev = cur
            cur  = nxt
          enddo
        enddo
      enddo
    enddo

  contains

    subroutine add_nbr ( u, v )
      integer, intent(in) :: u, v
      if ( any( nbr(1:nnb(u), u) == v ) ) return
      if ( nnb(u) == MAXNB ) then
        write(*,'(A,I0,A)') ' [ERROR] smoothing: surface node ', u, ' has too many neighbours (not a quad surface)'
        call fuss_abort ( 1 )
      endif
      nnb(u) = nnb(u) + 1
      nbr(nnb(u), u) = v
    end subroutine add_nbr

    !> the neighbour of cur that shares no 306 face with prev (0 = none)
    integer function straight ( prev, cur ) result ( nx )
      integer, intent(in) :: prev, cur
      integer :: k
      nx = 0
      if ( prev == 0 .or. cur == 0 ) return
      do k = 1, nnb(cur)
        if ( nbr(k,cur) == prev ) cycle
        if ( share_face( nbr(k,cur), prev ) ) cycle
        nx = nbr(k,cur)
        return
      enddo
    end function straight

    logical function share_face ( u, v )
      integer, intent(in) :: u, v
      integer :: k
      share_face = .false.
      do k = 1, smap(u)%nface
        if ( any( smap(v)%face(1:smap(v)%nface) == smap(u)%face(k) ) ) then
          share_face = .true.
          return
        endif
      enddo
    end function share_face

  end subroutine Setup_Smoothing


  !> One smoothing pass along surface direction f: nu_n + mean over the
  !> window (n itself and its existing window nodes) of nu_j - nu_n.
  subroutine Smooth_Pass ( v, f )
    real(R8), intent(inout) :: v(:)
    integer,  intent(in)    :: f
    real(R8) :: w(size(v)), dsum
    integer  :: n, m, cnt
    do n = 1, size(v)
      dsum = 0.0d0 ; cnt = 1
      do m = 1, 2*hwin
        if ( win(m,f,n) == 0 ) cycle
        dsum = dsum + ( v(win(m,f,n)) - v(n) )
        cnt  = cnt + 1
      enddo
      w(n) = v(n) + dsum / real( cnt, R8 )
    enddo
    v = w
  end subroutine Smooth_Pass


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
            call fuss_abort ( 2 )
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
      call fuss_abort ( 1 )
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
