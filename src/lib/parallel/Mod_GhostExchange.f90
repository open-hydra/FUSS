!>@brief Ghost cell MPI communication for inter-block boundaries.
!> Builds a communication schedule from the BC connectivity array,
!> then provides non-blocking exchange routines for T and R arrays.
!> Messages are aggregated per remote rank to minimize latency.
!> In serial mode (USE_MPI not defined), all routines are no-ops.
module FUSS_Mod_GhostExchange
  use iso_fortran_env, only: R8 => real64, I4 => int32
  use FUSS_Mod_MPI

  implicit none
  private

  !> Single ghost cell exchange entry
  type :: ghost_entry_type
    integer :: bc_idx          !< Index into domain%bc(:)
    integer :: local_block     !< Block on this rank
    integer :: remote_block    !< Block on the remote rank
    integer :: remote_rank     !< MPI rank of the remote block
    integer :: src_i, src_j, src_k, src_f  !< Source cell coords and face
    integer :: dst_i, dst_j, dst_k, dst_f  !< Destination cell coords and face
  end type ghost_entry_type

  !> Per-rank aggregation info
  type :: rank_group_type
    integer :: rank = -1       !< Remote MPI rank
    integer :: offset = 0      !< Start index in entry list (1-based)
    integer :: count  = 0      !< Number of entries for this rank
  end type rank_group_type

  !> Single chimera donor-cell exchange entry (BC type 102).
  !> Identifies one donor cell (bs,is,js,ks) whose T must travel from the
  !> donor's owner rank to the receiver block's owner rank.
  type :: chim_entry_type
    integer :: remote_rank     !< Peer rank (receiver owner on send side, donor owner on recv side)
    integer :: bs              !< Donor block
    integer :: is, js, ks      !< Donor cell indices
  end type chim_entry_type

  !> Face group: contiguous run of entries sharing (block, face).
  !> Precomputes guide offsets; for rectangular faces enables direct coordinate sweep.
  type :: face_group_type
    integer :: block            !< Local block index
    integer :: face             !< Face number (1-6)
    integer :: entry_start      !< First entry index in send/recv list
    integer :: entry_count      !< Number of entries
    integer :: di, dj, dk       !< Precomputed guide(face,:)
  end type face_group_type

  !> The halo of one field-list row (plan 11 section 2.2): gc layers of ncomp values per
  !> entry of the schedule's lists, one persistent send and receive per rank group, tag
  !> = the sender's rank + tag_base. Every exchanged field travels this one path; the row's
  !> code names its array to Field_Cell_Get / Field_Cell_Set.
  type :: halo_field_type
    integer :: code = 0, ncomp = 0, tag_base = 0
    logical :: init = .false.
    real(R8), allocatable :: send_buf(:)     !< (ncomp * gc * n_send)
    real(R8), allocatable :: recv_buf(:)     !< (ncomp * gc * n_recv)
    integer,  allocatable :: send_req(:)     !< (n_send_ranks), persistent
    integer,  allocatable :: recv_req(:)     !< (n_recv_ranks), persistent
  end type halo_field_type

  !> Communication schedule: entries grouped by remote rank for aggregated messaging
  type :: ghost_schedule_type
    integer :: n_send = 0, n_recv = 0
    type(ghost_entry_type), allocatable :: send_list(:)
    type(ghost_entry_type), allocatable :: recv_list(:)
    logical :: built = .false.

    ! Per-rank aggregation for sends
    integer :: n_send_ranks = 0
    type(rank_group_type), allocatable :: send_groups(:)

    ! Per-rank aggregation for receives
    integer :: n_recv_ranks = 0
    type(rank_group_type), allocatable :: recv_groups(:)

    ! Face groups for aggregated pack/unpack (reduce per-cell overhead)
    integer :: n_send_faces = 0
    type(face_group_type), allocatable :: send_faces(:)
    integer :: n_recv_faces = 0
    type(face_group_type), allocatable :: recv_faces(:)

    ! Pre-allocated MPI buffers (one contiguous buffer per rank group); the T halo
    ! itself travels in halo(:), these carry R (IRS) and the one-shot matID
    real(R8), allocatable :: T_send_buf(:)   !< (gc * n_send)
    real(R8), allocatable :: T_recv_buf(:)   !< (gc * n_recv)
    real(R8), allocatable :: Tg_send_buf(:)  !< (4  * n_send)
    real(R8), allocatable :: Tg_recv_buf(:)  !< (4  * n_recv)
    integer, allocatable  :: send_req(:)     !< (n_send_ranks)
    integer, allocatable  :: recv_req(:)     !< (n_recv_ranks)
    integer, allocatable  :: mpi_stat(:,:)   !< (MPI_STATUS_SIZE, max(n_send_ranks, n_recv_ranks))
    ! Persistent MPI requests for R exchange (initialized once, started each step)
    integer, allocatable  :: R_send_req_pers(:)  !< (n_send_ranks)
    integer, allocatable  :: R_recv_req_pers(:)  !< (n_recv_ranks)
    logical :: persistent_R_init = .false.

    ! The halo of every field-list row with exchange and has_ghost (T, chi, p_pg), by row
    type(halo_field_type), allocatable :: halo(:)

    ! Chimera (type 102) donor-cell exchange: entries, per-rank groups, buffers
    integer :: n_chim_send = 0, n_chim_recv = 0
    type(chim_entry_type), allocatable :: chim_send_list(:)
    type(chim_entry_type), allocatable :: chim_recv_list(:)
    integer :: n_chim_send_ranks = 0
    integer :: n_chim_recv_ranks = 0
    type(rank_group_type), allocatable :: chim_send_groups(:)
    type(rank_group_type), allocatable :: chim_recv_groups(:)
    real(R8), allocatable :: chim_send_buf(:)   !< (n_chim_send)
    real(R8), allocatable :: chim_recv_buf(:)   !< (n_chim_recv)
    integer, allocatable  :: chim_send_req(:)   !< (n_chim_send_ranks)
    integer, allocatable  :: chim_recv_req(:)   !< (n_chim_recv_ranks)
    integer, allocatable  :: chim_stat(:,:)     !< (MPI_STATUS_SIZE, max ranks)
  end type ghost_schedule_type

  !> One schedule per multigrid level: each level has its own BC file (and so
  !> its own connections and chimera donors). ghost_sched points at the level
  !> being worked on; every entry point that exchanges data selects it first.
  type(ghost_schedule_type), allocatable, target :: ghost_scheds(:)
  type(ghost_schedule_type), pointer, public :: ghost_sched => null()

  public :: allocate_exchange_schedules, set_active_mg_level

  ! Node planes of the block connections between two ranks (plan 11 stage 2), built at the first
  ! exchange: transfer t sends the planes 0..gc of face np_face(t) of block np_block(t) from that
  ! block's rank to rank np_to(t)
  logical :: np_built = .false.
  integer :: np_n = 0
  integer, allocatable :: np_block(:), np_face(:), np_to(:)
  integer, parameter :: NODE_GATHER_TAG = 70000, NODE_PLANE_TAG = 80000

  ! The halo of field code c is tagged HALO_TAG_STEP * (c - 1) + the sender's rank (T keeps
  ! the rank alone, as before plan 11 stage 3); the field-list row of T, set with the halos
  integer, parameter :: HALO_TAG_STEP = 100000
  integer :: fl_T = 0
  ! gather_field_to_root: block b's message is tagged FIELD_GATHER_TAG + b (Intel MPI 2021.9:
  ! MPI_TAG_UB = 524287)
  integer, parameter :: FIELD_GATHER_TAG = 400000
  ! the gas-pressure stencil pg(3:6) (Tg(3:6) keeps 10000)
  integer, parameter :: PG_TAG = 300000

  public :: build_ghost_schedule
  public :: build_local_bc_index
  public :: cleanup_ghost_schedule
  public :: exchange_ghost_T_post_recv, exchange_ghost_T_pack
  public :: exchange_ghost_T_post_send, exchange_ghost_T_wait_unpack
  public :: exchange_ghost_T_wait_send
  public :: exchange_ghost_post_recv, exchange_ghost_pack, exchange_ghost_post_send
  public :: exchange_ghost_wait_unpack, exchange_ghost_wait_send
  public :: exchange_ghost_R
  public :: exchange_ghost_Tg, exchange_ghost_pg
  public :: entry_send_count, entry_recv_count, entry_send_cell, entry_recv_cell, exchange_entry_values
  public :: exchange_ghost_chimera_post_recv, exchange_ghost_chimera_pack
  public :: exchange_ghost_chimera_post_send, exchange_ghost_chimera_wait_recv
  public :: exchange_ghost_chimera_unpack, exchange_ghost_chimera_wait_send
  public :: gather_T_to_root, gather_matID_to_root, gather_qvol_to_root
  public :: gather_node_to_root, exchange_node_planes
  public :: gather_field_to_root
  public :: gather_diagnostic_to_root
  public :: scatter_T_from_root
  public :: mpi_io_barrier
  public :: Ghost_Interrank

contains


  subroutine allocate_exchange_schedules(nlev)
    integer, intent(in) :: nlev

    if (allocated(ghost_scheds)) deallocate(ghost_scheds)
    allocate(ghost_scheds(max(nlev, 1)))
    ghost_sched => ghost_scheds(1)
  end subroutine allocate_exchange_schedules


  subroutine set_active_mg_level(lev)
    integer, intent(in) :: lev

    if (.not. allocated(ghost_scheds)) call allocate_exchange_schedules(lev)
    if (lev < 1 .or. lev > size(ghost_scheds)) then
      error stop 'set_active_mg_level: level outside allocated exchange schedules'
    end if
    ghost_sched => ghost_scheds(lev)
  end subroutine set_active_mg_level


  !> Build the communication schedule by scanning all BC entries of type 1 (connection).
  !> Entries are sorted by remote rank for aggregated MPI messaging.
  !> Must be called after partition_blocks and domain setup.
  subroutine build_ghost_schedule(domain)
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, ns, nr, bm, bs

    call set_active_mg_level(domain%mg_level)

    if (mpi_size_ <= 1) then
      ghost_sched%built = .true.
      return
    end if

#ifdef USE_MPI
    ! Count send and recv entries
    ns = 0; nr = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 101) cycle
      bm = domain%bc(i)%b   ! destination block
      bs = domain%bc(i)%bs  ! source block
      if (.not. allocated(block_owner)) cycle
      ! If I own the destination but not the source => I need to RECV
      if (is_local_block(bm) .and. (.not. is_local_block(bs))) then
        nr = nr + 1
      end if
      ! If I own the source but not the destination => I need to SEND
      if (is_local_block(bs) .and. (.not. is_local_block(bm))) then
        ns = ns + 1
      end if
    end do

    ghost_sched%n_send = ns
    ghost_sched%n_recv = nr
    if (allocated(ghost_sched%send_list)) deallocate(ghost_sched%send_list)
    if (allocated(ghost_sched%recv_list)) deallocate(ghost_sched%recv_list)
    allocate(ghost_sched%send_list(ns))
    allocate(ghost_sched%recv_list(nr))

    ! Fill send and recv lists
    ns = 0; nr = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 101) cycle
      bm = domain%bc(i)%b
      bs = domain%bc(i)%bs

      if (is_local_block(bm) .and. (.not. is_local_block(bs))) then
        nr = nr + 1
        ghost_sched%recv_list(nr)%bc_idx = i
        ghost_sched%recv_list(nr)%local_block  = bm
        ghost_sched%recv_list(nr)%remote_block = bs
        ghost_sched%recv_list(nr)%remote_rank  = block_owner(bs)
        ghost_sched%recv_list(nr)%src_i = domain%bc(i)%is
        ghost_sched%recv_list(nr)%src_j = domain%bc(i)%js
        ghost_sched%recv_list(nr)%src_k = domain%bc(i)%ks
        ghost_sched%recv_list(nr)%src_f = domain%bc(i)%fs
        ghost_sched%recv_list(nr)%dst_i = domain%bc(i)%i
        ghost_sched%recv_list(nr)%dst_j = domain%bc(i)%j
        ghost_sched%recv_list(nr)%dst_k = domain%bc(i)%k
        ghost_sched%recv_list(nr)%dst_f = domain%bc(i)%f
      end if

      if (is_local_block(bs) .and. (.not. is_local_block(bm))) then
        ns = ns + 1
        ghost_sched%send_list(ns)%bc_idx = i
        ghost_sched%send_list(ns)%local_block  = bs
        ghost_sched%send_list(ns)%remote_block = bm
        ghost_sched%send_list(ns)%remote_rank  = block_owner(bm)
        ghost_sched%send_list(ns)%src_i = domain%bc(i)%is
        ghost_sched%send_list(ns)%src_j = domain%bc(i)%js
        ghost_sched%send_list(ns)%src_k = domain%bc(i)%ks
        ghost_sched%send_list(ns)%src_f = domain%bc(i)%fs
        ghost_sched%send_list(ns)%dst_i = domain%bc(i)%i
        ghost_sched%send_list(ns)%dst_j = domain%bc(i)%j
        ghost_sched%send_list(ns)%dst_k = domain%bc(i)%k
        ghost_sched%send_list(ns)%dst_f = domain%bc(i)%f
      end if
    end do

    ! Sort entries by remote rank and build rank groups
    call sort_entries_by_rank(ghost_sched%send_list, ns)
    call sort_entries_by_rank(ghost_sched%recv_list, nr)
    call build_rank_groups(ghost_sched%send_list, ns, ghost_sched%send_groups, ghost_sched%n_send_ranks)
    call build_rank_groups(ghost_sched%recv_list, nr, ghost_sched%recv_groups, ghost_sched%n_recv_ranks)

    ! Build face groups for aggregated pack/unpack
    call build_face_groups_send(ghost_sched%send_list, ns, &
                                ghost_sched%send_faces, ghost_sched%n_send_faces)
    call build_face_groups_recv(ghost_sched%recv_list, nr, &
                                ghost_sched%recv_faces, ghost_sched%n_recv_faces)

    ! Pre-allocate MPI buffers (reused every RK step)
    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv
    if (ns > 0 .or. nr > 0) then
      block
        use mpi, only: MPI_STATUS_SIZE
        integer :: max_ranks
        max_ranks = max(ghost_sched%n_send_ranks, ghost_sched%n_recv_ranks, 1)
        allocate(ghost_sched%T_send_buf(gc * max(ns,1)))
        allocate(ghost_sched%T_recv_buf(gc * max(nr,1)))
        allocate(ghost_sched%Tg_send_buf(4 * max(ns,1)))
        allocate(ghost_sched%Tg_recv_buf(4 * max(nr,1)))
        allocate(ghost_sched%send_req(max(ghost_sched%n_send_ranks,1)))
        allocate(ghost_sched%recv_req(max(ghost_sched%n_recv_ranks,1)))
        allocate(ghost_sched%mpi_stat(MPI_STATUS_SIZE, max_ranks))
      end block
    end if

    ! Persistent MPI requests: the halo of every exchanged field (T first), then R
    call init_halo_fields()
    call init_persistent_R_requests()

    ! One-shot exchange of matID (static field) so ghost cells have correct material IDs
    call exchange_matID_once(domain)

    ! Build chimera donor-cell exchange schedule (BC type 102)
    call build_chimera_schedule(domain)

    if (mpi_is_root) then
      write(*,'(A,I0,A,I0,A)') ' Ghost schedule: ', ghost_sched%n_send_ranks, &
        ' send ranks, ', ghost_sched%n_recv_ranks, ' recv ranks'
    end if
#endif

    ghost_sched%built = .true.
  end subroutine build_ghost_schedule


  !> Build pre-filtered list of BC indices where Bm is local to this rank.
  !> Avoids O(nbound) loop scanning in hot paths (Fill_Ghost_Cell, etc.).
  subroutine build_local_bc_index(domain)
    use FUSS_Advanced_Types_m

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, n, ns

    ! Count local BC entries (Bm is local) and local source entries (Bs is local, type-1)
    n = 0; ns = 0
    do i = 1, domain%nbound
      if (is_local_block(domain%bc(i)%b)) n = n + 1
      if (domain%bc(i)%type == 101) then
        if (is_local_block(domain%bc(i)%bs)) ns = ns + 1
      end if
    end do

    ! Build local_bc_idx (Bm is local)
    domain%n_local_bc = n
    if (allocated(domain%local_bc_idx)) deallocate(domain%local_bc_idx)
    allocate(domain%local_bc_idx(n))

    ! Build local_bs_idx (Bs is local, type-1 only — for Tg computation)
    domain%n_local_bs = ns
    if (allocated(domain%local_bs_idx)) deallocate(domain%local_bs_idx)
    allocate(domain%local_bs_idx(ns))

    ! Fill both index arrays in a single pass
    n = 0; ns = 0
    do i = 1, domain%nbound
      if (is_local_block(domain%bc(i)%b)) then
        n = n + 1
        domain%local_bc_idx(n) = i
      end if
      if (domain%bc(i)%type == 101) then
        if (is_local_block(domain%bc(i)%bs)) then
          ns = ns + 1
          domain%local_bs_idx(ns) = i
        end if
      end if
    end do
  end subroutine build_local_bc_index


  !> Free persistent MPI requests and deallocate schedule arrays.
  subroutine cleanup_ghost_schedule()
    implicit none
#ifdef USE_MPI
    integer :: l

    if (allocated(ghost_scheds)) then
      do l = 1, size(ghost_scheds)
        call set_active_mg_level(l)
        call cleanup_persistent_requests()
      end do
    end if
#endif
  end subroutine cleanup_ghost_schedule


#ifdef USE_MPI
  !> The halos of the active schedule (plan 11 section 2.2): one per field-list row with
  !> exchange and has_ghost (T always; chi with pyrolysis, p_pg under type 2), each with
  !> its buffers and persistent requests. Called once per level from build_ghost_schedule.
  !> A row without an accessor code stops the setup: it would be exchanged by no one.
  subroutine init_halo_fields()
    use FUSS_Field_List_m, only: Field_List_Count, Field_List_Get, Field_List_Lookup, field_entry_t
    use mpi
    implicit none
    integer :: fl, ierr, top
    integer(kind=MPI_ADDRESS_KIND) :: tag_ub
    logical :: flag
    type(field_entry_t) :: e

    fl_T = Field_List_Lookup('T')
    if (allocated(ghost_sched%halo)) deallocate(ghost_sched%halo)
    allocate(ghost_sched%halo(Field_List_Count()))

    top = 0
    do fl = 1, Field_List_Count()
      e = Field_List_Get(fl)
      if (.not. (e%exchange .and. e%has_ghost)) cycle
      if (e%code == 0) then
        write(*,'(A)') ' [ERROR] init_halo_fields: no accessor code for the exchanged field '//trim(e%name)
        call fuss_abort(1)
      end if
      call init_halo_field(ghost_sched%halo(fl), e%code, e%ncomp)
      top = max(top, ghost_sched%halo(fl)%tag_base + mpi_size_ - 1)
    end do

    ! every tag must fit MPI_TAG_UB (the standard guarantees only 32767; Intel MPI 2021.9
    ! gives 524287), the gathers' too
    top = max(top, FIELD_GATHER_TAG + size(block_owner), BLOCK_GATHER_TAG + size(block_owner))
    call MPI_COMM_GET_ATTR(MPI_COMM_WORLD, MPI_TAG_UB, tag_ub, flag, ierr)
    call check_mpi_error(ierr)
    if (flag .and. tag_ub < top) then
      write(*,'(A,I0,A,I0)') ' [ERROR] MPI_TAG_UB = ', tag_ub, ' is below the largest halo tag ', top
      call fuss_abort(1)
    end if
  end subroutine init_halo_fields


  !> Buffers and persistent requests of one field's halo (MPI_SEND_INIT / MPI_RECV_INIT).
  subroutine init_halo_field(h, code, ncomp)
    use FUSS_Global_m, only: gc
    use mpi

    implicit none
    type(halo_field_type), intent(inout) :: h
    integer, intent(in) :: code, ncomp
    integer :: r, ierr, tag, buf_pos, entry_size

    h%code = code
    h%ncomp = ncomp
    h%tag_base = HALO_TAG_STEP * (code - 1)
    entry_size = ncomp * gc

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return

    allocate(h%send_buf(entry_size * max(ghost_sched%n_send, 1)))
    allocate(h%recv_buf(entry_size * max(ghost_sched%n_recv, 1)))
    allocate(h%recv_req(max(ghost_sched%n_recv_ranks, 1)))
    allocate(h%send_req(max(ghost_sched%n_send_ranks, 1)))

    ! Initialize persistent receives
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank + h%tag_base
      call MPI_RECV_INIT(h%recv_buf(buf_pos), &
                         ghost_sched%recv_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%recv_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, h%recv_req(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Initialize persistent sends
    do r = 1, ghost_sched%n_send_ranks
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_ + h%tag_base
      call MPI_SEND_INIT(h%send_buf(buf_pos), &
                         ghost_sched%send_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%send_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, h%send_req(r), ierr)
      call check_mpi_error(ierr)
    end do

    h%init = .true.
  end subroutine init_halo_field


  !> Initialize persistent MPI requests for R field exchange.
  !> Reuses T_send_buf/T_recv_buf (same size). Uses distinct tags (offset 30000).
  subroutine init_persistent_R_requests()
    use FUSS_Global_m, only: gc
    use mpi

    implicit none
    integer :: r, ierr, tag, buf_pos, entry_size
    integer, parameter :: TAG_OFFSET = 30000

    entry_size = gc

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return

    allocate(ghost_sched%R_recv_req_pers(max(ghost_sched%n_recv_ranks, 1)))
    allocate(ghost_sched%R_send_req_pers(max(ghost_sched%n_send_ranks, 1)))

    ! Initialize persistent receives for R
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank + TAG_OFFSET
      call MPI_RECV_INIT(ghost_sched%T_recv_buf(buf_pos), &
                         ghost_sched%recv_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%recv_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%R_recv_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Initialize persistent sends for R
    do r = 1, ghost_sched%n_send_ranks
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_ + TAG_OFFSET
      call MPI_SEND_INIT(ghost_sched%T_send_buf(buf_pos), &
                         ghost_sched%send_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%send_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%R_send_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ghost_sched%persistent_R_init = .true.
  end subroutine init_persistent_R_requests


  !> Free all persistent MPI requests.
  subroutine cleanup_persistent_requests()
    use mpi
    implicit none
    integer :: ierr, r, fl

    if (allocated(ghost_sched%halo)) then
      do fl = 1, size(ghost_sched%halo)
        if (.not. ghost_sched%halo(fl)%init) cycle
        do r = 1, ghost_sched%n_send_ranks
          call MPI_REQUEST_FREE(ghost_sched%halo(fl)%send_req(r), ierr)
        end do
        do r = 1, ghost_sched%n_recv_ranks
          call MPI_REQUEST_FREE(ghost_sched%halo(fl)%recv_req(r), ierr)
        end do
        ghost_sched%halo(fl)%init = .false.
      end do
    end if

    if (ghost_sched%persistent_R_init) then
      do r = 1, ghost_sched%n_send_ranks
        call MPI_REQUEST_FREE(ghost_sched%R_send_req_pers(r), ierr)
      end do
      do r = 1, ghost_sched%n_recv_ranks
        call MPI_REQUEST_FREE(ghost_sched%R_recv_req_pers(r), ierr)
      end do
      ghost_sched%persistent_R_init = .false.
    end if
  end subroutine cleanup_persistent_requests


  !> Sort ghost entries by remote_rank using insertion sort (lists are small).
  subroutine sort_entries_by_rank(list, n)
    implicit none
    type(ghost_entry_type), intent(inout) :: list(:)
    integer, intent(in) :: n
    ! Local
    integer :: i, j
    type(ghost_entry_type) :: tmp

    do i = 2, n
      tmp = list(i)
      j = i - 1
      do while (j >= 1 .and. list(j)%remote_rank > tmp%remote_rank)
        list(j+1) = list(j)
        j = j - 1
      end do
      list(j+1) = tmp
    end do
  end subroutine sort_entries_by_rank


  !> Build rank groups from a sorted entry list.
  subroutine build_rank_groups(list, n, groups, n_groups)
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(rank_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    ! Local
    integer :: i, ng, cur_rank

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count unique ranks
    ng = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%rank   = list(1)%remote_rank
    groups(1)%offset = 1
    groups(1)%count  = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) then
        ng = ng + 1
        groups(ng)%rank   = list(i)%remote_rank
        groups(ng)%offset = i
        groups(ng)%count  = 1
      else
        groups(ng)%count = groups(ng)%count + 1
      end if
    end do
  end subroutine build_rank_groups


  !> Build face groups from send entry list.
  !> Groups consecutive entries that share (local_block, src_f).
  subroutine build_face_groups_send(list, n, groups, n_groups)
    use FUSS_Parameters_m, only: guide
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(face_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    integer :: i, ng

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count face groups
    ng = 1
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%src_f /= list(i-1)%src_f) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%block = list(1)%local_block
    groups(1)%face  = list(1)%src_f
    groups(1)%entry_start = 1
    groups(1)%entry_count = 1
    groups(1)%di = guide(list(1)%src_f, 1)
    groups(1)%dj = guide(list(1)%src_f, 2)
    groups(1)%dk = guide(list(1)%src_f, 3)
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%src_f /= list(i-1)%src_f) then
        ng = ng + 1
        groups(ng)%block = list(i)%local_block
        groups(ng)%face  = list(i)%src_f
        groups(ng)%entry_start = i
        groups(ng)%entry_count = 1
        groups(ng)%di = guide(list(i)%src_f, 1)
        groups(ng)%dj = guide(list(i)%src_f, 2)
        groups(ng)%dk = guide(list(i)%src_f, 3)
      else
        groups(ng)%entry_count = groups(ng)%entry_count + 1
      end if
    end do
  end subroutine build_face_groups_send


  !> Build face groups from recv entry list.
  !> Groups consecutive entries that share (local_block, dst_f).
  subroutine build_face_groups_recv(list, n, groups, n_groups)
    use FUSS_Parameters_m, only: guide
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(face_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    integer :: i, ng

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count face groups
    ng = 1
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%dst_f /= list(i-1)%dst_f) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%block = list(1)%local_block
    groups(1)%face  = list(1)%dst_f
    groups(1)%entry_start = 1
    groups(1)%entry_count = 1
    groups(1)%di = guide(list(1)%dst_f, 1)
    groups(1)%dj = guide(list(1)%dst_f, 2)
    groups(1)%dk = guide(list(1)%dst_f, 3)
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%dst_f /= list(i-1)%dst_f) then
        ng = ng + 1
        groups(ng)%block = list(i)%local_block
        groups(ng)%face  = list(i)%dst_f
        groups(ng)%entry_start = i
        groups(ng)%entry_count = 1
        groups(ng)%di = guide(list(i)%dst_f, 1)
        groups(ng)%dj = guide(list(i)%dst_f, 2)
        groups(ng)%dk = guide(list(i)%dst_f, 3)
      else
        groups(ng)%entry_count = groups(ng)%entry_count + 1
      end if
    end do
  end subroutine build_face_groups_recv


  !> Build the chimera donor-cell exchange schedule from type-102 BC entries.
  !> Each distinct donor cell travels once per peer rank: entries are sorted by
  !> (peer rank, donor cell) and repeats dropped. Both sides hold the same set
  !> per rank pair, so the same full-key order keeps send and recv matched.
  subroutine build_chimera_schedule(domain)
    use FUSS_Advanced_Types_m
    use mpi, only: MPI_STATUS_SIZE

    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    integer :: i, c, ns, nr, bm, bs, max_ranks

    ! Count send and recv donor references
    ns = 0; nr = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 102) cycle
      if (.not. allocated(domain%bc(i)%donorID)) cycle
      bm = domain%bc(i)%b
      do c = 1, sum(domain%bc(i)%ni)
        bs = domain%bc(i)%donorID(c,1)
        if (is_local_block(bm) .and. (.not. is_local_block(bs))) nr = nr + 1
        if (is_local_block(bs) .and. (.not. is_local_block(bm))) ns = ns + 1
      end do
    end do

    if (allocated(ghost_sched%chim_send_list)) deallocate(ghost_sched%chim_send_list)
    if (allocated(ghost_sched%chim_recv_list)) deallocate(ghost_sched%chim_recv_list)
    if (allocated(ghost_sched%chim_send_groups)) deallocate(ghost_sched%chim_send_groups)
    if (allocated(ghost_sched%chim_recv_groups)) deallocate(ghost_sched%chim_recv_groups)
    if (allocated(ghost_sched%chim_send_buf)) deallocate(ghost_sched%chim_send_buf)
    if (allocated(ghost_sched%chim_recv_buf)) deallocate(ghost_sched%chim_recv_buf)
    if (allocated(ghost_sched%chim_send_req)) deallocate(ghost_sched%chim_send_req)
    if (allocated(ghost_sched%chim_recv_req)) deallocate(ghost_sched%chim_recv_req)
    if (allocated(ghost_sched%chim_stat)) deallocate(ghost_sched%chim_stat)
    allocate(ghost_sched%chim_send_list(ns))
    allocate(ghost_sched%chim_recv_list(nr))

    ! Fill lists
    ns = 0; nr = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 102) cycle
      if (.not. allocated(domain%bc(i)%donorID)) cycle
      bm = domain%bc(i)%b
      do c = 1, sum(domain%bc(i)%ni)
        bs = domain%bc(i)%donorID(c,1)
        if (is_local_block(bm) .and. (.not. is_local_block(bs))) then
          nr = nr + 1
          ghost_sched%chim_recv_list(nr)%remote_rank = block_owner(bs)
          ghost_sched%chim_recv_list(nr)%bs = bs
          ghost_sched%chim_recv_list(nr)%is = domain%bc(i)%donorID(c,2)
          ghost_sched%chim_recv_list(nr)%js = domain%bc(i)%donorID(c,3)
          ghost_sched%chim_recv_list(nr)%ks = domain%bc(i)%donorID(c,4)
        end if
        if (is_local_block(bs) .and. (.not. is_local_block(bm))) then
          ns = ns + 1
          ghost_sched%chim_send_list(ns)%remote_rank = block_owner(bm)
          ghost_sched%chim_send_list(ns)%bs = bs
          ghost_sched%chim_send_list(ns)%is = domain%bc(i)%donorID(c,2)
          ghost_sched%chim_send_list(ns)%js = domain%bc(i)%donorID(c,3)
          ghost_sched%chim_send_list(ns)%ks = domain%bc(i)%donorID(c,4)
        end if
      end do
    end do

    call sort_unique_chim_entries(ghost_sched%chim_send_list, ns)
    call sort_unique_chim_entries(ghost_sched%chim_recv_list, nr)
    ghost_sched%n_chim_send = ns
    ghost_sched%n_chim_recv = nr

    ! Build per-rank groups
    call build_chim_rank_groups(ghost_sched%chim_send_list, ns, &
                                ghost_sched%chim_send_groups, ghost_sched%n_chim_send_ranks)
    call build_chim_rank_groups(ghost_sched%chim_recv_list, nr, &
                                ghost_sched%chim_recv_groups, ghost_sched%n_chim_recv_ranks)

    if (ns > 0 .or. nr > 0) then
      max_ranks = max(ghost_sched%n_chim_send_ranks, ghost_sched%n_chim_recv_ranks, 1)
      allocate(ghost_sched%chim_send_buf(max(ns,1)))
      allocate(ghost_sched%chim_recv_buf(max(nr,1)))
      allocate(ghost_sched%chim_send_req(max(ghost_sched%n_chim_send_ranks,1)))
      allocate(ghost_sched%chim_recv_req(max(ghost_sched%n_chim_recv_ranks,1)))
      allocate(ghost_sched%chim_stat(MPI_STATUS_SIZE, max_ranks))
    end if

    if (mpi_is_root .and. (ns > 0 .or. nr > 0)) then
      write(*,'(A,I0,A,I0,A)') ' Chimera exchange schedule: ', ns, &
        ' donor cells to send, ', nr, ' to receive (rank 0)'
    end if
  end subroutine build_chimera_schedule


  !> Sort chimera entries by (remote_rank, bs, ks, js, is) and remove repeats.
  !> Bottom-up merge sort (the lists hold one entry per (bc, donor) pair); on
  !> return the list is shrunk to the n distinct entries.
  subroutine sort_unique_chim_entries(list, n)
    implicit none
    type(chim_entry_type), allocatable, intent(inout) :: list(:)
    integer, intent(inout) :: n
    ! Local
    type(chim_entry_type), allocatable :: tmp(:)
    integer :: width, lo, mid, hi, i, j, k, m

    if (n <= 1) return

    allocate(tmp(n))
    width = 1
    do while (width < n)
      do lo = 1, n, 2*width
        mid = min(lo + width, n + 1)
        hi  = min(lo + 2*width, n + 1)
        i = lo; j = mid; k = lo
        do while (i < mid .and. j < hi)
          if (chim_key_le(list(i), list(j))) then
            tmp(k) = list(i); i = i + 1
          else
            tmp(k) = list(j); j = j + 1
          end if
          k = k + 1
        end do
        do while (i < mid)
          tmp(k) = list(i); i = i + 1; k = k + 1
        end do
        do while (j < hi)
          tmp(k) = list(j); j = j + 1; k = k + 1
        end do
      end do
      list(1:n) = tmp(1:n)
      width = 2*width
    end do
    deallocate(tmp)

    ! Keep the first of each run of equal keys
    m = 1
    do i = 2, n
      if (chim_key_le(list(i), list(m)) .and. chim_key_le(list(m), list(i))) cycle
      m = m + 1
      list(m) = list(i)
    end do
    n = m
    list = list(1:n)
  end subroutine sort_unique_chim_entries


  !> Lexicographic a <= b on (remote_rank, bs, ks, js, is).
  pure logical function chim_key_le(a, b)
    implicit none
    type(chim_entry_type), intent(in) :: a, b

    if (a%remote_rank /= b%remote_rank) then
      chim_key_le = a%remote_rank < b%remote_rank
    elseif (a%bs /= b%bs) then
      chim_key_le = a%bs < b%bs
    elseif (a%ks /= b%ks) then
      chim_key_le = a%ks < b%ks
    elseif (a%js /= b%js) then
      chim_key_le = a%js < b%js
    else
      chim_key_le = a%is <= b%is
    end if
  end function chim_key_le


  !> Build rank groups from a sorted chimera entry list.
  subroutine build_chim_rank_groups(list, n, groups, n_groups)
    implicit none
    type(chim_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(rank_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    ! Local
    integer :: i, ng

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ng = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%rank   = list(1)%remote_rank
    groups(1)%offset = 1
    groups(1)%count  = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) then
        ng = ng + 1
        groups(ng)%rank   = list(i)%remote_rank
        groups(ng)%offset = i
        groups(ng)%count  = 1
      else
        groups(ng)%count = groups(ng)%count + 1
      end if
    end do
  end subroutine build_chim_rank_groups
#endif


  !> The halo of field-list row fl, split so that packing and unpacking can run over the
  !> whole thread team (plan 11 section 2.2; the T path is the template):
  !>   post_recv (single) -> pack (omp do over face groups) -> post_send (single)
  !>   wait_unpack (single) -> wait_send (single)
  !> A row without a halo on this rank (no block connection to another rank) is a no-op.

  !> Post the persistent receives of row fl. Call from !$omp single.
  subroutine exchange_ghost_post_recv(domain, fl)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in) :: fl
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call halo_post_recv(ghost_sched%halo(fl))
#endif
  end subroutine exchange_ghost_post_recv


  !> Pack face groups fg_start..fg_end of row fl. Call from !$omp do over
  !> 1..ghost_sched%n_send_faces.
  subroutine exchange_ghost_pack(domain, fl, fg_start, fg_end)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    integer, intent(in) :: fl, fg_start, fg_end
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call halo_pack_faces(ghost_sched%halo(fl), domain, fg_start, fg_end)
#endif
  end subroutine exchange_ghost_pack


  !> Post the persistent sends of row fl. Call from !$omp single after the pack.
  subroutine exchange_ghost_post_send(domain, fl)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in) :: fl
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call halo_post_send(ghost_sched%halo(fl))
#endif
  end subroutine exchange_ghost_post_send


  !> Wait for the receives of row fl and unpack them into the ghost cells. Call from
  !> !$omp single.
  subroutine exchange_ghost_wait_unpack(domain, fl)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in) :: fl
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call halo_wait_recv(ghost_sched%halo(fl))
    call halo_unpack_faces(ghost_sched%halo(fl), domain, 1, ghost_sched%n_recv_faces)
#endif
  end subroutine exchange_ghost_wait_unpack


  !> Wait for the sends of row fl. Call from !$omp single.
  subroutine exchange_ghost_wait_send(domain, fl)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in) :: fl
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call halo_wait_send(ghost_sched%halo(fl))
#endif
  end subroutine exchange_ghost_wait_send


  !> The T halo: the calls above on T's row (upstream's names, kept for Fill_Ghost_Cell).
  subroutine exchange_ghost_T_post_recv(domain)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
    call exchange_ghost_post_recv(domain, fl_T)
  end subroutine exchange_ghost_T_post_recv

  subroutine exchange_ghost_T_pack(domain, fg_start, fg_end)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    integer, intent(in) :: fg_start, fg_end
    if (mpi_size_ <= 1) return
    call exchange_ghost_pack(domain, fl_T, fg_start, fg_end)
  end subroutine exchange_ghost_T_pack

  subroutine exchange_ghost_T_post_send(domain)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
    call exchange_ghost_post_send(domain, fl_T)
  end subroutine exchange_ghost_T_post_send

  subroutine exchange_ghost_T_wait_unpack(domain)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
    call exchange_ghost_wait_unpack(domain, fl_T)
  end subroutine exchange_ghost_T_wait_unpack

  subroutine exchange_ghost_T_wait_send(domain)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
    call exchange_ghost_wait_send(domain, fl_T)
  end subroutine exchange_ghost_T_wait_send


  !> Exchange ghost cell R (residual) data between MPI ranks.
  subroutine exchange_ghost_R(domain)
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc
    use FUSS_Parameters_m, only: guide

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain

    if (mpi_size_ <= 1) return
    call set_active_mg_level(domain%mg_level)
#ifdef USE_MPI
    call exchange_R_field_begin(domain)
    call exchange_R_field_end(domain)
#endif
  end subroutine exchange_ghost_R


  !> Gather every remote block's nodes to root (plan 11 stage 2). A moving mesh moves each
  !> block on its own rank only, and root writes every block's geometry: in the solution
  !> file, which is also the restart, and in the wall output. Root's copies of the remote
  !> nodes (kept since setup) are overwritten. Collective: every rank calls it.
  subroutine gather_node_to_root(domain)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer :: b

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, off, n, i, j, k
      integer, allocatable :: reqs(:), stats(:,:), boff(:), bn(:)
      real(R8), allocatable :: buf(:)

      ! root holds the remote blocks, every other rank sends its own
      allocate(boff(domain%nb), bn(domain%nb))
      off = 0
      do b = 1, domain%nb
        bn(b) = 3 * (domain%blk(b)%dim(1)+1) * (domain%blk(b)%dim(2)+1) * (domain%blk(b)%dim(3)+1)
        boff(b) = off
        if (mpi_is_root .neqv. is_local_block(b)) off = off + bn(b)
      end do
      allocate(buf(max(off, 1)), reqs(domain%nb), stats(MPI_STATUS_SIZE, domain%nb))

      nreq = 0
      do b = 1, domain%nb
        if (mpi_is_root .and. .not. is_local_block(b)) then
          nreq = nreq + 1
          call MPI_IRECV(buf(boff(b)+1), bn(b), MPI_DOUBLE_PRECISION, block_owner(b), &
                         NODE_GATHER_TAG + b, MPI_COMM_WORLD, reqs(nreq), ierr)
          call check_mpi_error(ierr)
        else if (.not. mpi_is_root .and. is_local_block(b)) then
          n = boff(b)
          do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2); do i = 0, domain%blk(b)%dim(1)
            buf(n+1:n+3) = domain%blk(b)%node(i,j,k)%c
            n = n + 3
          end do; end do; end do
          nreq = nreq + 1
          call MPI_ISEND(buf(boff(b)+1), bn(b), MPI_DOUBLE_PRECISION, 0, &
                         NODE_GATHER_TAG + b, MPI_COMM_WORLD, reqs(nreq), ierr)
          call check_mpi_error(ierr)
        end if
      end do
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      if (mpi_is_root) then
        do b = 1, domain%nb
          if (is_local_block(b)) cycle
          n = boff(b)
          do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2); do i = 0, domain%blk(b)%dim(1)
            domain%blk(b)%node(i,j,k)%c = buf(n+1:n+3)
            n = n + 3
          end do; end do; end do
        end do
      end if
      deallocate(buf, reqs, stats, boff, bn)
    end block
#endif
  end subroutine gather_node_to_root


  !> The node planes of every block connection between two ranks (plan 11 stage 2).
  !> BC_Connect_Metrics builds a 101 face's ghost metric from the neighbour block's node
  !> planes beyond the shared one (gc of them, Lib_Metrics). With the neighbour on another
  !> rank this rank's copy of those nodes would be the setup geometry, so after every mesh
  !> update each rank sends, for every face of its own blocks that a block of another rank
  !> connects to, the planes 0..gc of that face, and receives its neighbours' into its
  !> copies of their nodes. One message per (block, face, receiving rank); both ends list
  !> the transfers in the same order, from the full bc list every rank holds.
  subroutine exchange_node_planes(domain)
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, t, i, b, f, nreq, n, off, ii, jj, kk
      integer :: lo(3), hi(3)
      logical, allocatable :: seen(:,:,:)
      integer, allocatable :: reqs(:), stats(:,:), toff(:), tn(:)
      real(R8), allocatable :: buf(:)

      if (.not. np_built) then
        allocate(seen(6, domain%nb, 0:mpi_size_-1))
        seen = .false.
        do i = 1, domain%nbound
          if (domain%bc(i)%type /= 101) cycle
          b = domain%bc(i)%bs
          f = domain%bc(i)%fs
          if (block_owner(b) == block_owner(domain%bc(i)%b)) cycle
          seen(f, b, block_owner(domain%bc(i)%b)) = .true.
        end do
        np_n = count(seen)
        allocate(np_block(np_n), np_face(np_n), np_to(np_n))
        t = 0
        do b = 1, domain%nb; do f = 1, 6; do i = 0, mpi_size_-1
          if (.not. seen(f, b, i)) cycle
          t = t + 1
          np_block(t) = b; np_face(t) = f; np_to(t) = i
        end do; end do; end do
        deallocate(seen)
        np_built = .true.
      end if
      if (np_n == 0) return

      allocate(toff(np_n), tn(np_n), reqs(np_n), stats(MPI_STATUS_SIZE, np_n))
      off = 0
      do t = 1, np_n
        call plane_range(np_block(t), np_face(t), lo, hi)
        tn(t) = 3 * (hi(1)-lo(1)+1) * (hi(2)-lo(2)+1) * (hi(3)-lo(3)+1)
        toff(t) = off
        if (block_owner(np_block(t)) == mpi_rank_ .or. np_to(t) == mpi_rank_) off = off + tn(t)
      end do
      allocate(buf(max(off, 1)))

      nreq = 0
      do t = 1, np_n
        b = np_block(t)
        if (np_to(t) == mpi_rank_) then
          nreq = nreq + 1
          call MPI_IRECV(buf(toff(t)+1), tn(t), MPI_DOUBLE_PRECISION, block_owner(b), &
                         NODE_PLANE_TAG + 6*(b-1) + np_face(t), MPI_COMM_WORLD, reqs(nreq), ierr)
          call check_mpi_error(ierr)
        else if (block_owner(b) == mpi_rank_) then
          call plane_range(b, np_face(t), lo, hi)
          n = toff(t)
          do kk = lo(3), hi(3); do jj = lo(2), hi(2); do ii = lo(1), hi(1)
            buf(n+1:n+3) = domain%blk(b)%node(ii,jj,kk)%c
            n = n + 3
          end do; end do; end do
          nreq = nreq + 1
          call MPI_ISEND(buf(toff(t)+1), tn(t), MPI_DOUBLE_PRECISION, np_to(t), &
                         NODE_PLANE_TAG + 6*(b-1) + np_face(t), MPI_COMM_WORLD, reqs(nreq), ierr)
          call check_mpi_error(ierr)
        end if
      end do
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      do t = 1, np_n
        if (np_to(t) /= mpi_rank_) cycle
        b = np_block(t)
        call plane_range(b, np_face(t), lo, hi)
        n = toff(t)
        do kk = lo(3), hi(3); do jj = lo(2), hi(2); do ii = lo(1), hi(1)
          domain%blk(b)%node(ii,jj,kk)%c = buf(n+1:n+3)
          n = n + 3
        end do; end do; end do
      end do
      deallocate(buf, toff, tn, reqs, stats)
    end block
#endif

  contains

    !> The node index ranges of the planes 0..gc of face f of block b (clamped to the block)
    subroutine plane_range(b, f, lo, hi)
      integer, intent(in)  :: b, f
      integer, intent(out) :: lo(3), hi(3)
      integer :: d
      lo = 0
      hi = domain%blk(b)%dim
      d = (f + 1) / 2
      if (mod(f, 2) == 1) then
        hi(d) = min(gc, domain%blk(b)%dim(d))
      else
        lo(d) = max(0, domain%blk(b)%dim(d) - gc)
      end if
    end subroutine plane_range

  end subroutine exchange_node_planes


  !> Gather the interior cells of one field (the array with code `code`, ncomp values per
  !> cell) of every block to root, for the solution file (plan 11 stage 3: the field-list
  !> output rows and the state the auxiliary columns are computed from). Root's copy of a
  !> remote block's array, freed at setup, is allocated first. Exact copies, one message
  !> per block. Collective: every rank calls it, outside any parallel region.
  subroutine gather_field_to_root(domain, code, ncomp)
    use FUSS_Advanced_Types_m
    use FUSS_Field_List_m, only: Field_Cell_Get, Field_Cell_Set, Field_Allocate
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in) :: code, ncomp
    integer :: b

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    if (code == 0) then
      write(*,'(A)') ' [ERROR] gather_field_to_root: an output field without an accessor code'
      call fuss_abort(1)
    end if
    block
      use mpi
      integer :: ierr, nreq, off, n, i, j, k
      integer, allocatable :: reqs(:), stats(:,:), boff(:), bn(:)
      real(R8), allocatable :: buf(:)

      ! root holds the remote blocks, every other rank sends its own
      allocate(boff(domain%nb), bn(domain%nb))
      off = 0
      do b = 1, domain%nb
        bn(b) = ncomp * domain%blk(b)%dim(1) * domain%blk(b)%dim(2) * domain%blk(b)%dim(3)
        boff(b) = off
        if (mpi_is_root .neqv. is_local_block(b)) off = off + bn(b)
      end do
      allocate(buf(max(off, 1)), reqs(domain%nb), stats(MPI_STATUS_SIZE, domain%nb))

      nreq = 0
      do b = 1, domain%nb
        if (mpi_is_root .and. .not. is_local_block(b)) then
          nreq = nreq + 1
          call MPI_IRECV(buf(boff(b)+1), bn(b), MPI_DOUBLE_PRECISION, block_owner(b), &
                         FIELD_GATHER_TAG + b, MPI_COMM_WORLD, reqs(nreq), ierr)
          call check_mpi_error(ierr)
        else if (.not. mpi_is_root .and. is_local_block(b)) then
          n = boff(b)
          do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
            call Field_Cell_Get(domain%blk(b), code, i, j, k, buf(n+1:n+ncomp))
            n = n + ncomp
          end do; end do; end do
          nreq = nreq + 1
          call MPI_ISEND(buf(boff(b)+1), bn(b), MPI_DOUBLE_PRECISION, 0, &
                         FIELD_GATHER_TAG + b, MPI_COMM_WORLD, reqs(nreq), ierr)
          call check_mpi_error(ierr)
        end if
      end do
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      if (mpi_is_root) then
        do b = 1, domain%nb
          if (is_local_block(b)) cycle
          call Field_Allocate(domain%blk(b), code, ncomp)
          n = boff(b)
          do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
            call Field_Cell_Set(domain%blk(b), code, i, j, k, buf(n+1:n+ncomp))
            n = n + ncomp
          end do; end do; end do
        end do
      end if
      deallocate(buf, reqs, stats, boff, bn)
    end block
#endif
  end subroutine gather_field_to_root


  !> Gather all blocks' T interior data to rank 0 for I/O.
  !> Uses non-blocking MPI_ISEND/IRECV + WAITALL for overlap.
  !> Temporarily allocates T on root for remote blocks if not already allocated.
  subroutine gather_T_to_root(domain)
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, n1, n2, n3, ncells

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, r, total_buf_size, offset
      integer, allocatable :: reqs(:), blk_offset(:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      ! Compute per-block sizes and total buffer size
      allocate(blk_offset(domain%nb), blk_ncells(domain%nb))
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        n3 = domain%blk(b)%dim(3)
        blk_ncells(b) = n1 * n2 * n3
        blk_offset(b) = total_buf_size
        if (mpi_is_root) then
          if (.not. is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        else
          if (is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        end if
      end do
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      ! Recompute offsets contiguously
      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_ncells(b)
        else
          if (is_local_block(b)) offset = offset + blk_ncells(b)
        end if
      end do

      nreq = 0
      if (mpi_is_root) then
        ! Root: post all receives for remote blocks
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: pack local blocks and post all sends
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)) = &
              reshape(domain%blk(b)%T(1:n1, 1:n2, 1:n3), [blk_ncells(b)])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Root: unpack received data (allocate T if it was previously freed)
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            if (.not. allocated(domain%blk(b)%T)) &
              allocate(domain%blk(b)%T(1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc))
            domain%blk(b)%T(1:n1, 1:n2, 1:n3) = &
              reshape(buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)), [n1, n2, n3])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_ncells)
    end block
#endif
  end subroutine gather_T_to_root


  !> Gather all blocks' T interior data to rank 0 for I/O.
  !> Uses non-blocking MPI_ISEND/IRECV + WAITALL for overlap.
  !> Temporarily allocates T on root for remote blocks if not already allocated.
  subroutine gather_matID_to_root(domain)
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, n1, n2, n3, ncells

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, r, total_buf_size, offset
      integer, allocatable :: reqs(:), blk_offset(:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      ! Compute per-block sizes and total buffer size
      allocate(blk_offset(domain%nb), blk_ncells(domain%nb))
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        n3 = domain%blk(b)%dim(3)
        blk_ncells(b) = n1 * n2 * n3
        blk_offset(b) = total_buf_size
        if (mpi_is_root) then
          if (.not. is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        else
          if (is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        end if
      end do
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      ! Recompute offsets contiguously
      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_ncells(b)
        else
          if (is_local_block(b)) offset = offset + blk_ncells(b)
        end if
      end do

      nreq = 0
      if (mpi_is_root) then
        ! Root: post all receives for remote blocks
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: pack local blocks and post all sends
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)) = &
              reshape(domain%blk(b)%matID(1:n1, 1:n2, 1:n3), [blk_ncells(b)])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Root: unpack received data (allocate T if it was previously freed)
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            if (.not. allocated(domain%blk(b)%matID)) &
              allocate(domain%blk(b)%matID(1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc))
            domain%blk(b)%matID(1:n1, 1:n2, 1:n3) = &
              reshape(buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)), [n1, n2, n3])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_ncells)
    end block
#endif
  end subroutine gather_matID_to_root


  !> Gather all blocks' T interior data to rank 0 for I/O.
  !> Uses non-blocking MPI_ISEND/IRECV + WAITALL for overlap.
  !> Temporarily allocates T on root for remote blocks if not already allocated.
  subroutine gather_qvol_to_root(domain)
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, n1, n2, n3, ncells

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, r, total_buf_size, offset
      integer, allocatable :: reqs(:), blk_offset(:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      ! Compute per-block sizes and total buffer size
      allocate(blk_offset(domain%nb), blk_ncells(domain%nb))
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        n3 = domain%blk(b)%dim(3)
        blk_ncells(b) = n1 * n2 * n3
        blk_offset(b) = total_buf_size
        if (mpi_is_root) then
          if (.not. is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        else
          if (is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        end if
      end do
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      ! Recompute offsets contiguously
      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_ncells(b)
        else
          if (is_local_block(b)) offset = offset + blk_ncells(b)
        end if
      end do

      nreq = 0
      if (mpi_is_root) then
        ! Root: post all receives for remote blocks
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: pack local blocks and post all sends
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)) = &
              reshape(domain%blk(b)%qvol(1:n1, 1:n2, 1:n3), [blk_ncells(b)])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Root: unpack received data (allocate T if it was previously freed)
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            if (.not. allocated(domain%blk(b)%qvol)) &
              allocate(domain%blk(b)%qvol(1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc))
            domain%blk(b)%qvol(1:n1, 1:n2, 1:n3) = &
              reshape(buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)), [n1, n2, n3])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_ncells)
    end block
#endif
  end subroutine gather_qvol_to_root


  !> Gather diagnostic arrays (R, dtlocal, beta) from all ranks to root.
  !> On root, remote blocks may have had these arrays deallocated by
  !> deallocate_remote_computation_data, so they are re-allocated here
  !> to receive the data.  All ranks must call this (collective).
  subroutine gather_diagnostic_to_root(domain)
    use FUSS_Advanced_Types_m

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer :: b, n1, n2, n3

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, offset, total_buf_size
      integer :: ncells_R, ncells_dt
      integer, allocatable :: reqs(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)
      ! Per-block: n1*n2*n3 for R, n1*n2*n3 for dtlocal
      integer, allocatable :: blk_offset(:), blk_size(:)

      allocate(blk_offset(domain%nb), blk_size(domain%nb))

      ! Compute buffer sizes
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        n3 = domain%blk(b)%dim(3)
        blk_size(b) = n1 * n2 * n3 + n1 * n2 * n3  ! R + dtlocal
      end do

      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_size(b)
        else
          if (is_local_block(b)) offset = offset + blk_size(b)
        end if
      end do
      total_buf_size = offset
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      ! On root: ensure R, dtlocal are allocated for remote blocks
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            if (.not. allocated(domain%blk(b)%R)) &
              allocate(domain%blk(b)%R(n1, n2, n3))
            if (.not. allocated(domain%blk(b)%dtlocal)) &
              allocate(domain%blk(b)%dtlocal(n1, n2, n3))
          end if
        end do
      end if

      nreq = 0
      if (mpi_is_root) then
        ! Root: post receives for remote blocks
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_size(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), 1000+b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: pack local blocks and send
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            ncells_R  = n1 * n2 * n3
            ncells_dt = n1 * n2 * n3
            offset = blk_offset(b)
            buf(offset+1 : offset+ncells_R) = &
              reshape(domain%blk(b)%R(1:n1, 1:n2, 1:n3), [ncells_R])
            offset = offset + ncells_R
            buf(offset+1 : offset+ncells_dt) = &
              reshape(domain%blk(b)%dtlocal(1:n1, 1:n2, 1:n3), [ncells_dt])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_size(b), MPI_DOUBLE_PRECISION, &
                           0, 1000+b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Root: unpack received data
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            ncells_R  = n1 * n2 * n3
            ncells_dt = n1 * n2 * n3
            offset = blk_offset(b)
            domain%blk(b)%R(1:n1, 1:n2, 1:n3) = &
              reshape(buf(offset+1 : offset+ncells_R), [n1, n2, n3])
            offset = offset + ncells_R
            domain%blk(b)%dtlocal(1:n1, 1:n2, 1:n3) = &
              reshape(buf(offset+1 : offset+ncells_dt), [n1, n2, n3])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_size)
    end block
#endif
  end subroutine gather_diagnostic_to_root


  !> Scatter blocks' T data from rank 0 to owning ranks after reading.
  !> Uses non-blocking MPI_ISEND/IRECV + WAITALL for overlap.
  subroutine scatter_T_from_root(domain)
    use FUSS_Advanced_Types_m

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, n1, n2, n3, ncells

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, r, total_buf_size, offset
      integer, allocatable :: reqs(:), blk_offset(:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      ! Compute per-block sizes and total buffer size
      allocate(blk_offset(domain%nb), blk_ncells(domain%nb))
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        n3 = domain%blk(b)%dim(3)
        blk_ncells(b) = n1 * n2 * n3
      end do

      ! Compute contiguous buffer offsets
      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_ncells(b)
        else
          if (is_local_block(b)) offset = offset + blk_ncells(b)
        end if
      end do
      total_buf_size = offset
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      nreq = 0
      if (mpi_is_root) then
        ! Root: pack remote blocks and post all sends
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)) = &
              reshape(domain%blk(b)%T(1:n1, 1:n2, 1:n3), [blk_ncells(b)])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: post all receives for local blocks
        do b = 1, domain%nb
          if (is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Non-root: unpack received data
      if (.not. mpi_is_root) then
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            domain%blk(b)%T(1:n1, 1:n2, 1:n3) = &
              reshape(buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)), [n1, n2, n3])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_ncells)
    end block
#endif
  end subroutine scatter_T_from_root


  !> Exchange Tg(3:6) stencil data for inter-rank viscous connection BCs.
  !> Messages are aggregated per remote rank.
  subroutine exchange_ghost_Tg(domain)
    use FUSS_Advanced_Types_m

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_Tg_field(domain, 1)
#endif
  end subroutine exchange_ghost_Tg


  !> The gas-pressure stencil pg(3:6) of the inter-rank connections (plan 11 stage 3), as
  !> Tg(3:6): the source block's tangential pressures, computed on its rank (Fill_Ghost_p).
  subroutine exchange_ghost_pg(domain)
    use FUSS_Advanced_Types_m

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_Tg_field(domain, 2)
#endif
  end subroutine exchange_ghost_pg


  !> The entries of the active schedule (plan 11 stage 3): a send entry is a block
  !> connection whose source block lives here and whose receiving block on another rank, a
  !> recv entry the converse. Both sides list them in bc order, so send entry e to a rank
  !> pairs with that rank's recv entry e from here. 0 on one rank.
  integer function entry_send_count() result(n)
    n = 0
    if (mpi_size_ <= 1) return
    n = ghost_sched%n_send
  end function entry_send_count

  integer function entry_recv_count() result(n)
    n = 0
    if (mpi_size_ <= 1) return
    n = ghost_sched%n_recv
  end function entry_recv_count

  !> The source cell of send entry e: block b (on this rank), cell (i, j, k) next to the face.
  subroutine entry_send_cell(e, b, i, j, k)
    integer, intent(in)  :: e
    integer, intent(out) :: b, i, j, k
    b = ghost_sched%send_list(e)%local_block
    i = ghost_sched%send_list(e)%src_i
    j = ghost_sched%send_list(e)%src_j
    k = ghost_sched%send_list(e)%src_k
  end subroutine entry_send_cell

  !> The source cell of recv entry e: block b (on another rank), cell (i, j, k), and the bc
  !> index ii of the connection (its receiving block lives here).
  subroutine entry_recv_cell(e, b, i, j, k, ii)
    integer, intent(in)  :: e
    integer, intent(out) :: b, i, j, k, ii
    b  = ghost_sched%recv_list(e)%remote_block
    i  = ghost_sched%recv_list(e)%src_i
    j  = ghost_sched%recv_list(e)%src_j
    k  = ghost_sched%recv_list(e)%src_k
    ii = ghost_sched%recv_list(e)%bc_idx
  end subroutine entry_recv_cell

  !> nv values per entry (plan 11 stage 3): sbuf(:, e) of every send entry e travels to the
  !> rank of its receiving block, into rbuf(:, e) of the paired recv entry there. One
  !> message per rank pair, tagged tag_base + the sender's rank. The caller packs and
  !> unpacks (entry_send_cell / entry_recv_cell). Master thread only (FUNNELED); every rank
  !> that has entries calls it with the same tag_base.
  subroutine exchange_entry_values(nv, sbuf, rbuf, tag_base)
    integer,  intent(in)    :: nv, tag_base
    real(R8), intent(in)    :: sbuf(nv,*)
    real(R8), intent(inout) :: rbuf(nv,*)

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: r, ierr
      integer, allocatable :: sreq(:), rreq(:), stats(:,:)
      allocate(sreq(max(ghost_sched%n_send_ranks,1)), rreq(max(ghost_sched%n_recv_ranks,1)), &
               stats(MPI_STATUS_SIZE, max(ghost_sched%n_send_ranks, ghost_sched%n_recv_ranks, 1)))
      do r = 1, ghost_sched%n_recv_ranks
        call MPI_IRECV(rbuf(1, ghost_sched%recv_groups(r)%offset), nv * ghost_sched%recv_groups(r)%count, &
                       MPI_DOUBLE_PRECISION, ghost_sched%recv_groups(r)%rank, &
                       tag_base + ghost_sched%recv_groups(r)%rank, MPI_COMM_WORLD, rreq(r), ierr)
        call check_mpi_error(ierr)
      end do
      do r = 1, ghost_sched%n_send_ranks
        call MPI_ISEND(sbuf(1, ghost_sched%send_groups(r)%offset), nv * ghost_sched%send_groups(r)%count, &
                       MPI_DOUBLE_PRECISION, ghost_sched%send_groups(r)%rank, &
                       tag_base + mpi_rank_, MPI_COMM_WORLD, sreq(r), ierr)
        call check_mpi_error(ierr)
      end do
      if (ghost_sched%n_recv_ranks > 0) then
        call MPI_WAITALL(ghost_sched%n_recv_ranks, rreq, stats, ierr)
        call check_mpi_error(ierr)
      end if
      if (ghost_sched%n_send_ranks > 0) then
        call MPI_WAITALL(ghost_sched%n_send_ranks, sreq, stats, ierr)
        call check_mpi_error(ierr)
      end if
      deallocate(sreq, rreq, stats)
    end block
#endif
  end subroutine exchange_entry_values


  !> Chimera donor-cell exchange, split like the T exchange so that packing
  !> and unpacking run over the whole thread team:
  !>   post_recv (single) -> pack (omp do) -> post_send (single)
  !>   wait_recv (single) -> unpack (omp do) -> wait_send (single)
  !> Unpacked cells land in the kept-alive T arrays of remote donor blocks,
  !> where Ghost_Chimera reads them.

  !> Post chimera receives. Call from !$omp single.
  subroutine exchange_ghost_chimera_post_recv(domain)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call chimera_post_recv()
#endif
  end subroutine exchange_ghost_chimera_post_recv


  !> Pack chimera send entries i_start..i_end. Call from !$omp do over
  !> 1..ghost_sched%n_chim_send.
  subroutine exchange_ghost_chimera_pack(domain, i_start, i_end)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    integer, intent(in) :: i_start, i_end
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call chimera_pack(domain, i_start, i_end)
#endif
  end subroutine exchange_ghost_chimera_pack


  !> Post chimera sends. Call from !$omp single after packing is complete.
  subroutine exchange_ghost_chimera_post_send(domain)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call chimera_post_send()
#endif
  end subroutine exchange_ghost_chimera_post_send


  !> Wait for chimera receives. Call from !$omp single.
  subroutine exchange_ghost_chimera_wait_recv(domain)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call chimera_wait_recv()
#endif
  end subroutine exchange_ghost_chimera_wait_recv


  !> Unpack chimera recv entries i_start..i_end into remote donor blocks' T.
  !> Call from !$omp do over 1..ghost_sched%n_chim_recv. Entries are
  !> distinct cells, so threads never write the same location.
  subroutine exchange_ghost_chimera_unpack(domain, i_start, i_end)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in) :: i_start, i_end
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call chimera_unpack(domain, i_start, i_end)
#endif
  end subroutine exchange_ghost_chimera_unpack


  !> Wait for chimera sends before the send buffer is reused. Call from !$omp single.
  subroutine exchange_ghost_chimera_wait_send(domain)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call chimera_wait_send()
#endif
  end subroutine exchange_ghost_chimera_wait_send


#ifdef USE_MPI
  !> One-shot exchange of matID into ghost cells. Called once during setup
  !> since matID is static. Reuses T_send_buf/T_recv_buf and send_req/recv_req.
  subroutine exchange_matID_once(domain)
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc
    use mpi

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: ns, nr, r, i, g, b, fg, ierr, tag, buf_pos
    integer :: Is, Js, Ks, Id, Jd, Kd, di, dj, dk
    integer :: entry_size
    integer, parameter :: TAG_OFFSET = 40000

    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv
    if (ns == 0 .and. nr == 0) return

    entry_size = gc

    ! Post one receive per remote rank (aggregated)
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank + TAG_OFFSET
      call MPI_IRECV(ghost_sched%T_recv_buf(buf_pos), &
                     ghost_sched%recv_groups(r)%count * entry_size, &
                     MPI_DOUBLE_PRECISION, &
                     ghost_sched%recv_groups(r)%rank, tag, &
                     MPI_COMM_WORLD, ghost_sched%recv_req(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Pack matID into send buffer using face groups (same indexing as T pack)
    do fg = 1, ghost_sched%n_send_faces
      b  = ghost_sched%send_faces(fg)%block
      di = ghost_sched%send_faces(fg)%di
      dj = ghost_sched%send_faces(fg)%dj
      dk = ghost_sched%send_faces(fg)%dk
      do i = ghost_sched%send_faces(fg)%entry_start, &
             ghost_sched%send_faces(fg)%entry_start + ghost_sched%send_faces(fg)%entry_count - 1
        Is = ghost_sched%send_list(i)%src_i
        Js = ghost_sched%send_list(i)%src_j
        Ks = ghost_sched%send_list(i)%src_k
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          ghost_sched%T_send_buf(buf_pos + g) = &
            domain%blk(b)%matID(Is+di*(g-1), Js+dj*(g-1), Ks+dk*(g-1))
        end do
      end do
    end do

    ! Post one send per remote rank
    do r = 1, ghost_sched%n_send_ranks
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_ + TAG_OFFSET
      call MPI_ISEND(ghost_sched%T_send_buf(buf_pos), &
                     ghost_sched%send_groups(r)%count * entry_size, &
                     MPI_DOUBLE_PRECISION, &
                     ghost_sched%send_groups(r)%rank, tag, &
                     MPI_COMM_WORLD, ghost_sched%send_req(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Wait for all receives
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%recv_req, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if

    ! Unpack into ghost cells of blk%matID (same indexing as T unpack)
    do fg = 1, ghost_sched%n_recv_faces
      b  = ghost_sched%recv_faces(fg)%block
      di = ghost_sched%recv_faces(fg)%di
      dj = ghost_sched%recv_faces(fg)%dj
      dk = ghost_sched%recv_faces(fg)%dk
      do i = ghost_sched%recv_faces(fg)%entry_start, &
             ghost_sched%recv_faces(fg)%entry_start + ghost_sched%recv_faces(fg)%entry_count - 1
        Id = ghost_sched%recv_list(i)%dst_i
        Jd = ghost_sched%recv_list(i)%dst_j
        Kd = ghost_sched%recv_list(i)%dst_k
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          domain%blk(b)%matID(Id-di*g, Jd-dj*g, Kd-dk*g) = &
            ghost_sched%T_recv_buf(buf_pos + g)
        end do
      end do
    end do

    ! Wait for all sends to complete
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%send_req, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_matID_once

  
  !> Internal: exchange a tangential stencil using aggregated messaging: which = 1 the
  !> temperature's Tg(3:6), 2 the gas pressure's pg(3:6) (plan 11 stage 3).
  subroutine exchange_Tg_field(domain, which)
    use FUSS_Advanced_Types_m
    use mpi

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in) :: which
    ! Local
    integer :: ns, nr, i, r, ierr, tag, bc_idx, buf_pos, tag_off
    integer :: entry_size
    integer, parameter :: TG_NCOLS = 4   ! Tg columns 3..6

    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv
    if (ns == 0 .and. nr == 0) return

    entry_size = TG_NCOLS
    tag_off = 10000                       ! Tg
    if (which == 2) tag_off = PG_TAG      ! pg

    ! Post one receive per remote rank (aggregated)
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank + tag_off  ! unique tag per rank and stencil
      call MPI_IRECV(ghost_sched%Tg_recv_buf(buf_pos), &
                     ghost_sched%recv_groups(r)%count * entry_size, &
                     MPI_DOUBLE_PRECISION, &
                     ghost_sched%recv_groups(r)%rank, tag, &
                     MPI_COMM_WORLD, ghost_sched%recv_req(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Pack all entries per rank into contiguous buffer, then one send per rank
    do r = 1, ghost_sched%n_send_ranks
      do i = ghost_sched%send_groups(r)%offset, &
             ghost_sched%send_groups(r)%offset + ghost_sched%send_groups(r)%count - 1
        bc_idx = ghost_sched%send_list(i)%bc_idx
        buf_pos = (i - 1) * entry_size
        if (which == 2) then
          ghost_sched%Tg_send_buf(buf_pos+1:buf_pos+4) = domain%bc(bc_idx)%pg(3:6)
        else
          ghost_sched%Tg_send_buf(buf_pos+1) = domain%bc(bc_idx)%Tg(3)
          ghost_sched%Tg_send_buf(buf_pos+2) = domain%bc(bc_idx)%Tg(4)
          ghost_sched%Tg_send_buf(buf_pos+3) = domain%bc(bc_idx)%Tg(5)
          ghost_sched%Tg_send_buf(buf_pos+4) = domain%bc(bc_idx)%Tg(6)
        end if
      end do
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_ + tag_off  ! sender's rank as tag
      call MPI_ISEND(ghost_sched%Tg_send_buf(buf_pos), &
                     ghost_sched%send_groups(r)%count * entry_size, &
                     MPI_DOUBLE_PRECISION, &
                     ghost_sched%send_groups(r)%rank, tag, &
                     MPI_COMM_WORLD, ghost_sched%send_req(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Wait for all receives
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%recv_req, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if

    ! Unpack into bc%Tg(3:6) (bc%pg(3:6))
    do i = 1, nr
      bc_idx = ghost_sched%recv_list(i)%bc_idx
      buf_pos = (i - 1) * entry_size
      if (which == 2) then
        domain%bc(bc_idx)%pg(3:6) = ghost_sched%Tg_recv_buf(buf_pos+1:buf_pos+4)
        cycle
      end if
      domain%bc(bc_idx)%Tg(3) = ghost_sched%Tg_recv_buf(buf_pos+1)
      domain%bc(bc_idx)%Tg(4) = ghost_sched%Tg_recv_buf(buf_pos+2)
      domain%bc(bc_idx)%Tg(5) = ghost_sched%Tg_recv_buf(buf_pos+3)
      domain%bc(bc_idx)%Tg(6) = ghost_sched%Tg_recv_buf(buf_pos+4)
    end do

    ! Wait for all sends to complete
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%send_req, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_Tg_field


  !> Internal: post one chimera receive per donor-owner rank (aggregated).
  subroutine chimera_post_recv()
    use mpi
    implicit none
    ! Local
    integer :: r, ierr, tag, buf_pos
    integer, parameter :: TAG_OFFSET = 20000

    if (ghost_sched%n_chim_send == 0 .and. ghost_sched%n_chim_recv == 0) return

    do r = 1, ghost_sched%n_chim_recv_ranks
      buf_pos = ghost_sched%chim_recv_groups(r)%offset
      tag = ghost_sched%chim_recv_groups(r)%rank + TAG_OFFSET
      call MPI_IRECV(ghost_sched%chim_recv_buf(buf_pos), &
                     ghost_sched%chim_recv_groups(r)%count, &
                     MPI_DOUBLE_PRECISION, &
                     ghost_sched%chim_recv_groups(r)%rank, tag, &
                     MPI_COMM_WORLD, ghost_sched%chim_recv_req(r), ierr)
      call check_mpi_error(ierr)
    end do
  end subroutine chimera_post_recv


  !> Internal: pack local donor cells i_start..i_end. Safe from !$omp do.
  subroutine chimera_pack(domain, i_start, i_end)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    integer, intent(in) :: i_start, i_end
    ! Local
    integer :: i

    do i = i_start, i_end
      ghost_sched%chim_send_buf(i) = &
        domain%blk(ghost_sched%chim_send_list(i)%bs)%T( &
          ghost_sched%chim_send_list(i)%is, &
          ghost_sched%chim_send_list(i)%js, &
          ghost_sched%chim_send_list(i)%ks)
    end do
  end subroutine chimera_pack


  !> Internal: post one chimera send per receiver rank (aggregated).
  subroutine chimera_post_send()
    use mpi
    implicit none
    ! Local
    integer :: r, ierr, tag, buf_pos
    integer, parameter :: TAG_OFFSET = 20000

    if (ghost_sched%n_chim_send == 0 .and. ghost_sched%n_chim_recv == 0) return

    do r = 1, ghost_sched%n_chim_send_ranks
      buf_pos = ghost_sched%chim_send_groups(r)%offset
      tag = mpi_rank_ + TAG_OFFSET
      call MPI_ISEND(ghost_sched%chim_send_buf(buf_pos), &
                     ghost_sched%chim_send_groups(r)%count, &
                     MPI_DOUBLE_PRECISION, &
                     ghost_sched%chim_send_groups(r)%rank, tag, &
                     MPI_COMM_WORLD, ghost_sched%chim_send_req(r), ierr)
      call check_mpi_error(ierr)
    end do
  end subroutine chimera_post_send


  !> Internal: wait for all chimera receives.
  subroutine chimera_wait_recv()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_chim_send == 0 .and. ghost_sched%n_chim_recv == 0) return

    if (ghost_sched%n_chim_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_chim_recv_ranks, ghost_sched%chim_recv_req, &
                       ghost_sched%chim_stat(:,1:ghost_sched%n_chim_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine chimera_wait_recv


  !> Internal: unpack donor cells i_start..i_end into the remote donor blocks'
  !> T arrays (kept allocated by deallocate_remote_computation_data). Safe
  !> from !$omp do.
  subroutine chimera_unpack(domain, i_start, i_end)
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in) :: i_start, i_end
    ! Local
    integer :: i

    do i = i_start, i_end
      domain%blk(ghost_sched%chim_recv_list(i)%bs)%T( &
        ghost_sched%chim_recv_list(i)%is, &
        ghost_sched%chim_recv_list(i)%js, &
        ghost_sched%chim_recv_list(i)%ks) = ghost_sched%chim_recv_buf(i)
    end do
  end subroutine chimera_unpack


  !> Internal: wait for all chimera sends.
  subroutine chimera_wait_send()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_chim_send == 0 .and. ghost_sched%n_chim_recv == 0) return

    if (ghost_sched%n_chim_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_chim_send_ranks, ghost_sched%chim_send_req, &
                       ghost_sched%chim_stat(:,1:ghost_sched%n_chim_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine chimera_wait_send


  !> Start the persistent receives of a halo. Call from !$omp single.
  subroutine halo_post_recv(h)
    use mpi
    implicit none
    type(halo_field_type), intent(inout) :: h
    integer :: ierr

    if (.not. h%init) return
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_recv_ranks, h%recv_req, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine halo_post_recv


  !> Pack a halo's send buffer: the gc cells behind every sending face, ncomp values each
  !> (Field_Cell_Get). Safe to call from !$omp do over face groups.
  subroutine halo_pack_faces(h, domain, fg_start, fg_end)
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc
    use FUSS_Field_List_m, only: Field_Cell_Get
    implicit none
    type(halo_field_type), intent(inout) :: h
    type(FUSS_domain_type), intent(in) :: domain
    integer, intent(in) :: fg_start, fg_end
    ! Local
    integer :: fg, i, g, b, buf_pos, nc
    integer :: Is, Js, Ks, di, dj, dk

    if (.not. h%init) return
    nc = h%ncomp

    do fg = fg_start, fg_end
      b  = ghost_sched%send_faces(fg)%block
      di = ghost_sched%send_faces(fg)%di
      dj = ghost_sched%send_faces(fg)%dj
      dk = ghost_sched%send_faces(fg)%dk
      do i = ghost_sched%send_faces(fg)%entry_start, &
             ghost_sched%send_faces(fg)%entry_start + ghost_sched%send_faces(fg)%entry_count - 1
        Is = ghost_sched%send_list(i)%src_i
        Js = ghost_sched%send_list(i)%src_j
        Ks = ghost_sched%send_list(i)%src_k
        buf_pos = (i - 1) * gc * nc
        do g = 1, gc
          call Field_Cell_Get(domain%blk(b), h%code, Is+di*(g-1), Js+dj*(g-1), Ks+dk*(g-1), &
                              h%send_buf(buf_pos + (g-1)*nc + 1 : buf_pos + g*nc))
        end do
      end do
    end do
  end subroutine halo_pack_faces


  !> Start the persistent sends of a halo. Call from !$omp single after the pack.
  subroutine halo_post_send(h)
    use mpi
    implicit none
    type(halo_field_type), intent(inout) :: h
    integer :: ierr

    if (.not. h%init) return
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_send_ranks, h%send_req, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine halo_post_send


  !> Wait for a halo's receives. Call from !$omp single.
  subroutine halo_wait_recv(h)
    use mpi
    implicit none
    type(halo_field_type), intent(inout) :: h
    integer :: ierr

    if (.not. h%init) return
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, h%recv_req, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine halo_wait_recv


  !> Unpack a halo's receive buffer into the gc ghost cells of every receiving face
  !> (Field_Cell_Set). Safe to call from !$omp do over face groups.
  subroutine halo_unpack_faces(h, domain, fg_start, fg_end)
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc
    use FUSS_Field_List_m, only: Field_Cell_Set
    implicit none
    type(halo_field_type), intent(inout) :: h
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in) :: fg_start, fg_end
    ! Local
    integer :: fg, i, g, b, buf_pos, nc
    integer :: Id, Jd, Kd, di, dj, dk

    if (.not. h%init) return
    nc = h%ncomp

    do fg = fg_start, fg_end
      b  = ghost_sched%recv_faces(fg)%block
      di = ghost_sched%recv_faces(fg)%di
      dj = ghost_sched%recv_faces(fg)%dj
      dk = ghost_sched%recv_faces(fg)%dk
      do i = ghost_sched%recv_faces(fg)%entry_start, &
             ghost_sched%recv_faces(fg)%entry_start + ghost_sched%recv_faces(fg)%entry_count - 1
        Id = ghost_sched%recv_list(i)%dst_i
        Jd = ghost_sched%recv_list(i)%dst_j
        Kd = ghost_sched%recv_list(i)%dst_k
        buf_pos = (i - 1) * gc * nc
        do g = 1, gc
          call Field_Cell_Set(domain%blk(b), h%code, Id-di*g, Jd-dj*g, Kd-dk*g, &
                              h%recv_buf(buf_pos + (g-1)*nc + 1 : buf_pos + g*nc))
        end do
      end do
    end do
  end subroutine halo_unpack_faces


  !> Wait for a halo's sends. Call from !$omp single.
  subroutine halo_wait_send(h)
    use mpi
    implicit none
    type(halo_field_type), intent(inout) :: h
    integer :: ierr

    if (.not. h%init) return
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, h%send_req, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine halo_wait_send


  !> Internal: begin non-blocking exchange of R field using persistent requests.
  subroutine exchange_R_field_begin(domain)
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc
    use mpi

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: ns, nr, i, g, b, fg, ierr, buf_pos
    integer :: Is, Js, Ks, di, dj, dk
    integer :: entry_size

    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv
    if (ns == 0 .and. nr == 0) return

    entry_size = gc

    ! Start persistent receives for R
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_recv_ranks, ghost_sched%R_recv_req_pers, ierr)
      call check_mpi_error(ierr)
    end if

    ! Pack using face groups (reads R instead of T)
    do fg = 1, ghost_sched%n_send_faces
      b  = ghost_sched%send_faces(fg)%block
      di = ghost_sched%send_faces(fg)%di
      dj = ghost_sched%send_faces(fg)%dj
      dk = ghost_sched%send_faces(fg)%dk
      do i = ghost_sched%send_faces(fg)%entry_start, &
             ghost_sched%send_faces(fg)%entry_start + ghost_sched%send_faces(fg)%entry_count - 1
        Is = ghost_sched%send_list(i)%src_i
        Js = ghost_sched%send_list(i)%src_j
        Ks = ghost_sched%send_list(i)%src_k
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          ghost_sched%T_send_buf(buf_pos + g) = &
            domain%blk(b)%R(Is+di*(g-1), Js+dj*(g-1), Ks+dk*(g-1))
        end do
      end do
    end do

    ! Start persistent sends for R
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_send_ranks, ghost_sched%R_send_req_pers, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_R_field_begin


  !> Internal: complete non-blocking exchange and unpack R into ghost cells.
  subroutine exchange_R_field_end(domain)
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: gc
    use mpi

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: nr, i, g, b, fg, ierr, buf_pos
    integer :: Id, Jd, Kd, di, dj, dk
    integer :: entry_size

    nr = ghost_sched%n_recv
    if (ghost_sched%n_send == 0 .and. nr == 0) return

    entry_size = gc

    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%R_recv_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if

    ! Unpack using face groups
    do fg = 1, ghost_sched%n_recv_faces
      b  = ghost_sched%recv_faces(fg)%block
      di = ghost_sched%recv_faces(fg)%di
      dj = ghost_sched%recv_faces(fg)%dj
      dk = ghost_sched%recv_faces(fg)%dk
      do i = ghost_sched%recv_faces(fg)%entry_start, &
             ghost_sched%recv_faces(fg)%entry_start + ghost_sched%recv_faces(fg)%entry_count - 1
        Id = ghost_sched%recv_list(i)%dst_i
        Jd = ghost_sched%recv_list(i)%dst_j
        Kd = ghost_sched%recv_list(i)%dst_k
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          domain%blk(b)%R(Id-di*g, Jd-dj*g, Kd-dk*g) = &
            ghost_sched%T_recv_buf(buf_pos + g)
        end do
      end do
    end do

    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%R_send_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_R_field_end
#endif


  !> MPI barrier for synchronizing ranks after I/O operations.
  !> Prevents non-root ranks from advancing to the next MPI collective
  !> while root is still writing files.
  subroutine mpi_io_barrier()
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
    end block
#endif
  end subroutine mpi_io_barrier

  subroutine Ghost_Interrank( Im, Jm, Km, Fm, blkm, Tg )
    use FUSS_Global_m, only: gc
    use FUSS_Parameters_m, only: guide
    use FUSS_Advanced_Types_m
    implicit none
    integer, intent(in)                  :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout) :: blkm
    real(R8), intent(inout)              :: Tg(6)
    integer                              :: g, Ig, Jg, Kg

    do g = 1, gc
        Ig = Im - guide(Fm,1)*g
        Jg = Jm - guide(Fm,2)*g
        Kg = Km - guide(Fm,3)*g
        Tg(g)  = blkm % T (Ig,Jg,Kg)
    enddo
  end subroutine Ghost_Interrank


end module FUSS_Mod_GhostExchange