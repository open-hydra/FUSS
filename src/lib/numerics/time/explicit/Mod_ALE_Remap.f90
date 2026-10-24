!>@brief Conservative remap of the state onto the moved mesh (ALE).
!>
!> WHY A REMAP RATHER THAN AN EXTRA TERM IN THE RK STAGES
!> -----------------------------------------------------
!> The mesh moves ONCE per time step, but a Runge-Kutta step evaluates its
!> residual once per STAGE. Folding the swept-volume flux and the cell-volume
!> change into the stage residual therefore applies a per-step quantity three
!> times, each with a different RK weight and a different intermediate state.
!>
!> That was measured, not guessed: with the terms inside the RK stages, global
!> energy was off by 5.7e-4 relative and -- the diagnostic that settles it --
!> the error did NOT shrink when dt was halved or quartered at fixed end time
!> (5.678e-4, 5.719e-4, 5.740e-4). A temporal-consistency artifact would have
!> converged; a constant error means the accounting itself is wrong.
!>
!> Splitting the two operations fixes it exactly:
!>
!>   1. REMAP (here). Transfer the state from the old mesh to the new one:
!>
!>        h_new = [ h_old * V_old  +  sum_f h_upwind * dV_swept(f) ] / V_new
!>
!>      applied ONCE per step, outside the RK stages.
!>
!>   2. CONDUCT. The ordinary static-mesh Runge-Kutta step, on the new mesh,
!>      completely unchanged. This is why Lib_Newstate needs no ALE branch.
!>
!> PROPERTIES
!>
!>   Exactly conservative. Summing over cells, sum(V_new h_new) = sum(V_old
!>   h_old) + sum_cells sum_f h_upwind dV. Two cells sharing a face see equal
!>   and opposite dV (verified in test_morph_block to 2.4e-15) and, because the
!>   upwind value is selected on the SIGN of dV, they select the SAME h. The
!>   double sum therefore telescopes to zero over interior faces.
!>
!>   Free-stream preserving. For uniform h the bracket is h*(V_old + (V_new -
!>   V_old)) = h*V_new by the discrete GCL, so h_new = h exactly.
!>
!>   Inert when the mesh does not move: dV == 0 and V_new == V_old give
!>   h_new = h_old, bit-for-bit.
module FUSS_Mod_ALE_Remap
  use FUSS_Mod_MPI, only: fuss_abort_in_thread
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Remap_State_ALE, Remap_Pyrolysis_ALE

  ! ---- the pyrolysis remap (plan 08 section 2.9): fixed slots of the
  ! per-cell remapped vector -- energy E, stored gas G (quasi-steady: the
  ! carry C), the intensive p_pg, eta_prev, T_prev, then the nsub sub-phase
  ! masses. Type 1 uses QE and the masses only.
  integer, parameter :: QE = 1, QG = 2, QP = 3, QEP = 4, QTP = 5, QM = 6
  integer, parameter :: MAXSUB = 64                  ! sub-phases the fixed-size work vectors hold
  real(R8), allocatable, save :: qo(:,:,:,:), qn(:,:,:,:), qwall(:,:)
  logical,  allocatable, save :: moved(:,:,:)
  !> per block, the bc index of the 306 face of a cell on the receding face
  type :: map306_t
    integer, allocatable :: idx(:,:,:)
  end type map306_t
  type(map306_t), allocatable, save :: map306(:)
  integer, allocatable, save :: f306(:)              ! the receding face id of a block, 0 = none

contains

  !> Conservative transfer of the state onto the moved mesh.
  !>
  !> Call ONCE per time step, after the mesh has moved and before Copy_State,
  !> from inside an OpenMP parallel region (it uses `!$omp do`).
  !>
  !> Ghost temperatures must already be valid: the upwind value at a block
  !> interface is read from the ghost cell.
  !>
  !> A receding ablating face (306, dV < 0) sweeps out the wall state, as the
  !> pyrolysis remap does (P-D3; plan 12 U19): the solid at the face's T_w of
  !> its last balance (Tw_prev), the cell's own before any balance.
  subroutine Remap_State_ALE ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_mesh_motion
    use FUSS_Lib_Solid,      only: co_H, co_T
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: i, j, k, b, f, ii, jj, kk, ib
    real(R8) :: h, h_up, hsum, dV, vnew, vold

    if ( .not. obj_mesh_motion%enabled ) return

    !$omp single
    call Map306_Setup ( domain )
    !$omp end single

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle

      ! PASS 1 -- compute the remapped volumetric enthalpy for every cell,
      ! parking it in TO. Reading T and writing TO keeps every cell's upwind
      ! lookup on the SAME pre-remap state; writing T in place here would let a
      ! cell see its neighbour's already-remapped value, the two sides of a face
      ! would disagree on the upwind enthalpy, and conservation would break.
      ! (TO is free: Copy_State overwrites it immediately after this routine.)
      !$omp do collapse(3) schedule(static) private(i,j,k,f,ii,jj,kk,ib,h,h_up,hsum,dV,vnew,vold)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)

        vnew = domain%blk(b)%vol(i,j,k)
        vold = domain%blk(b)%vol_old(i,j,k)

        call co_H( domain%blk(b)%matID(i,j,k), domain%blk(b)%T(i,j,k), h )
        hsum = h * vold

        do f = 1, 6
          dV = domain%blk(b)%dV_swept(f,i,j,k)
          if ( dV == 0.0d0 ) cycle

          if ( dV > 0.0d0 ) then
            ! Face swept outward: the volume gained used to belong to the
            ! neighbour, so it arrives carrying the neighbour's enthalpy.
            ii = i; jj = j; kk = k
            select case (f)
              case (1); ii = i - 1
              case (2); ii = i + 1
              case (3); jj = j - 1
              case (4); jj = j + 1
              case (5); kk = k - 1
              case (6); kk = k + 1
            end select
            call co_H( domain%blk(b)%matID(ii,jj,kk), domain%blk(b)%T(ii,jj,kk), h_up )
          else
            ! Face swept inward: this cell loses volume, carrying its own;
            ! through a receding 306 face, the wall's.
            h_up = h
            if ( f == f306(b) ) then
              ib = map306(b)%idx(i,j,k)
              if ( ib > 0 ) then
                if ( domain%bc(ib)%Tw_prev > 0.0d0 ) &
                  call co_H( domain%blk(b)%matID(i,j,k), domain%bc(ib)%Tw_prev, h_up )
              endif
            endif
          endif

          hsum = hsum + h_up * dV
        enddo

        domain%blk(b)%TO(i,j,k) = hsum / vnew

      enddo; enddo; enddo
      !$omp end do

      ! PASS 2 -- convert back to temperature.
      !$omp do collapse(3) schedule(static) private(i,j,k)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        call co_T( domain%blk(b)%matID(i,j,k), domain%blk(b)%TO(i,j,k), &
                   domain%blk(b)%T(i,j,k) )
      enddo; enddo; enddo
      !$omp end do

    enddo

  end subroutine Remap_State_ALE


  !> The pyrolysis state onto the moved mesh (plan 08 section 2.9, phases 1-2;
  !> phase 0, the T^(n-1) shift, is Darcy_Shift_History before the call).
  !> Inside the parallel region, once per step, before Copy_State; law =
  !> nodeshift only (Check_Pyrolysis_Compatibility guard 6).
  !>
  !> Phase 1 evaluates the remapped vector of every MOVED cell from the
  !> pre-remap state only: q_new = [q_old V_old + sum_f q_up dV_f] / V_new,
  !> q_up = the neighbour's q_old for dV_f > 0 across an interior face, the
  !> cell's own q_old at every block-boundary face (no ghost value enters,
  !> decision E-D11), and on the receding 306 face (dV_f < 0) the WALL state
  !> (decision P-D3, as built): the cell's masses and histories with E taken
  !> at the wall temperature T_w and the env's p_w, and the stored gas
  !> eta_st p_w (quasi-steady: (eta_st - eta_prev) p_w). A cell whose six
  !> dV_swept are exactly 0 and whose volume did not change is skipped and
  !> keeps its state bit for bit. Phase 2 recovers, in this order: chi from
  !> the masses, p_pg, eta_prev, T_prev, eta_st = G/p_pg (quasi-steady
  !> eta_prev + C/p_pg), T last (one bracketed Newton), then the solid caches.
  !>
  !> Between the phases, single-threaded and in the fixed order (k, j, i,
  !> face), the conservation audit (decision E-D9): for E, the total sub-phase
  !> mass and (type 2) G, sum q_new V_new - sum q_old V_old must equal the sum
  !> of q_up dV_f over the block-boundary faces, the values the remap used, to
  !> 1e-14 of sum |q_old| V_old -- else error stop 2.
  subroutine Remap_Pyrolysis_ALE ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_pyrolysis, obj_darcy
    use FUSS_Mod_Pyrolysis,  only: Convert_T_To, Convert_T_From, Convert_Chi_To, Convert_Chi_From
    use FUSS_Lib_Pyrolysis,  only: Energy_Of, Energy_Of2, Solid_State
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b, i, j, k, f, ii, jj, kk, nq, nsub, ib, m
    real(R8) :: acc(QM - 1 + MAXSUB), qup(QM - 1 + MAXSUB), dV, vnew, vold, tau, eps_pg
    logical  :: type2, qs, interior

    nsub  = obj_pyrolysis%nsub
    nq    = QM - 1 + nsub
    type2 = obj_darcy%enabled
    qs    = type2 .and. obj_darcy%quasi_steady

    !$omp single
    call Remap_Setup ( domain, nq )
    !$omp end single

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      associate ( blk => domain%blk(b) )

      ! PHASE 1a -- the pre-remap vector of every cell
      !$omp do collapse(3) schedule(static) private(i,j,k)
      do k = 1, blk%dim(3)
      do j = 1, blk%dim(2)
      do i = 1, blk%dim(1)
        call State_To_Q ( blk, i, j, k, qo(1:nq,i,j,k) )
      enddo ; enddo ; enddo
      !$omp end do

      ! the wall state of every 306 face of this block (P-D3)
      !$omp single
      do ib = 1, domain%nbound
        if ( domain%bc(ib)%type /= 306 .or. domain%bc(ib)%b /= b ) cycle
        call Wall_Q ( blk, domain%bc(ib), qwall(1:nq,ib) )
      enddo
      !$omp end single

      ! PHASE 1b -- the donor-cell sums of the moved cells
      !$omp do collapse(3) schedule(static) private(i,j,k,f,ii,jj,kk,acc,qup,dV,vnew,vold,interior,ib)
      do k = 1, blk%dim(3)
      do j = 1, blk%dim(2)
      do i = 1, blk%dim(1)
        vnew = blk%vol(i,j,k)
        vold = blk%vol_old(i,j,k)
        moved(i,j,k) = any( blk%dV_swept(1:6,i,j,k) /= 0.0d0 ) .or. ( vnew /= vold )
        if ( .not. moved(i,j,k) ) cycle
        acc(1:nq) = qo(1:nq,i,j,k) * vold
        do f = 1, 6
          dV = blk%dV_swept(f,i,j,k)
          if ( dV == 0.0d0 ) cycle
          call Face_Neighbour ( blk%dim, i, j, k, f, ii, jj, kk, interior )
          ib = 0
          if ( dV > 0.0d0 .and. interior ) then
            qup(1:nq) = qo(1:nq,ii,jj,kk)
          else
            qup(1:nq) = qo(1:nq,i,j,k)
            if ( dV < 0.0d0 .and. .not. interior ) then
              if ( f == f306(b) ) ib = map306(b)%idx(i,j,k)
              if ( ib > 0 ) qup(1:nq) = qwall(1:nq,ib)
            endif
          endif
          acc(1:nq) = acc(1:nq) + qup(1:nq) * dV
        enddo
        qn(1:nq,i,j,k) = acc(1:nq) / vnew
      enddo ; enddo ; enddo
      !$omp end do

      ! the conservation audit and the swept-out totals (serial, fixed order)
      !$omp single
      call Remap_Audit ( domain, b, nq )
      !$omp end single

      ! PHASE 2 -- recovery of the moved cells, in the fixed order
      !$omp do collapse(3) schedule(static) private(i,j,k,m,tau,eps_pg)
      do k = 1, blk%dim(3)
      do j = 1, blk%dim(2)
      do i = 1, blk%dim(1)
        if ( .not. moved(i,j,k) ) cycle
        call Convert_Chi_From ( blk, i, j, k, qn(QM:nq,i,j,k) )
        if ( type2 ) then
          blk%p_pg(i,j,k)     = qn(QP,i,j,k)
          blk%eta_prev(i,j,k) = qn(QEP,i,j,k)
          blk%T_prev(i,j,k)   = qn(QTP,i,j,k)
          if ( qs ) then
            blk%eta_st(i,j,k) = blk%eta_prev(i,j,k) + qn(QG,i,j,k) / qn(QP,i,j,k)
          else
            blk%eta_st(i,j,k) = qn(QG,i,j,k) / qn(QP,i,j,k)
          endif
        endif
        call Convert_T_From ( blk, i, j, k, qn(QE:QE,i,j,k) )
        m = nint( blk%matID(i,j,k) )
        call Solid_State ( m, blk%T(i,j,k), blk%chi(:,i,j,k), blk%rho_s(i,j,k), blk%cp_s(i,j,k), &
                           blk%k_s(i,j,k), tau, eps_pg )
      enddo ; enddo ; enddo
      !$omp end do

      end associate
    enddo

  contains

    !> The pre-remap vector of cell (i, j, k): energy, masses and (type 2) the gas state.
    subroutine State_To_Q ( bk, ic, jc, kc, q )
      type(FUSS_block_type), intent(inout) :: bk
      integer,  intent(in)  :: ic, jc, kc
      real(R8), intent(out) :: q(:)
      q = 0.0d0
      call Convert_Chi_To ( bk, ic, jc, kc, q(QM:nq) )
      call Convert_T_To   ( bk, ic, jc, kc, q(QE:QE) )
      if ( type2 ) then
        q(QP)  = bk%p_pg(ic,jc,kc)
        q(QEP) = bk%eta_prev(ic,jc,kc)
        q(QTP) = bk%T_prev(ic,jc,kc)
        if ( qs ) then
          q(QG) = ( bk%eta_st(ic,jc,kc) - bk%eta_prev(ic,jc,kc) ) * bk%p_pg(ic,jc,kc)
        else
          q(QG) = bk%eta_st(ic,jc,kc) * bk%p_pg(ic,jc,kc)
        endif
      endif
    end subroutine State_To_Q

    !> The state a receding 306 face sweeps out of its cell (P-D3): the cell's
    !> vector with E at (T_w, p_w) and the stored gas at p_w.
    subroutine Wall_Q ( bk, bcf, q )
      type(FUSS_block_type), intent(inout) :: bk
      type(FUSS_bc_type),    intent(in)    :: bcf
      real(R8), intent(out) :: q(:)
      real(R8) :: Tw
      integer  :: mm
      q = qo(1:nq,bcf%i,bcf%j,bcf%k)
      mm = nint( bk%matID(bcf%i,bcf%j,bcf%k) )
      Tw = bcf%Tw_prev
      if ( .not. ( Tw > 0.0d0 ) ) Tw = bk%T(bcf%i,bcf%j,bcf%k)
      if ( type2 ) then
        q(QP) = bcf%pw
        if ( qs ) then
          q(QG) = ( bk%eta_st(bcf%i,bcf%j,bcf%k) - bk%eta_prev(bcf%i,bcf%j,bcf%k) ) * bcf%pw
        else
          q(QG) = bk%eta_st(bcf%i,bcf%j,bcf%k) * bcf%pw
        endif
        q(QE) = Energy_Of2( mm, Tw, bk%chi(:,bcf%i,bcf%j,bcf%k), bcf%pw, bk%eta_st(bcf%i,bcf%j,bcf%k) * bcf%pw )
      else
        q(QE) = Energy_Of( mm, Tw, bk%chi(:,bcf%i,bcf%j,bcf%k), obj_pyrolysis%p_gas )
      endif
    end subroutine Wall_Q

  end subroutine Remap_Pyrolysis_ALE


  !> (ii, jj, kk) across face f of cell (i, j, k), and whether that cell is
  !> in the block (an interior face) or beyond it (a block-boundary face).
  pure subroutine Face_Neighbour ( d, i, j, k, f, ii, jj, kk, interior )
    integer, intent(in)  :: d(3), i, j, k, f
    integer, intent(out) :: ii, jj, kk
    logical, intent(out) :: interior
    ii = i ; jj = j ; kk = k
    select case ( f )
      case (1) ; ii = i - 1
      case (2) ; ii = i + 1
      case (3) ; jj = j - 1
      case (4) ; jj = j + 1
      case (5) ; kk = k - 1
      case (6) ; kk = k + 1
    end select
    interior = ii >= 1 .and. ii <= d(1) .and. jj >= 1 .and. jj <= d(2) .and. kk >= 1 .and. kk <= d(3)
  end subroutine Face_Neighbour


  !> Once (single thread): the work arrays and the 306 face map.
  subroutine Remap_Setup ( domain, nq )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_ablation
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    integer, intent(in) :: nq
    integer :: b, n(3)

    if ( allocated(qo) ) return
    if ( nq > QM - 1 + MAXSUB ) then
      write(*,'(A,I0,A,I0)') ' [ERROR] pyrolysis remap: more than ', MAXSUB, ' sub-phases (reactions): ', nq - QM + 1
#ifdef USE_MPI
      call fuss_abort_in_thread ( 1 )   ! with MPI (plan 11); the line below is then not reached
#endif
      error stop 1
    endif
    n = 0
    do b = 1, domain%nb
      n = max( n, domain%blk(b)%dim )
    enddo
    allocate( qo(nq, n(1), n(2), n(3)), qn(nq, n(1), n(2), n(3)), moved(n(1), n(2), n(3)) )
    allocate( qwall(nq, domain%nbound) )
    qwall = 0.0d0
    call Map306_Setup ( domain )
    if ( .not. allocated(obj_ablation%abl) ) then
      allocate( obj_ablation%abl(3, domain%nb) )
      obj_ablation%abl = 0.0d0
    endif
  end subroutine Remap_Setup


  !> Once (single thread): per block, the receding 306 face id (0 = none)
  !> and, for a block that has one, the bc index of the 306 face of each cell
  !> on it (both remaps; idx is read only where f == f306(b) > 0).
  subroutine Map306_Setup ( domain )
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    integer :: b, ib

    if ( allocated(f306) ) return
    allocate( map306(domain%nb), f306(domain%nb) )
    f306 = 0
    do ib = 1, domain%nbound
      if ( domain%bc(ib)%type == 306 ) f306(domain%bc(ib)%b) = domain%bc(ib)%f
    enddo
    do b = 1, domain%nb
      if ( f306(b) == 0 ) cycle
      allocate( map306(b)%idx(domain%blk(b)%dim(1), domain%blk(b)%dim(2), domain%blk(b)%dim(3)) )
      map306(b)%idx = 0
    enddo
    do ib = 1, domain%nbound
      if ( domain%bc(ib)%type /= 306 ) cycle
      b = domain%bc(ib)%b
      map306(b)%idx(domain%bc(ib)%i, domain%bc(ib)%j, domain%bc(ib)%k) = ib
    enddo
  end subroutine Map306_Setup


  !> The conservation audit of one block's remap (E-D9) and the swept-out
  !> totals of its 306 faces. Single thread, fixed order. The audit checks
  !> what the remap carries -- E, the reactive masses F (1 - chi) and G --
  !> exactly; the swept-out solid M_abl adds the non-reacting part rho_c of
  !> the material (rho_s = rho_c + sum F (1 - chi), Solid_State), so that
  !> M_s + M_abl of energy-budget.txt closes. Under quasi-steady G is the
  !> carry, the only stored gas the mass balance holds.
  subroutine Remap_Audit ( domain, b, nq )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_darcy, obj_ablation, obj_pyrolysis
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    integer, intent(in) :: b, nq
    integer  :: i, j, k, f, ii, jj, kk, ib, a
    real(R8) :: lhs(3), rhs(3), scale(3), qa(3), qb(3), qu(3), dV, tol
    logical  :: interior
    character(len=14), parameter :: names(3) = [ 'energy        ', 'reactive solid', 'gas           ' ]

    lhs = 0.0d0 ; rhs = 0.0d0 ; scale = 0.0d0
    associate ( blk => domain%blk(b) )
    do k = 1, blk%dim(3)
    do j = 1, blk%dim(2)
    do i = 1, blk%dim(1)
      qb = [ qo(QE,i,j,k), sum( qo(QM:nq,i,j,k) ), qo(QG,i,j,k) ]
      qa = qb
      if ( moved(i,j,k) ) qa = [ qn(QE,i,j,k), sum( qn(QM:nq,i,j,k) ), qn(QG,i,j,k) ]
      lhs   = lhs + qa * blk%vol(i,j,k) - qb * blk%vol_old(i,j,k)
      scale = scale + abs( qb ) * blk%vol_old(i,j,k)
      if ( .not. moved(i,j,k) ) cycle
      do f = 1, 6
        dV = blk%dV_swept(f,i,j,k)
        if ( dV == 0.0d0 ) cycle
        call Face_Neighbour ( blk%dim, i, j, k, f, ii, jj, kk, interior )
        if ( interior ) cycle
        qu = qb
        ib = 0
        if ( dV < 0.0d0 .and. f == f306(b) ) ib = map306(b)%idx(i,j,k)
        if ( ib > 0 ) then
          qu = [ qwall(QE,ib), sum( qwall(QM:nq,ib) ), qwall(QG,ib) ]
          obj_ablation%abl(1,b) = obj_ablation%abl(1,b) - qu(1) * dV
          obj_ablation%abl(2,b) = obj_ablation%abl(2,b) - &
            ( obj_pyrolysis%mat(nint(blk%matID(i,j,k)))%rho_c + qu(2) ) * dV
          obj_ablation%abl(3,b) = obj_ablation%abl(3,b) - qu(3) * dV
        endif
        rhs = rhs + qu * dV
      enddo
    enddo ; enddo ; enddo
    end associate

    do a = 1, 3
      if ( a == 3 .and. .not. obj_darcy%enabled ) cycle
      tol = 1.0d-14 * scale(a)
      if ( abs( lhs(a) - rhs(a) ) > tol ) then
        write(*,'(A)') ' [ERROR] remap conservation: the '//trim(names(a))//' remapped in block '// &
          trim(adjustl(itoa(b)))//' does not balance its block-boundary sweeps'
        write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') '         sum q_new V_new - sum q_old V_old = ', lhs(a), &
          ', sum of boundary sweeps = ', rhs(a), ', tolerance ', tol
#ifdef USE_MPI
        call fuss_abort_in_thread ( 2 )   ! with MPI (plan 11); the line below is then not reached
#endif
        error stop 2
      endif
    enddo

  contains

    function itoa ( n ) result ( s )
      integer, intent(in) :: n
      character(len=12) :: s
      write(s,'(I0)') n
    end function itoa

  end subroutine Remap_Audit

end module FUSS_Mod_ALE_Remap
