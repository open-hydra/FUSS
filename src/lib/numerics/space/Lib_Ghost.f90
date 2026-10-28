module FUSS_Lib_Ghost
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Parameters_m
  use FUSS_Global_m

  implicit none

contains

  subroutine Fill_Ghost_Cell ( domain )
    use FUSS_Mod_MPI, only: is_local_block, mpi_size_
    use FUSS_Mod_GhostExchange, only: exchange_ghost_T_post_recv, exchange_ghost_T_pack, &
                                      exchange_ghost_T_post_send, exchange_ghost_T_wait_unpack, &
                                      exchange_ghost_T_wait_send, &
                                      exchange_ghost_chimera_post_recv, exchange_ghost_chimera_pack, &
                                      exchange_ghost_chimera_post_send, exchange_ghost_chimera_wait_recv, &
                                      exchange_ghost_chimera_unpack, exchange_ghost_chimera_wait_send, &
                                      Ghost_Interrank, exchange_ghost_Tg, ghost_sched, &
                                      set_active_mg_level
    use FUSS_Mod_Timers, only: timer_comm_begin, timer_comm_end
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: ii, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    integer :: fg
    logical :: mpi_on

    ! On a single rank the exchange has nothing to do, and running it anyway
    ! costs a chain of empty master blocks, each one followed by a barrier.  The
    ! test is uniform over the team, so the worksharing constructs and barriers
    ! it guards are still encountered by every thread.
    mpi_on = mpi_size_ > 1

    if (mpi_on) then
      ! MPI: post persistent receives, pack buffer in parallel, then post sends.
      ! Chimera donor cells travel in their own non-blocking exchange, started
      ! here so it overlaps with the local BC processing below.
      ! MPI is initialised with MPI_THREAD_FUNNELED: every MPI call is made by
      ! the primary thread, in a master construct (named masked since OpenMP
      ! 5.1; master is kept because older compilers lack masked), followed by
      ! an explicit barrier wherever the other threads need its result. A
      ! single construct would let any thread make the call, which FUNNELED
      ! forbids.
      !$omp master
      call set_active_mg_level(domain%mg_level)
      call exchange_ghost_T_post_recv(domain)
      call exchange_ghost_chimera_post_recv(domain)
      !$omp end master
      !$omp barrier

      ! Pack send buffers in parallel: T over face groups, chimera over donor cells
      !$omp do schedule(static) private(fg)
      do fg = 1, ghost_sched%n_send_faces
        call exchange_ghost_T_pack(domain, fg, fg)
      end do
      !$omp end do nowait
      !$omp do schedule(static) private(ii)
      do ii = 1, ghost_sched%n_chim_send
        call exchange_ghost_chimera_pack(domain, ii, ii)
      end do

      ! Post sends (must wait for all packing to complete — implicit barrier from !$omp do)
      !$omp master
      call exchange_ghost_T_post_send(domain)
      call exchange_ghost_chimera_post_send(domain)
      !$omp end master
    end if

    ! Process LOCAL BC entries while MPI communication is in flight
    ! Uses pre-filtered local_bc_idx to avoid scanning all nbound entries
    !$omp do schedule (dynamic) private(ii, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs)
    do ii = 1, domain % n_local_bc
      i  = domain % local_bc_idx(ii)
      Bm = domain % bc(i) % b
      Im = domain % bc(i) % i 
      Jm = domain % bc(i) % j 
      Km = domain % bc(i) % k 
      Fm = domain % bc(i) % f
      select case ( domain % bc(i) % type)
        case(101) ! block connection
          Bs = domain % bc(i) % bs
          if (.not. is_local_block(Bs)) cycle  ! inter-rank handled after MPI completes
          Is = domain % bc(i) % is 
          Js = domain % bc(i) % js 
          Ks = domain % bc(i) % ks 
          Fs = domain % bc(i) % fs
          call Ghost_Connection ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                  Is, Js, Ks, Fs, domain % blk(Bs), &
                                  domain % bc(i) % Tg, domain % bc(i) % matIDg )
        case(301:305)
          call Ghost_Wall ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i), domain % time )
        case(102)
          ! Chimera: processed after the donor-cell MPI exchange completes (below)
          continue
        case default
          call Ghost_Symmetry ( Im, Jm, Km, Fm, domain % blk(Bm) )
      end select
    enddo

    ! MPI: wait for T receives to complete; wait for chimera donor cells and
    ! unpack them into the kept-alive T arrays of remote donor blocks. What is
    ! timed here is the communication the local BC work did not hide, plus the
    ! wait on slower neighbours.
    if (mpi_on) then
      !$omp master
      call timer_comm_begin()
      call exchange_ghost_T_wait_unpack(domain)
      call exchange_ghost_chimera_wait_recv(domain)
      call timer_comm_end()
      !$omp end master
      !$omp barrier

      ! Chimera donor cells are distinct, so the unpack has no write conflicts
      !$omp do schedule(static) private(ii)
      do ii = 1, ghost_sched%n_chim_recv
        call exchange_ghost_chimera_unpack(domain, ii, ii)
      end do
    end if

    ! Chimera ghost fill: all donor data (local and remote) is now current
    !$omp do schedule (dynamic) private(ii, i)
    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      if (domain % bc(i) % type /= 102) cycle
      call Ghost_Chimera ( domain % nb, domain % blk, domain % bc(i) )
    enddo

    ! Chimera Tg(3:6): reads the ghost row written just above, so it must stay
    ! a separate loop (implicit barrier in between). Each entry writes its own Tg
    ! only: beyond the edges of a block face the ghost row is extrapolated into Tg,
    ! not written in the block (Fill_BC_Ghost_Chimera), so the entries of two chimera
    ! faces of a block that meet at an edge write no common location, and one loop
    ! over all the entries has no data race
    !$omp do schedule (dynamic) private(ii, i, Bm, Im, Jm, Km, Fm)
    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      if (domain % bc(i) % type /= 102) cycle
      Bm = domain % bc(i) % b
      Im = domain % bc(i) % i
      Jm = domain % bc(i) % j
      Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
      call Fill_BC_Ghost_Chimera ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Tg )
    enddo

    ! Tg(3:6) of connections whose source block is local. On a face edge they are ghost
    ! cells of Bs, set by the MPI unpack or by Ghost_Chimera (the ghost cells of a chimera
    ! face of Bs; Fill_BC_Ghost_Chimera writes nothing in the block): so this loop comes
    ! after both, and needs no other chimera step. (MOSE fills its connection Pg(:,3:6)
    ! before Ghost_Chimera: a connection face that ends on a chimera face reads that ghost
    ! one stage old there.)
    !$omp do schedule (dynamic) private(ii, i, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do ii = 1, domain % n_local_bs
      i  = domain % local_bs_idx(ii)
      Bs = domain % bc(i) % bs
      Is = domain % bc(i) % is
      Js = domain % bc(i) % js
      Ks = domain % bc(i) % ks
      Fs = domain % bc(i) % fs
      d11s = domain % bc(i) % d11
      d12s = domain % bc(i) % d12
      d21s = domain % bc(i) % d21
      d22s = domain % bc(i) % d22
      call Fill_BC_Ghost_Connection ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, &
                                      domain % blk(Bs), domain % bc(i) % Tg )
    enddo

    ! Process INTER-RANK type-1 entries (Bm local, Bs remote)
    if (mpi_on) then
      !$omp do schedule (dynamic) private(ii, i, Bm, Im, Jm, Km, Fm, Bs)
      do ii = 1, domain % n_local_bc
        i  = domain % local_bc_idx(ii)
        if (domain%bc(i)%type /= 101) cycle
        Bs = domain % bc(i) % bs
        if (is_local_block(Bs)) cycle  ! already processed above
        Bm = domain % bc(i) % b
        Im = domain % bc(i) % i
        Jm = domain % bc(i) % j
        Km = domain % bc(i) % k
        Fm = domain % bc(i) % f
        call Ghost_Interrank(Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Tg)
      enddo

      ! Wait for T sends to complete before reusing buffers
      !$omp master
      call exchange_ghost_T_wait_send(domain)
      call exchange_ghost_chimera_wait_send(domain)
      call exchange_ghost_Tg(domain)
      !$omp end master
      !$omp barrier
    end if

  end subroutine Fill_Ghost_Cell


  subroutine Ghost_Connection ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, Tg, matIDg )

    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm, Is, Js, Ks, Fs
    type(FUSS_block_type), intent(in) :: blks
    type(FUSS_block_type), intent(inout) :: blkm
    real(R8), intent(inout)              :: Tg(6), matIDg
    ! Local
    integer :: g, Ig, Jg, Kg

    do g = 1, gc
      Ig = Im - guide(Fm,1)*(g)
      Jg = Jm - guide(Fm,2)*(g)
      Kg = Km - guide(Fm,3)*(g)
      blkm % T(Ig,Jg,Kg) = blks % T (Is + guide(Fs,1)*(g-1), &
                                     Js + guide(Fs,2)*(g-1), &
                                     Ks + guide(Fs,3)*(g-1)  )
      blkm % matID(Ig,Jg,Kg) = blks % matID (Is + guide(Fs,1)*(g-1), &
                                             Js + guide(Fs,2)*(g-1), &
                                             Ks + guide(Fs,3)*(g-1)  )
      Tg (g) = blkm % T (Ig,Jg,Kg)
      if (g == 1) matIDg = blkm % matID (Ig,Jg,Kg)
    enddo

  end subroutine Ghost_Connection


  subroutine Ghost_Symmetry ( Im, Jm, Km, Fm, blk )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout) :: blk
    ! Local
    integer :: g, gi

    ! The mirrored source stays inside the block (gi clamped): on a 1-cell-wide block the
    ! second layer would otherwise read the opposite face's ghost, written concurrently
    do g = 1, gc
      select case (Fm)
        case(1)
            gi = min ( Im+g-1, blk % dim(1) )
            blk % T(Im-g,Jm,Km) = blk % T (gi,Jm,Km)
            blk % matID(Im-g,Jm,Km) = blk % matID (gi,Jm,Km)
        case(2)
            gi = max ( Im-g+1, 1 )
            blk % T(Im+g,Jm,Km) = blk % T (gi,Jm,Km)
            blk % matID(Im+g,Jm,Km) = blk % matID (gi,Jm,Km)
        case(3)
            gi = min ( Jm+g-1, blk % dim(2) )
            blk % T(Im,Jm-g,Km) = blk % T (Im,gi,Km)
            blk % matID(Im,Jm-g,Km) = blk % matID (Im,gi,Km)
        case(4)
            gi = max ( Jm-g+1, 1 )
            blk % T(Im,Jm+g,Km) = blk % T (Im,gi,Km)
            blk % matID(Im,Jm+g,Km) = blk % matID (Im,gi,Km)
        case(5)
            gi = min ( Km+g-1, blk % dim(3) )
            blk % T(Im,Jm,Km-g) = blk % T (Im,Jm,gi)
            blk % matID(Im,Jm,Km-g) = blk % matID (Im,Jm,gi)
        case(6)
            gi = max ( Km-g+1, 1 )
            blk % T(Im,Jm,Km+g) = blk % T (Im,Jm,gi)
            blk % matID(Im,Jm,Km+g) = blk % matID (Im,Jm,gi)
      end select
    enddo

  end subroutine Ghost_Symmetry


  !> Wall ghosts: the cells mirrored through the wall temperature, 2 Tw - T, exact for a linear field.
  !> Time tables of 301/302 are evaluated here, before the fluxes of the stage.
  subroutine Ghost_Wall ( Im, Jm, Km, Fm, blk, bc, time )
    use FUSS_Lib_BC_Fluxes_Wall, only: Wall_Temperature
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout) :: blk
    type(FUSS_bc_type), intent(inout) :: bc
    real(R8), intent(in) :: time
    ! Local
    integer :: g
    real(R8) :: Tw, q_in

    if (bc % BCtime % exists) then
      if (bc % type == 301) bc % qw = bc % BCtime % update(time)
      if (bc % type == 302) bc % Tw = bc % BCtime % update(time)
    endif
    call Wall_Temperature ( Im, Jm, Km, Fm, blk, bc, .false., Tw, q_in )

    call Ghost_Symmetry ( Im, Jm, Km, Fm, blk )
    do g = 1, gc
      select case (Fm)
        case(1)
          blk % T(Im-g,Jm,Km) = 2d0*Tw - blk % T(Im-g,Jm,Km)
        case(2)
          blk % T(Im+g,Jm,Km) = 2d0*Tw - blk % T(Im+g,Jm,Km)
        case(3)
          blk % T(Im,Jm-g,Km) = 2d0*Tw - blk % T(Im,Jm-g,Km)
        case(4)
          blk % T(Im,Jm+g,Km) = 2d0*Tw - blk % T(Im,Jm+g,Km)
        case(5)
          blk % T(Im,Jm,Km-g) = 2d0*Tw - blk % T(Im,Jm,Km-g)
        case(6)
          blk % T(Im,Jm,Km+g) = 2d0*Tw - blk % T(Im,Jm,Km+g)
      end select
    enddo

  end subroutine Ghost_Wall


  subroutine Ghost_Chimera ( nb, Blk, bc )
    use FUSS_Lib_Solid
    implicit none
    integer, intent(in)                  :: nb
    type(FUSS_block_type), intent(inout) :: Blk(nb)
    type(FUSS_bc_type), intent(inout)    :: bc
    ! Local
    integer  :: Bm, Im, Jm, Km, Fm, Ig, Jg, Kg, Ig2, Jg2, Kg2, Bs, Is, Js, Ks, c 
    real(R8) :: consi, consg

    ! Preliminary assignments
    Bm = bc % b
    Im = bc % i
    Jm = bc % j
    Km = bc % k
    Fm = bc % f

    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)
    Ig2 = Im - guide(Fm,1)*2
    Jg2 = Jm - guide(Fm,2)*2
    Kg2 = Km - guide(Fm,3)*2
          
    ! Ghost T from the volume-weighted enthalpy of the donors. First ghost row:
    consg = 0.d0
    do c = 1, bc % ni(1)
      Bs = bc % donorID(c,1)
      Is = bc % donorID(c,2)
      Js = bc % donorID(c,3)
      Ks = bc % donorID(c,4)
      call co_H ( blk(Bs) % matID (Is,Js,Ks), blk(Bs) % T (Is,Js,Ks), consi )
      consg = consg + consi * bc % volume_fraction(c)
      blk(Bm) % matID (Ig,Jg,Kg) = blk(Bs) % matID (Is,Js,Ks)  ! the donors of a ghost must share one material
    enddo
    call co_T ( blk(Bm) % matID (Ig,Jg,Kg), consg, blk(Bm) % T (Ig,Jg,Kg) )
    bc % Tg (1) = blk(Bm) % T (Ig,Jg,Kg)
    bc % matIDg = blk(Bm) % matID (Ig,Jg,Kg)
          
    ! Second row of ghost cell coordinates
    consg = 0.d0
    do c = bc % ni(1)+1, sum ( bc % ni )
      Bs = bc % donorID(c,1)
      Is = bc % donorID(c,2)
      Js = bc % donorID(c,3)
      Ks = bc % donorID(c,4)
      call co_H ( blk(Bs) % matID (Is,Js,Ks), blk(Bs) % T (Is,Js,Ks), consi )
      consg = consg + consi * bc % volume_fraction(c)
      blk(Bm) % matID (Ig2,Jg2,Kg2) = blk(Bs) % matID (Is,Js,Ks)  ! as for the first row: nothing else sets it
    enddo
    call co_T ( blk(Bm) % matID (Ig2,Jg2,Kg2), consg, blk(Bm) % T (Ig2,Jg2,Kg2) )
    bc % Tg (2) = blk(Bm) % T (Ig2,Jg2,Kg2)
          
  end subroutine Ghost_Chimera



  subroutine Fill_BC_Ghost_Connection ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, blks, Tg )
    implicit none
    integer, intent(in) :: Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(FUSS_block_type), intent(in) :: blks
    real(R8), intent(inout)           :: Tg(6)
    ! Local
    integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4

    select case(Fs)
      case(1:2)
        i1 = Is
        j1 = Js - d11s
        k1 = Ks - d12s
        i2 = Is
        j2 = Js + d11s
        k2 = Ks + d12s
        i3 = Is
        j3 = Js - d21s
        k3 = Ks - d22s
        i4 = Is
        j4 = Js + d21s
        k4 = Ks + d22s
      case(3:4)
        i1 = Is - d11s
        j1 = Js 
        k1 = Ks - d12s
        i2 = Is + d11s
        j2 = Js 
        k2 = Ks + d12s
        i3 = Is - d21s
        j3 = Js
        k3 = Ks - d22s
        i4 = Is + d21s
        j4 = Js 
        k4 = Ks + d22s
      case(5:6)
        i1 = Is - d11s
        j1 = Js - d12s
        k1 = Ks 
        i2 = Is + d11s
        j2 = Js + d12s
        k2 = Ks
        i3 = Is - d21s
        j3 = Js - d22s
        k3 = Ks
        i4 = Is + d21s
        j4 = Js + d22s
        k4 = Ks
    end select

    Tg (3) = blks % T (i1,j1,k1)
    Tg (4) = blks % T (i2,j2,k2)
    Tg (5) = blks % T (i3,j3,k3)
    Tg (6) = blks % T (i4,j4,k4)

  end subroutine Fill_BC_Ghost_Connection


  !> Tg(3:6) of a chimera face: the tangential neighbours of its ghost cell, for the tangential differences of the face flux
  !> (BC_Connection), which takes on the side of the ghost cells the central differences along the face of the row of the
  !> ghost cells, as on the side of the cell those of the row of the cell. Tg(3:4) are the lower and the upper neighbour in
  !> the first tangential direction of the face, Tg(5:6) in the second (j and k for faces 1-2, i and k for 3-4, i and j for
  !> 5-6). Within the block face a neighbour is the ghost cell of the row, filled from the donors (Ghost_Chimera). Beyond an
  !> edge of the block face the row has no ghost cell: the neighbour is extrapolated linearly along the row, 2 a - b from
  !> the ghost cell a and the next one b, so that the central difference there is the one-sided one, exact for a linear
  !> field as that of the row of the cell is (the ghost cell copied there halves it, and the face takes 3/4 of the
  !> derivative of a linear field along that direction). A block one cell wide in a direction has no row to extrapolate
  !> along (the k direction of the 2-D meshes): the neighbour is the ghost cell itself. T is a scalar: no rotation on the
  !> axisymmetric wedge, no property lookup, so the extrapolated value needs no admissibility check. Nothing is written in
  !> the block: two chimera faces that meet at an edge share no written location.
  subroutine Fill_BC_Ghost_Chimera ( Im, Jm, Km, Fm, blk, Tg )
    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(in) :: blk
    real(R8), intent(inout)           :: Tg(6)
    ! Local
    integer :: c(3), g(3), e(3), t(2), n, s, m

    c = [ Im, Jm, Km ]
    g = c - guide(Fm,:)
    select case (Fm)
      case (1:2)
        t = [ 2, 3 ]
      case (3:4)
        t = [ 1, 3 ]
      case default
        t = [ 1, 2 ]
    end select

    do n = 1, 2
      e = 0
      e(t(n)) = 1
      do s = -1, 1, 2
        m = 2*n + 1 + (s+1)/2   ! Tg(3) and Tg(4) for n = 1, Tg(5) and Tg(6) for n = 2
        if ( c(t(n)) + s >= 1 .and. c(t(n)) + s <= blk % dim(t(n)) ) then
          ! Within the block face: the ghost cell of the row
          Tg(m) = blk % T (g(1)+s*e(1),g(2)+s*e(2),g(3)+s*e(3))
        elseif ( blk % dim(t(n)) > 1 ) then
          ! Beyond an edge of the block face: extrapolated along the row
          Tg(m) = 2d0 * blk % T (g(1),g(2),g(3)) - blk % T (g(1)-s*e(1),g(2)-s*e(2),g(3)-s*e(3))
        else
          ! One cell in this direction: the ghost cell
          Tg(m) = blk % T (g(1),g(2),g(3))
        end if
      enddo
    enddo

  end subroutine Fill_BC_Ghost_Chimera


  subroutine Fill_matIDg( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Parameters_m
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: f, lower, upper, i, b
    integer  :: Bm, Im, Jm, Km, Fm, Ig, Jg, Kg
    
    upper = 0
    blocks: do b = 1, domain % nb
      faces: do f = 1, 6
        lower = upper + 1                   ! Update lower bound
        upper = upper + domain % n_bf(b,f)  ! Upper bound: add number of cells on face f of block b
        !$omp do schedule ( dynamic ) private(i, Bm, Im, Jm, Km, Fm)
        do i = lower, upper
          Bm = domain % bc(i) % b
          if (.not. is_local_block(Bm)) cycle
          Im = domain % bc(i) % i 
          Jm = domain % bc(i) % j 
          Km = domain % bc(i) % k 
          Fm = domain % bc(i) % f
          select case ( domain % bc(i) % type )
            case (101,102) ! connection & chimera
              Ig = Im - guide(Fm,1)
              Jg = Jm - guide(Fm,2)
              Kg = Km - guide(Fm,3)
              domain % bc(i) % matIDg = domain % blk(Bm) % matID (Ig,Jg,Kg)
          end select
        enddo
      enddo faces
    enddo blocks

  end subroutine Fill_matIDg
  
end module FUSS_Lib_Ghost