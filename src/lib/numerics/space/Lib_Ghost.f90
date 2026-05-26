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
                                      Ghost_Interrank, exchange_ghost_Tg, ghost_sched
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: ii, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    integer :: fg


    ! MPI: post persistent receives, pack buffer in parallel, then post sends
    !$omp single
    call exchange_ghost_T_post_recv(domain)
    !$omp end single

    ! Pack send buffer in parallel over face groups
    !$omp do schedule(static) private(fg)
    do fg = 1, ghost_sched%n_send_faces
      call exchange_ghost_T_pack(domain, fg, fg)
    end do

    ! Post sends (must wait for all packing to complete — implicit barrier from !$omp do)
    !$omp single
    call exchange_ghost_T_post_send(domain)
    !$omp end single nowait

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
        case(300)
          call Ghost_Symmetry ( Im, Jm, Km, Fm, domain % blk(Bm) )
        case(102)
          call Ghost_Chimera ( domain % nb, domain % blk, domain % bc(i) )
        case default
          call Ghost_Symmetry ( Im, Jm, Km, Fm, domain % blk(Bm) )
      end select
    enddo

    ! Compute Tg(3:6) for type-1 connections where source block is local.
    ! This only reads local blk%T data, so it can run before the MPI exchange.
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

    ! Compute Tg for chimera (666) connections where Bm is local.
    !$omp do schedule (dynamic) private(ii, i, Bm, Im, Jm, Km, Fm)
    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      Bm = domain % bc(i) % b
      Im = domain % bc(i) % i
      Jm = domain % bc(i) % j
      Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
      select case (domain % bc(i) % type)
        case(102)
          call Fill_BC_Ghost_Chimera ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Tg )
      end select
    enddo

    ! MPI: wait for P receives to complete
    !$omp single
    call exchange_ghost_T_wait_unpack(domain)
    !$omp end single

    ! Process INTER-RANK type-1 entries (Bm local, Bs remote)
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

    ! Wait for P sends to complete before reusing buffers
    !$omp single
    call exchange_ghost_T_wait_send(domain)
    call exchange_ghost_Tg(domain)
    !$omp end single

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
    integer :: g

    do g = 1, gc
      select case (Fm)
        case(1)
            blk % T(Im-g,Jm,Km) = blk % T (Im+g-1,Jm,Km)
            blk % matID(Im-g,Jm,Km) = blk % matID (Im+g-1,Jm,Km)
        case(2)
            blk % T(Im+g,Jm,Km) = blk % T (Im-g+1,Jm,Km)
            blk % matID(Im+g,Jm,Km) = blk % matID (Im-g+1,Jm,Km)
        case(3)
            blk % T(Im,Jm-g,Km) = blk % T (Im,Jm+g-1,Km)
            blk % matID(Im,Jm-g,Km) = blk % matID (Im,Jm+g-1,Km)
        case(4)
            blk % T(Im,Jm+g,Km) = blk % T (Im,Jm-g+1,Km)
            blk % matID(Im,Jm+g,Km) = blk % matID (Im,Jm-g+1,Km)
        case(5)
            blk % T(Im,Jm,Km-g) = blk % T (Im,Jm,Km+g-1)
            blk % matID(Im,Jm,Km-g) = blk % matID (Im,Jm,Km+g-1)
        case(6)
            blk % T(Im,Jm,Km+g) = blk % T (Im,Jm,Km-g+1)
            blk % matID(Im,Jm,Km+g) = blk % matID (Im,Jm,Km-g+1)
      end select
    enddo

  end subroutine Ghost_Symmetry


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
          
    ! RIEMPIMPENTO GHOST CON VARIABILI CONSERVATE
    ! First row of ghost cell coordinates
    consg = 0.d0
    do c = 1, bc % ni(1)
      Bs = bc % donorID(c,1)
      Is = bc % donorID(c,2)
      Js = bc % donorID(c,3)
      Ks = bc % donorID(c,4)
      call co_H ( blk(Bs) % matID (Is,Js,Ks), blk(Bs) % T (Is,Js,Ks), consi )
      consg = consg + consi * bc % volume_fraction(c)
      blk(Bm) % matID (Ig,Jg,Kg) = blk(Bs) % matID (Is,Js,Ks)  ! WARNING: non possono essere connessi contemporaneamente due blocchi di materiali diversi!
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
    enddo
    call co_T ( blk(Bm) % matID (Ig2,Jg2,Kg2), consg, blk(Bm) % T (Ig2,Jg2,Kg2) )
    bc % Tg (2) = blk(Bm) % T (Ig2,Jg2,Kg2)
          
  end subroutine Ghost_Chimera


  subroutine Ghost_Extrapolate ( Im, Jm, Km, Fm, blk )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout) :: blk
    ! Local
    integer :: g

    ! Extrapolation with 2nd order accuracy. Theory: suppose we have a stencil x1,x2,x3,x4 with x4 unknown.
    ! Taylor 2nd order: x4 = x3 + f1(x3) + 1/2 f2(x3) with f1 and f2 1st and 2nd derivative. Approximate 
    ! f1(x3) = (x4 - x2)/2 central and f2(x3) = x3 - 2x2 + x1 left sided => x4 = 3*x3 - 3*x2 + x1
    do g = 1, gc
      select case (Fm)
        case(1)
          blk % T(Im-g,Jm,Km) = 3d0 * blk % T (Im-g+1,Jm,Km) - &
                                3d0 * blk % T (Im-g+2,Jm,Km) + &
                                      blk % T (Im-g+3,Jm,Km)
          blk % matID(Im-g,Jm,Km) = blk % matID (Im,Jm,Km)
        case(2)
          blk % T(Im+g,Jm,Km) = 3d0 * blk % T (Im+g-1,Jm,Km) - &
                                3d0 * blk % T (Im+g-2,Jm,Km) + &
                                      blk % T (Im+g-3,Jm,Km)
          blk % matID(Im+g,Jm,Km) = blk % matID (Im,Jm,Km)
        case(3)
          blk % T(Im,Jm-g,Km) = 3d0 * blk % T (Im,Jm-g+1,Km) - &
                                3d0 * blk % T (Im,Jm-g+2,Km) + &
                                      blk % T (Im,Jm-g+3,Km)
          blk % matID(Im,Jm-g,Km) = blk % matID (Im,Jm,Km)
        case(4)
          blk % T(Im,Jm+g,Km) = 3d0 * blk % T (Im,Jm+g-1,Km) - &
                                3d0 * blk % T (Im,Jm+g-2,Km) + &
                                      blk % T (Im,Jm+g-3,Km)
          blk % matID(Im,Jm+g,Km) = blk % matID (Im,Jm,Km)
        case(5)
          blk % T(Im,Jm,Km-g) = 3d0 * blk % T (Im,Jm,Km-g+1) - &
                                3d0 * blk % T (Im,Jm,Km-g+2) + &
                                      blk % T (Im,Jm,Km-g+3)
          blk % matID(Im,Jm,Km-g) = blk % matID (Im,Jm,Km)
        case(6)
          blk % T(Im,Jm,Km+g) = 3d0 * blk % T (Im,Jm,Km+g-1) - &
                                3d0 * blk % T (Im,Jm,Km+g-2) + &
                                      blk % T (Im,Jm,Km+g-3)
          blk % matID(Im,Jm,Km+g) = blk % matID (Im,Jm,Km)
      end select
    enddo
      
  end subroutine Ghost_Extrapolate


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


  subroutine Fill_BC_Ghost_Chimera ( Im, Jm, Km, Fm, blk, Tg )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout) :: blk
    real(R8), intent(inout)              :: Tg(6)
    ! Local
    integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4
    integer :: dim(3), Ig, Jg, Kg

    dim = blk % dim
    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)

    select case (Fm)
      case(1:2)
        if ( Jm == 1 )      blk % T (Ig, Jm-1, Km) = blk % T (Ig, Jm, Km)
        if ( Jm == dim(2) ) blk % T (Ig, Jm+1, Km) = blk % T (Ig, Jm, Km)
        if ( Km == 1 )      blk % T (Ig, Jm, Km-1) = blk % T (Ig, Jm, Km)
        if ( Km == dim(3) ) blk % T (Ig, Jm, Km+1) = blk % T (Ig, Jm, Km)
      case(3:4)
        if ( Im == 1 )      blk % T (Im-1, Jg, Km) = blk % T (Im, Jg, Km)
        if ( Im == dim(1) ) blk % T (Im+1, Jg, Km) = blk % T (Im, Jg, Km)
        if ( Km == 1 )      blk % T (Im, Jg, Km-1) = blk % T (Im, Jg, Km)
        if ( Km == dim(3) ) blk % T (Im, Jg, Km+1) = blk % T (Im, Jg, Km)
      case(5:6)
        if ( Im == 1 )      blk % T (Im-1, Jm, Kg) = blk % T (Im, Jm, Kg)
        if ( Im == dim(1) ) blk % T (Im+1, Jm, Kg) = blk % T (Im, Jm, Kg)
        if ( Jm == 1 )      blk % T (Im, Jm-1, Kg) = blk % T (Im, Jm, Kg)
        if ( Jm == dim(2) ) blk % T (Im, Jm+1, Kg) = blk % T (Im, Jm, Kg)
    end select

    select case(Fm)
      case(1:2)
        i1 = Im
        j1 = Jm - 1
        k1 = Km
        i2 = Im
        j2 = Jm + 1
        k2 = Km
        i3 = Im
        j3 = Jm
        k3 = Km - 1
        i4 = Im
        j4 = Jm
        k4 = Km + 1
      case(3:4)
        i1 = Im - 1
        j1 = Jm
        k1 = Km
        i2 = Im + 1
        j2 = Jm 
        k2 = Km
        i3 = Im
        j3 = Jm
        k3 = Km - 1
        i4 = Im
        j4 = Jm 
        k4 = Km + 1
      case(5:6)
        i1 = Im - 1
        j1 = Jm
        k1 = Km
        i2 = Im + 1
        j2 = Jm
        k2 = Km
        i3 = Im
        j3 = Jm - 1
        k3 = Km
        i4 = Im
        j4 = Jm + 1
        k4 = Km
    end select

    Tg (3) = Blk % T (i1,j1,k1)
    Tg (4) = Blk % T (i2,j2,k2)
    Tg (5) = Blk % T (i3,j3,k3)
    Tg (6) = Blk % T (i4,j4,k4)

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