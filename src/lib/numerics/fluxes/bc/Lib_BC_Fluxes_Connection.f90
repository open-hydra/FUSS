module FUSS_Lib_BC_Fluxes_Connection
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Lib_BC_Fluxes, only: Face_Index

  implicit none
  public

contains

  !> Chimera (102): two-sided balance with the cell metric and the ghost metric, residual update of the cell.
  subroutine BC_Connection ( Im, Jm, Km, Fm, Blk, Mg, Tg, matIDg )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_tensor_3D_type), intent(in) :: Mg
    real(R8), intent(in)         :: Tg(6), matIDg
    type(FUSS_block_type), intent(inout) :: Blk
    ! Local
    real(R8) :: Flux

    call Connection_Face_Flux_Tint ( Im, Jm, Km, Fm, Blk, Mg, Tg, matIDg, Flux )
    Blk % r(Im,Jm,Km) = Blk % r(Im,Jm,Km) - ( 1 - 2 * mod(Fm,2) ) * Flux

  end subroutine BC_Connection


  !> Flux of a face between blocks or with a chimera boundary from the two-sided balance of the temperature at the
  !> face (Tint) with the cell metric and the ghost metric Mg: the form of 102 for every metric, and of 101 with
  !> diffusive-metric = mean (the inverse-mean form of the former metric). Flux along the normal of growing index.
  subroutine Connection_Face_Flux_Tint ( Im, Jm, Km, Fm, Blk, Mg, Tg, matIDg, Flux )
    use FUSS_Lib_Solid
    use FUSS_Lib_Diffusive
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_tensor_3D_type), intent(in) :: Mg
    real(R8), intent(in)         :: Tg(6), matIDg
    type(FUSS_block_type), intent(in) :: Blk
    real(R8), intent(out)        :: Flux
    ! Local
    integer :: Dir, Face_im, Face_jm, Face_km, modfm2
    real(R8) :: T_loc, T_con, T_ip_loc, T_im_loc, T_jp_loc, T_jm_loc, T_kp_loc, T_km_loc
    real(R8) :: Normal(3), Area, M_loc(3,3), M_con(3,3), Gradient(3)
    real(R8) :: kappa_loc, kappa_con, Num_loc, Num_con, Den_loc, Den_con
    real(R8) :: Tint
    
    ! Boundary face index
    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_im, Face_jm, Face_km )
    modfm2 = 1 - 2 * mod (Fm,2)

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_im,Face_jm,Face_km) % N
    Area   = Blk % dir(Dir) % f(Face_im,Face_jm,Face_km) % A
    M_loc  = Blk % M(Im,Jm,Km) % c
    M_con  = Mg % c

    ! Stencil building
    T_loc    = Blk % T (Im,Jm,Km)
    T_ip_loc = Blk % T (Im+1,Jm,Km)
    T_im_loc = Blk % T (Im-1,Jm,Km)
    T_jp_loc = Blk % T (Im,Jm+1,Km)
    T_jm_loc = Blk % T (Im,Jm-1,Km)
    if (ndir == 2) then
      T_kp_loc = huge(0.0)
      T_km_loc = huge(0.0)
    else
      T_kp_loc = Blk % T (Im,Jm,Km+1)
      T_km_loc = Blk % T (Im,Jm,Km-1)
    endif
    
    T_con = Tg(1)

    call co_Kappa( Blk % matID(Im,Jm,Km), T_loc, kappa_loc )
    call co_Kappa( matIDg, T_con, kappa_con )
          
    ! Tint computations
    select case ( Fm )
      case ( 1 : 2 )
        Den_loc = 2.0d0 * dot_product( M_loc(1,:), Normal ) * modfm2
        Num_loc = 2.0d0*T_loc * dot_product( M_loc(1,:), Normal ) * modfm2 &
                - 0.5d0*(T_jp_loc-T_jm_loc) * dot_product( M_loc(2,:), Normal ) &
                - 0.5d0*(T_kp_loc-T_km_loc) * dot_product( M_loc(3,:), Normal )
        Den_con = 2.0d0 * dot_product( M_con(1,:), Normal ) * modfm2
        Num_con = 2.0d0*T_con * dot_product( M_con(1,:), Normal ) * modfm2 &
                + 0.5d0*(Tg(4)-Tg(3)) * dot_product( M_con(2,:), Normal ) &
                + 0.5d0*(Tg(6)-Tg(5)) * dot_product( M_con(3,:), Normal )
      case ( 3 : 4 )
        Den_loc = 2.0d0 * dot_product( M_loc(2,:), Normal ) * modfm2
        Num_loc = - 0.5d0*(T_ip_loc-T_im_loc) * dot_product( M_loc(1,:), Normal ) &
                + 2.0d0*T_loc * dot_product( M_loc(2,:), Normal ) * modfm2 &
                - 0.5d0*(T_kp_loc-T_km_loc) * dot_product( M_loc(3,:), Normal )
        Den_con = 2.0d0 * dot_product( M_con(2,:), Normal ) * modfm2
        Num_con = 0.5d0*(Tg(4)-Tg(3)) * dot_product( M_con(1,:), Normal ) &
                + 2.0d0*T_con * dot_product( M_con(2,:), Normal ) * modfm2 &
                + 0.5d0*(Tg(6)-Tg(5)) * dot_product( M_con(3,:), Normal )
      case ( 5 : 6 )
        Den_loc = 2.0d0 * dot_product( M_loc(3,:), Normal ) * modfm2
        Num_loc = - 0.5d0*(T_ip_loc-T_im_loc) * dot_product( M_loc(1,:), Normal ) &
                - 0.5d0*(T_jp_loc-T_jm_loc) * dot_product( M_loc(2,:), Normal ) &
                + 2.0d0*T_loc * dot_product( M_loc(3,:), Normal ) * modfm2
        Den_con = 2.0d0 * dot_product( M_con(3,:), Normal ) * modfm2
        Num_con = 0.5d0*(Tg(4)-Tg(3)) * dot_product( M_con(1,:), Normal ) &
                + 0.5d0*(Tg(6)-Tg(5)) * dot_product( M_con(2,:), Normal ) &
                + 2.0d0*T_con * dot_product( M_con(3,:), Normal ) * modfm2
    end select

    Tint = (kappa_loc*Num_loc + kappa_con*Num_con) / (kappa_loc*Den_loc + kappa_con*Den_con)
    
    ! Gradient computations
    select case ( Fm )
    case ( 1 : 2 )
      Gradient (1) = 2.0d0 * ( Tint - T_loc ) * modfm2
      Gradient (2) = ( T_jp_loc - T_jm_loc ) * 0.5d0
      Gradient (3) = ( T_kp_loc - T_km_loc ) * 0.5d0
    case ( 3 : 4 )
      Gradient (1) = ( T_ip_loc - T_im_loc ) * 0.5d0
      Gradient (2) = 2.0d0 * ( Tint - T_loc ) * modfm2
      Gradient (3) = ( T_kp_loc - T_km_loc ) * 0.5d0
    case ( 5 : 6 )
      Gradient (1) = ( T_ip_loc - T_im_loc ) * 0.5d0
      Gradient (2) = ( T_jp_loc - T_jm_loc ) * 0.5d0
      Gradient (3) = 2.0d0 * ( Tint - T_loc ) * modfm2
    end select

    if (ndir == 2) Gradient(3) = 0.0
    Gradient = matmul ( Gradient, M_loc )

    call Compute_Diffusive_Flux ( kappa_loc, Gradient, Area, Normal, Flux )

  end subroutine Connection_Face_Flux_Tint


  !> Block connection (101) with the stencil face metric: the flux of the face as an interior face, from one gradient of
  !> the two cells in the order of growing index of the block (L, R), the tangential neighbours of the other block
  !> in Tg(3:6), the metric of Lib_Stencil_Metrics (bc%Mf) and the conductivity of the interior faces
  !> (Diffusive_Flux): k at the mean temperature for one material, the two half cells in series between two. The
  !> chimera faces (102) keep BC_Connection. Flux: heat through the face along the normal of growing index.
  subroutine Connection_Face_Flux ( Im, Jm, Km, Fm, Blk, bc, Flux )
    use FUSS_Lib_Solid
    use FUSS_Lib_Diffusive
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(in) :: Blk
    type(FUSS_bc_type), intent(in)    :: bc
    real(R8), intent(out) :: Flux
    ! Local
    integer  :: Dir, Face_im, Face_jm, Face_km, aa, bb, cc, t, q, tb(2), e(3,3), c(3)
    real(R8) :: T_loc, T_con, T1, T2, Normal(3), Area, Gradient(3), kappa, dloc, dgh
    real(R8) :: Mi_loc(3,3), Mi_con(3,3), h_loc, h_con, mat1, mat2

    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_im, Face_jm, Face_km )
    Normal = Blk % dir(Dir) % f(Face_im,Face_jm,Face_km) % N
    Area   = Blk % dir(Dir) % f(Face_im,Face_jm,Face_km) % A
    select case (Dir)
      case (1) ; aa = 1 ; tb = [2, 3]
      case (2) ; aa = 2 ; tb = [1, 3]
      case default ; aa = 3 ; tb = [1, 2]
    end select
    e = reshape ( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    c = [Im, Jm, Km]

    T_loc = Blk % T (Im,Jm,Km)
    T_con = bc % Tg(1)

    ! One gradient: across the face, then the mean of the central differences of the two cells along each tangential
    ! direction (the neighbours of the other block are the cells of its side, Tg(3:6))
    Gradient = 0d0
    if ( mod(Fm,2) == 0 ) then
      T1 = T_loc ; T2 = T_con
    else
      T1 = T_con ; T2 = T_loc
    end if
    Gradient(aa) = T2 - T1
    do q = 1, 2
      t = tb(q)
      if ( t == 3 .and. ndir == 2 ) cycle
      dloc = Blk % T (Im+e(1,t),Jm+e(2,t),Km+e(3,t)) - Blk % T (Im-e(1,t),Jm-e(2,t),Km-e(3,t))
      dgh  = bc % Tg(2*q+2) - bc % Tg(2*q+1)
      if ( mod(Fm,2) == 0 ) then
        Gradient(t) = 0.25d0 * ( dloc + dgh )
      else
        Gradient(t) = 0.25d0 * ( dgh + dloc )
      end if
    enddo
    Gradient = matmul ( Gradient, bc % Mf )

    ! Conductivity as at an interior face (lower cell first)
    if ( Blk % matID(Im,Jm,Km) == bc % matIDg ) then
      call co_Kappa( Blk % matID(Im,Jm,Km), 0.5d0 * ( T1 + T2 ), kappa )
    else
      Mi_loc = Inverse3 ( Blk % M(Im,Jm,Km) % c )
      Mi_con = Inverse3 ( bc % Mg(1) % c )
      h_loc = dot_product ( Mi_loc(:,aa), Normal )
      h_con = dot_product ( Mi_con(:,aa), Normal )
      if ( mod(Fm,2) == 0 ) then
        kappa = Series_Kappa ( Blk % matID(Im,Jm,Km), T_loc, h_loc, bc % matIDg, T_con, h_con )
      else
        kappa = Series_Kappa ( bc % matIDg, T_con, h_con, Blk % matID(Im,Jm,Km), T_loc, h_loc )
      end if
    end if

    call Compute_Diffusive_Flux ( kappa, Gradient, Area, Normal, Flux )

  end subroutine Connection_Face_Flux


  !> Pairs the two records of each connection face (101): the record of the same face seen from the other block
  !> (partner), and which of the two computes the flux of the face for both cells (owner, Connection_Fluxes). The
  !> owner is the record whose face is the high-index one of its cell (2, 4, 6) when the other is a low-index one, as
  !> on an interior face, where the flux between cells n and n+1 is computed from cell n. When the two faces have the
  !> same parity (blocks that meet with their index directions permuted) the owner is the record with the smaller
  !> (block, face, k, j, i), a choice that does not depend on the order of the records in the BC file. A record whose
  !> partner is missing, does not point back at it or points at itself stops the run. Called on every rank (the records
  !> of all blocks exist on all ranks). Port of MOSE ec1e480 (Pair_Connections of IO_BC).
  subroutine Pair_Connections ( domain )
    use FUSS_Mod_MPI, only: mpi_abort_all
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    type :: face_records_type
      integer, allocatable :: r(:,:)
    end type face_records_type
    type(face_records_type), allocatable :: rec(:,:)
    integer :: b, f, i, p, nt(2), t(2)
    character(len=160) :: msg

    allocate ( rec(domain % nb, 6) )

    ! Map of the connection records of each block face, by the two tangential indices of the cell
    do i = 1, domain % nbound
      if ( domain % bc(i) % type /= 101 ) cycle
      b = domain % bc(i) % b
      f = domain % bc(i) % f
      if ( .not. allocated(rec(b,f) % r) ) then
        nt = Tangential ( f, domain % blk(b) % dim(1), domain % blk(b) % dim(2), domain % blk(b) % dim(3) )
        allocate ( rec(b,f) % r(nt(1), nt(2)) )
        rec(b,f) % r = 0
      end if
      t = Tangential ( f, domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k )
      rec(b,f) % r(t(1), t(2)) = i
    end do

    ! Partner and owner of each connection record
    do i = 1, domain % nbound
      if ( domain % bc(i) % type /= 101 ) cycle
      b = domain % bc(i) % bs
      f = domain % bc(i) % fs
      p = 0
      if ( b >= 1 .and. b <= domain % nb .and. f >= 1 .and. f <= 6 ) then
        if ( allocated(rec(b,f) % r) ) then
          t = Tangential ( f, domain % bc(i) % is, domain % bc(i) % js, domain % bc(i) % ks )
          if ( all(t >= 1) .and. all(t <= shape(rec(b,f) % r)) ) p = rec(b,f) % r(t(1), t(2))
        end if
      end if
      if ( p == 0 ) then
        write(msg,'(A,I0,A)') 'Connection record ', i, ': no connection record at the cell and face it points at'
        call mpi_abort_all ( trim(msg) )
      end if
      if ( p == i ) then
        write(msg,'(A,I0,A)') 'Connection record ', i, ' points at its own cell and face'
        call mpi_abort_all ( trim(msg) )
      end if
      if ( domain % bc(p) % type /= 101 .or. domain % bc(p) % bs /= domain % bc(i) % b .or. &
           domain % bc(p) % fs /= domain % bc(i) % f .or. domain % bc(p) % is /= domain % bc(i) % i .or. &
           domain % bc(p) % js /= domain % bc(i) % j .or. domain % bc(p) % ks /= domain % bc(i) % k ) then
        write(msg,'(A,I0,A,I0,A)') 'Connection records ', i, ' and ', p, ' do not point at each other'
        call mpi_abort_all ( trim(msg) )
      end if
      domain % bc(i) % partner = p
      if ( mod(domain % bc(i) % f, 2) == 0 .and. mod(domain % bc(p) % f, 2) == 1 ) then
        domain % bc(i) % owner = .true.
      else if ( mod(domain % bc(i) % f, 2) == 1 .and. mod(domain % bc(p) % f, 2) == 0 ) then
        domain % bc(i) % owner = .false.
      else
        domain % bc(i) % owner = Key_Less ( domain % bc(i) % b, domain % bc(i) % f, domain % bc(i) % k, &
                                            domain % bc(i) % j, domain % bc(i) % i, &
                                            domain % bc(p) % b, domain % bc(p) % f, domain % bc(p) % k, &
                                            domain % bc(p) % j, domain % bc(p) % i )
      end if
    end do

    contains

      !> Whether the record (b1, f1, k1, j1, i1) comes before (b2, f2, k2, j2, i2): block, face, then k, j, i
      pure logical function Key_Less ( b1, f1, k1, j1, i1, b2, f2, k2, j2, i2 )
        implicit none
        integer, intent(in) :: b1, f1, k1, j1, i1, b2, f2, k2, j2, i2
        integer :: key1(5), key2(5), n

        key1 = [ b1, f1, k1, j1, i1 ]
        key2 = [ b2, f2, k2, j2, i2 ]
        Key_Less = .false.
        do n = 1, 5
          if ( key1(n) /= key2(n) ) then
            Key_Less = key1(n) < key2(n)
            return
          end if
        end do

      end function Key_Less

      !> The two indices of a cell (or the two dimensions of a block) along the face f
      pure function Tangential ( f, i, j, k ) result ( t )
        implicit none
        integer, intent(in) :: f, i, j, k
        integer :: t(2)

        select case ( f )
          case (1:2)
            t = [ j, k ]
          case (3:4)
            t = [ i, k ]
          case default
            t = [ i, j ]
        end select

      end function Tangential

  end subroutine Pair_Connections


  !> The owner of a connection face computes its flux once and keeps the share of its cell in bc%Fc (heat into the
  !> cell, to be added to the residual by BC_Fluxes) and gives the opposite share to the record of the other cell.
  !> Called by Fluxes, from inside the parallel region: an orphaned loop over the owner records of the rank. The
  !> partner on another rank gets its share by exchange_connection_flux.
  subroutine Connection_Fluxes ( domain )
    use FUSS_Mod_MPI, only: mpi_size_
    use FUSS_Mod_GhostExchange, only: exchange_connection_flux
    use FUSS_Config_Types_m, only: obj_space_scheme
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: ii, i, p, Bm
    real(R8) :: Flux, s

    !$omp do schedule (dynamic, 64) private(ii, i, p, Bm, Flux, s)
    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      if ( domain % bc(i) % type /= 101 ) cycle
      if ( .not. domain % bc(i) % owner ) cycle
      Bm = domain % bc(i) % b
      p  = domain % bc(i) % partner
      if ( obj_space_scheme % stencil_metric ) then
        call Connection_Face_Flux ( domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k, domain % bc(i) % f, &
                                    domain % blk(Bm), domain % bc(i), Flux )
      else
        call Connection_Face_Flux_Tint ( domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k, domain % bc(i) % f, &
                                         domain % blk(Bm), domain % bc(i) % Mg(1), domain % bc(i) % Tg, &
                                         domain % bc(i) % matIDg, Flux )
      end if
      s = 1 - 2 * mod(domain % bc(i) % f, 2)     ! +1: the cell is on the low side of the face
      domain % bc(i) % Fc = - s * Flux
      domain % bc(p) % Fc = + s * Flux
    enddo

    ! Faces whose other cell is on another rank (MPI calls from the primary thread)
    if ( mpi_size_ > 1 ) then
      !$omp master
      call exchange_connection_flux ( domain )
      !$omp end master
      !$omp barrier
    end if

  end subroutine Connection_Fluxes


end module FUSS_Lib_BC_Fluxes_Connection