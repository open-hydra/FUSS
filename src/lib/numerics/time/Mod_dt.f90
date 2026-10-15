module FUSS_Mod_dt
  use iso_fortran_env, only: I1 => int8, I4 => int32, R8 => real64

  implicit none
  private
  public :: Compute_dt, Set_Global_dt

  ! dt-face-k: the faces of each local cell that are block connections (101) or chimera faces
  ! (102), bit f-1 for face f, per multigrid level; built at the first call on the level
  type :: conn_blk_t
    integer(I1), allocatable :: f(:,:,:)
  end type conn_blk_t
  type :: conn_level_t
    type(conn_blk_t), allocatable :: blk(:)
  end type conn_level_t
  type(conn_level_t), allocatable :: conn(:)

contains

  subroutine Compute_dt ( domain, vnn, rampa_iter )
    use FUSS_Advanced_Types_m
    use FUSS_Global_m
    use FUSS_Config_Types_m, only: obj_time_scheme
    use FUSS_Mod_MPI, only: is_local_block, mpi_allreduce_min_r8
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: vnn
    integer, intent(in)  :: rampa_iter
    ! Local
    integer  :: i, j, k, b
    real(R8) :: dtcell, dtglobal, dtglobal_mpi
    
    dtglobal = domain % dtglobal

    do b = 1, domain % nb ! loop over blocks
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(3) private ( dtcell ), reduction ( min : dtglobal )
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        
        ! Compute local cell dt according to VNN number
        call compute (T = domain % blk(b) % T(i,j,k), &
                      matID = domain % blk(b) % matID(i,j,k), &
                      M = domain % blk(b) % M(i,j,k) % c, &
                      dtmin = dtcell, &
                      vnn = vnn )

        ! Apply VNN reduction if required
        if ( domain%iter < rampa_iter ) dtcell = dtcell * domain%iter / rampa_iter

        ! Update local cell dt and global minimum dt
        domain % blk(b) % dtlocal(i,j,k) = dtcell
        dtglobal = min ( dtcell, dtglobal )

      enddo; enddo; enddo
      !$omp end parallel
    enddo ! end of loop over blocks

    ! dt-face-k: the cells whose face conductivities differ from their own get the face step
    if ( obj_time_scheme%dt_face_k ) call Face_k_dt ( domain, vnn, rampa_iter, dtglobal )

    ! MPI: global minimum across all ranks
    call mpi_allreduce_min_r8(dtglobal, dtglobal_mpi)
    domain % dtglobal = dtglobal_mpi

    contains
      
      ! Explicit diffusion limit alpha dt sum_d |grad xi_d|^2 <= vnn/2, grad xi_d = row d of M:
      ! exact on skewed cells too, where the edge lengths overstate the stable step by 1/sin^2(phi)
      subroutine compute ( T, matID, M, dtmin, vnn )
        use FUSS_Global_m
        use FUSS_Lib_Solid
        implicit none
        real(R8), intent(in)  :: T, matID, M(3,3), vnn
        real(R8), intent(out) :: dtmin
        ! Local
        integer :: d
        real(R8) :: kappa, rho, cp, alpha, G, dt

        ! Thermal diffusivity
        call co_Kappa( matID, T, kappa )
        call co_Rho( matID, T, rho )
        call co_Cp( matID, T, cp )
        alpha = kappa/rho/cp

        dtmin = 1d8

        G = 0.0d0
        do d = 1, ndir
          G = G + M(d,1)**2 + M(d,2)**2 + M(d,3)**2
        enddo
        dt = 0.5d0 * vnn / (alpha * G)

        dtmin = min ( dt, dtmin )

      end subroutine compute

  end subroutine Compute_dt


  !> dt-face-k = true: second pass of Compute_dt. The fluxes take k at the face: at a material
  !> jump the face conductivity is the series value of the two half cells, up to twice the k of
  !> the less conductive cell, which the cell k of Compute_dt does not see. Here direction d of
  !> the cell sum takes the larger conductivity k_d of its two faces, each evaluated as its flux
  !> evaluates it:
  !>   dt = vnn rho cp / (2 sum_d k_d |grad xi_d|^2)
  !> A cell whose k_d all equal its own k keeps the step of the cell loop untouched, so one material
  !> at constant k runs bit for bit as with the key off. dtglobal is recomputed over the final steps
  !> (a face step is longer than the cell one where both faces of a direction carry less k).
  subroutine Face_k_dt ( domain, vnn, rampa_iter, dtglobal )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_multigrid
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    real(R8), intent(in)  :: vnn
    integer, intent(in)   :: rampa_iter
    real(R8), intent(out) :: dtglobal
    ! Local
    integer  :: i, j, k, b, level
    real(R8) :: dtcell
    logical  :: face

    level = domain % mg_level
    if ( .not. allocated(conn) ) allocate( conn(obj_multigrid%MGL) )
    if ( .not. allocated(conn(level)%blk) ) call Connection_Faces ( domain, conn(level) )

    dtglobal = domain % dtglobal

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(3) private ( dtcell, face ), reduction ( min : dtglobal )
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)

        call Face_k_Cell_dt ( domain % blk(b), conn(level) % blk(b) % f(i,j,k), i, j, k, vnn, dtcell, face )

        if ( face ) then
          if ( domain%iter < rampa_iter ) dtcell = dtcell * domain%iter / rampa_iter
          domain % blk(b) % dtlocal(i,j,k) = dtcell
        else
          dtcell = domain % blk(b) % dtlocal(i,j,k)
        endif
        dtglobal = min ( dtcell, dtglobal )

      enddo; enddo; enddo
      !$omp end parallel
    enddo

  end subroutine Face_k_dt


  !> Face-conductivity step of cell (i,j,k); face = .false. (dtmin not set) when every k_d equals
  !> the cell k. Interior faces as Diffusive_Flux takes them, at the current T (one material: k at
  !> the mean temperature; two: Series_Kappa). Block-boundary faces: a connection (101) or chimera
  !> (102) face takes the k of BC_Connection's flux for equal half cells, the harmonic mean of the
  !> two sides' k, the neighbour's from the first ghost (its T and matID, which BC_Connection reads
  !> as Tg(1) and matIDg); every other face (walls, symmetry, coupling) the cell k, as the wall
  !> fluxes do.
  subroutine Face_k_Cell_dt ( blk, cf, i, j, k, vnn, dtmin, face )
    use FUSS_Advanced_Types_m
    use FUSS_Global_m
    use FUSS_Lib_Solid
    use FUSS_Lib_Diffusive, only: Series_Kappa
    implicit none
    type(FUSS_block_type), intent(in) :: blk
    integer(I1), intent(in) :: cf               ! connection faces of the cell (Connection_Faces)
    integer, intent(in)     :: i, j, k
    real(R8), intent(in)    :: vnn
    real(R8), intent(out)   :: dtmin
    logical, intent(out)    :: face
    ! Local
    integer  :: d, c(3), n(3)
    real(R8) :: kappa, rho, cp, Gk, dt, kd(3), klo, kup

    call co_Kappa( blk%matID(i,j,k), blk%T(i,j,k), kappa )

    c = [i, j, k]
    do d = 1, ndir
      ! Lower face (2d-1), between n and the cell
      n = c; n(d) = c(d) - 1
      if ( n(d) < 1 ) then
        klo = boundary_k ( 2*d-1 )
      else
        klo = interior_k ( n, c )
      endif
      ! Upper face (2d), between the cell and n
      n = c; n(d) = c(d) + 1
      if ( n(d) > blk%dim(d) ) then
        kup = boundary_k ( 2*d )
      else
        kup = interior_k ( c, n )
      endif
      kd(d) = max ( klo, kup )
    enddo

    face = .not. all ( kd(1:ndir) == kappa )
    if ( .not. face ) return

    call co_Rho( blk%matID(i,j,k), blk%T(i,j,k), rho )
    call co_Cp( blk%matID(i,j,k), blk%T(i,j,k), cp )

    dtmin = 1d8
    Gk = 0.0d0
    do d = 1, ndir
      Gk = Gk + kd(d) * ( blk%M(i,j,k)%c(d,1)**2 + blk%M(i,j,k)%c(d,2)**2 + blk%M(i,j,k)%c(d,3)**2 )
    enddo
    dt = 0.5d0 * vnn * rho * cp / Gk
    dtmin = min ( dt, dtmin )

  contains

    ! k of the interior face in direction d between cells lo and up = lo + e_d, as Diffusive_Flux
    ! evaluates it (lower cell first)
    real(R8) function interior_k ( lo, up )
      integer, intent(in) :: lo(3), up(3)

      if ( blk%matID(lo(1),lo(2),lo(3)) == blk%matID(up(1),up(2),up(3)) ) then
        call co_Kappa( blk%matID(lo(1),lo(2),lo(3)), &
                       0.5d0*( blk%T(lo(1),lo(2),lo(3)) + blk%T(up(1),up(2),up(3)) ), interior_k )
      else
        interior_k = Series_Kappa ( &
          blk%matID(lo(1),lo(2),lo(3)), blk%T(lo(1),lo(2),lo(3)), &
          dot_product ( blk%Minv(lo(1),lo(2),lo(3))%c(:,d), blk%dir(d)%f(lo(1),lo(2),lo(3))%n ), &
          blk%matID(up(1),up(2),up(3)), blk%T(up(1),up(2),up(3)), &
          dot_product ( blk%Minv(up(1),up(2),up(3))%c(:,d), blk%dir(d)%f(lo(1),lo(2),lo(3))%n ) )
      endif

    end function interior_k

    ! k of block-boundary face f, whose first ghost is n
    real(R8) function boundary_k ( f )
      integer, intent(in) :: f
      real(R8) :: kg

      if ( btest ( cf, f-1 ) ) then
        call co_Kappa( blk%matID(n(1),n(2),n(3)), blk%T(n(1),n(2),n(3)), kg )
        if ( kg == kappa ) then
          boundary_k = kappa
        else
          boundary_k = 2d0 * kappa * kg / ( kappa + kg )
        endif
      else
        boundary_k = kappa
      endif

    end function boundary_k

  end subroutine Face_k_Cell_dt


  !> Mark the block-connection (101) and chimera (102) faces of the local cells of one level
  subroutine Connection_Faces ( domain, lev )
    use FUSS_Advanced_Types_m
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    type(conn_level_t), intent(inout)  :: lev
    ! Local
    integer :: b, ii, i

    allocate( lev%blk(domain%nb) )
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      allocate( lev%blk(b)%f( domain%blk(b)%dim(1), domain%blk(b)%dim(2), domain%blk(b)%dim(3) ) )
      lev%blk(b)%f = 0_I1
    enddo

    do ii = 1, domain%n_local_bc
      i = domain%local_bc_idx(ii)
      if ( domain%bc(i)%type /= 101 .and. domain%bc(i)%type /= 102 ) cycle
      associate ( cf => lev%blk(domain%bc(i)%b)%f(domain%bc(i)%i, domain%bc(i)%j, domain%bc(i)%k) )
        cf = ibset ( cf, domain%bc(i)%f - 1 )
      end associate
    enddo

  end subroutine Connection_Faces


  subroutine Set_Global_dt ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        domain % blk(b) % dtlocal(i,j,k) = domain % dtglobal
      enddo; enddo; enddo
      !$omp end parallel
    enddo

  end subroutine Set_Global_dt

end module FUSS_Mod_dt