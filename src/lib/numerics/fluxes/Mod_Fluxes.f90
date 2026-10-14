module FUSS_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Fluxes

contains

  subroutine Fluxes ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_space_scheme
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: b

    do b = 1, domain % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle
      ! face-k = cell-mean: k of every interior cell at its own temperature, once per stage.
      ! The unallocated kc cannot be passed: the default mode hands T as a stand-in the
      ! callee never reads as a conductivity.
      if ( obj_space_scheme % cell_mean_k ) then
        call Fill_Cell_K ( domain % blk(b) % T, domain % blk(b) % matID, &
                           domain % blk(b) % kc, domain % blk(b) % dim )
      end if
      ! Cell metrics for the diffusive face metric: inverse metrics (inverse-mean) or metrics (mean)
      if ( obj_space_scheme % inverse_metric ) then
        call Fluxes_blk ( domain % blk(b) % T,   &
                          domain % blk(b) % matID, &
                          domain % blk(b) % r,   &
                          domain % blk(b) % dir, &
                          domain % blk(b) % Minv, &
                          domain % blk(b) % dim, &
                          .true., &
                          merge_kc ( domain % blk(b), obj_space_scheme % cell_mean_k ), &
                          obj_space_scheme % cell_mean_k )
      else
        call Fluxes_blk ( domain % blk(b) % T,   &
                          domain % blk(b) % matID, &
                          domain % blk(b) % r,   &
                          domain % blk(b) % dir, &
                          domain % blk(b) % m,   &
                          domain % blk(b) % dim, &
                          .false., &
                          merge_kc ( domain % blk(b), obj_space_scheme % cell_mean_k ), &
                          obj_space_scheme % cell_mean_k )
      end if
    enddo

  contains

    !> blk%kc when face-k = cell-mean, else blk%T as an inert stand-in (kc is unallocated)
    function merge_kc ( blk, cell_mean_k ) result ( kc )
      implicit none
      type(FUSS_block_type), intent(in), target :: blk
      logical, intent(in) :: cell_mean_k
      real(R8), pointer :: kc(:,:,:)
      if ( cell_mean_k ) then
        kc => blk % kc
      else
        kc => blk % T
      end if
    end function merge_kc

  end subroutine Fluxes


  !> face-k = cell-mean: k of every interior cell at its own temperature and material
  subroutine Fill_Cell_K ( T, matID, kc, n )
    use FUSS_Global_m, only: gc
    use FUSS_Lib_Solid, only: co_Kappa
    implicit none
    integer, intent(in) :: n(3)
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in)  :: T, matID
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(out) :: kc
    ! Local
    integer :: i, j, k

    !$omp do collapse (2)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      call co_Kappa ( matID(i,j,k), T(i,j,k), kc(i,j,k) )
    enddo; enddo; enddo

  end subroutine Fill_Cell_K


  subroutine Fluxes_blk ( T, matID, Res, Dir, M, n, inverse_metric, kc, cell_mean_k )
    use FUSS_Base_Types_m
    use FUSS_Global_m, only: gc
    use FUSS_Lib_Diffusive
    implicit none
    integer, intent(in) :: n(3)
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in)  :: matID
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in)  :: T
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(out) :: Res
    type(FUSS_d_metrics_type), dimension(3), intent(in) :: Dir
    type(FUSS_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M   ! cell metrics, or their inverses
    logical, intent(in) :: inverse_metric   ! M holds the inverse metrics (diffusive-metric = inverse-mean)
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: kc   ! per-cell k (face-k = cell-mean; otherwise never read)
    logical, intent(in) :: cell_mean_k
    ! Local
    integer :: i, j, k

    ! -----------------------------------------------------------------
    ! Reset residuals
    !$omp do collapse (3)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      Res(i,j,k) = 0d0
    enddo; enddo; enddo

    !$omp do collapse (2)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1) - 1
      call Diffusive_Flux ( matID(i,j,k), &
                            Dir(1) % f(i,j,k) % n,  &
                            Dir(1) % f(i,j,k) % A,  &
                            T(i  ,j,k),        &
                            T(i+1,j,k),        &
                            T(i  ,j-1,k),      &
                            T(i  ,j+1,k),      &
                            T(i+1,j-1,k),      &
                            T(i+1,j+1,k),      &
                            T(i  ,j,k-1),      &
                            T(i  ,j,k+1),      &
                            T(i+1,j,k-1),      &
                            T(i+1,j,k+1),      &
                            M(i  ,j,k) % c,    &
                            M(i+1,j,k) % c,    &
                            Res(i  ,j,k),      &
                            Res(i+1,j,k),      &
                            1, 2, 3, inverse_metric, &
                            cell_mean_k, 0.5d0 * ( kc(i,j,k) + kc(i+1,j,k) ) )
    enddo; enddo; enddo

    !$omp do collapse (2)
    do k = 1, n(3)
    do i = 1, n(1)
    do j = 1, n(2) - 1
      call Diffusive_Flux ( matID(i,j,k), &
                            Dir(2) % f(i,j,k) % n,  &
                            Dir(2) % f(i,j,k) % A,  &
                            T(i,j  ,k),        &
                            T(i,j+1,k),        &
                            T(i-1,j  ,k),      &
                            T(i+1,j  ,k),      &
                            T(i-1,j+1,k),      &
                            T(i+1,j+1,k),      &
                            T(i,j  ,k-1),      &
                            T(i,j  ,k+1),      &
                            T(i,j+1,k-1),      &
                            T(i,j+1,k+1),      &
                            M(i,j  ,k) % c,    &
                            M(i,j+1,k) % c,    &
                            Res(i,j  ,k),      &
                            Res(i,j+1,k),      &
                            2, 1, 3, inverse_metric, &
                            cell_mean_k, 0.5d0 * ( kc(i,j,k) + kc(i,j+1,k) ) )
    enddo; enddo; enddo

    !$omp do collapse (2)
    do j = 1, n(2)
    do i = 1, n(1)
    do k = 1, n(3) - 1
      call Diffusive_Flux ( matID(i,j,k), &
                            Dir(3) % f(i,j,k) % n,  &
                            Dir(3) % f(i,j,k) % A,  &
                            T(i,j,k  ),        &
                            T(i,j,k+1),        &
                            T(i-1,j,k  ),      &
                            T(i+1,j,k  ),      &
                            T(i-1,j,k+1),      &
                            T(i+1,j,k+1),      &
                            T(i,j-1,k  ),      &
                            T(i,j+1,k  ),      &
                            T(i,j-1,k+1),      &
                            T(i,j+1,k+1),      &
                            M(i,j,k  ) % c,    &
                            M(i,j,k+1) % c,    &
                            Res(i,j,k  ),      &
                            Res(i,j,k+1),      &
                            3, 1, 2, inverse_metric, &
                            cell_mean_k, 0.5d0 * ( kc(i,j,k) + kc(i,j,k+1) ) )
    enddo; enddo; enddo

  end subroutine Fluxes_blk

end module FUSS_Mod_Fluxes