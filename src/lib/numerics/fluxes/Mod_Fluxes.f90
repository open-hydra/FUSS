module FUSS_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Fluxes

contains

  subroutine Fluxes ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Mod_MPI, only: is_local_block
    use FUSS_Config_Types_m, only: obj_space_scheme
    use FUSS_Lib_BC_Fluxes_Connection, only: Connection_Fluxes
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: b

    do b = 1, domain % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle
      ! The metrics of the cells build the diffusive face metric
      call Fluxes_blk ( domain % blk(b) % T,   &
                        domain % blk(b) % matID, &
                        domain % blk(b) % r,   &
                        domain % blk(b) % dir, &
                        domain % blk(b) % M, &
                        domain % blk(b) % fm, &
                        obj_space_scheme % stencil_metric, &
                        domain % blk(b) % dim )
    enddo

    ! One flux per block connection (101), computed by the record that owns the face
    call Connection_Fluxes ( domain )

  end subroutine Fluxes


  subroutine Fluxes_blk ( T, matID, Res, Dir, M, Mfd, stencil, n )
    use FUSS_Base_Types_m
    use FUSS_Advanced_Types_m, only: face_metric_type
    use FUSS_Global_m, only: gc
    use FUSS_Lib_Diffusive
    implicit none
    integer, intent(in) :: n(3)
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in)  :: matID
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in)  :: T
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(out) :: Res
    type(FUSS_d_metrics_type), dimension(3), intent(in) :: Dir
    type(FUSS_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M   ! cell metrics
    type(face_metric_type), dimension(3), intent(in) :: Mfd   ! stencil face metrics
    logical, intent(in) :: stencil   ! diffusive-metric = stencil (else the mean of the two cell metrics)
    ! Local
    integer :: i, j, k
    real(R8) :: Mf(3,3)

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
      if ( stencil ) then
        Mf = Mfd(1) % Mf(i,j,k) % c
      else
        Mf = 0.5d0 * ( M(i  ,j,k) % c + M(i+1,j,k) % c )
      end if
      call Diffusive_Flux ( matID(i,j,k), matID(i+1,j,k), &
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
                            Mf, &
                            M(i  ,j,k) % c,    &
                            M(i+1,j,k) % c,    &
                            Res(i  ,j,k),      &
                            Res(i+1,j,k),      &
                            1, 2, 3 )
    enddo; enddo; enddo

    !$omp do collapse (2)
    do k = 1, n(3)
    do i = 1, n(1)
    do j = 1, n(2) - 1
      if ( stencil ) then
        Mf = Mfd(2) % Mf(i,j,k) % c
      else
        Mf = 0.5d0 * ( M(i,j  ,k) % c + M(i,j+1,k) % c )
      end if
      call Diffusive_Flux ( matID(i,j,k), matID(i,j+1,k), &
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
                            Mf, &
                            M(i,j  ,k) % c,    &
                            M(i,j+1,k) % c,    &
                            Res(i,j  ,k),      &
                            Res(i,j+1,k),      &
                            2, 1, 3 )
    enddo; enddo; enddo

    !$omp do collapse (2)
    do j = 1, n(2)
    do i = 1, n(1)
    do k = 1, n(3) - 1
      if ( stencil ) then
        Mf = Mfd(3) % Mf(i,j,k) % c
      else
        Mf = 0.5d0 * ( M(i,j,k  ) % c + M(i,j,k+1) % c )
      end if
      call Diffusive_Flux ( matID(i,j,k), matID(i,j,k+1), &
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
                            Mf, &
                            M(i,j,k  ) % c,    &
                            M(i,j,k+1) % c,    &
                            Res(i,j,k  ),      &
                            Res(i,j,k+1),      &
                            3, 1, 2 )
    enddo; enddo; enddo

  end subroutine Fluxes_blk

end module FUSS_Mod_Fluxes