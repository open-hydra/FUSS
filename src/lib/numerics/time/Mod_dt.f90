module FUSS_Mod_dt
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Compute_dt, Set_Global_dt

contains

  subroutine Compute_dt ( domain, vnn, rampa_iter )
    use FUSS_Advanced_Types_m
    use FUSS_Global_m
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

    ! MPI: global minimum across all ranks
    call mpi_allreduce_min_r8(dtglobal, dtglobal_mpi)
    domain % dtglobal = dtglobal_mpi

    contains
      
      ! Explicit limit of the diffusion operator: alpha dt sum_d |grad xi_d|^2 <= vnn / 2, with grad xi_d
      ! (xi_d the cell index along d) the row d of the cell metric M. The largest eigenvalue of the
      ! discrete operator on a uniform mesh is 4 alpha sum_d |grad xi_d|^2, also on a skewed one: the
      ! mixed-derivative terms vanish in that mode. On an orthogonal cell |grad xi_d| = 1/dl_d; on a cell
      ! skewed by an angle phi it is larger (1/(dl_d sin phi) for a parallelogram), and the edge lengths
      ! dl used before overstated the stable step by 1/sin^2(phi), 2 at 45 degrees.
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