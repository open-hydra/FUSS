module FUSS_Mod_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none

contains

  subroutine Setup_Metrics ( domain )
    use FUSS_Base_Types_m
    use FUSS_Advanced_Types_m
    use FUSS_Global_m
    use FUSS_Lib_Metrics
    use FUSS_Lib_Diffusive, only: Inverse3
    use FUSS_Config_Types_m, only: obj_space_scheme
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k
    logical :: inverse_metric
    integer :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(FUSS_vector_3D_type) :: N1, N2, N3, N4, N5, N6, N7, N8

    call Check_Mesh_Type ( domain )

    ! Inverse metric of the interior cells, for the diffusive face metric (diffusive-metric =
    ! inverse-mean): the face then needs one 3x3 inversion instead of three. Same bounds as M,
    ! because Fluxes_blk takes either array through the same explicit-shape argument.
    inverse_metric = obj_space_scheme % inverse_metric

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      if ( inverse_metric .and. .not. allocated(domain % blk(b) % Minv) ) &
        allocate( domain % blk(b) % Minv, mold = domain % blk(b) % M )
      !$omp parallel
      ! Compute metric tensor and cell dimension across i,j,k. => block % M, & block % dl      
      !$omp do collapse(3) private(i, j, k, N1, N2, N3, N4, N5, N6, N7, N8)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        N1 % c = domain % blk(b) % node(i-1,j-1,k-1) % c
        N2 % c = domain % blk(b) % node(i-1,j-1,k  ) % c
        N3 % c = domain % blk(b) % node(i-1,j  ,k-1) % c
        N4 % c = domain % blk(b) % node(i-1,j  ,k  ) % c
        N5 % c = domain % blk(b) % node(i  ,j-1,k-1) % c
        N6 % c = domain % blk(b) % node(i  ,j-1,k  ) % c
        N7 % c = domain % blk(b) % node(i  ,j  ,k-1) % c
        N8 % c = domain % blk(b) % node(i  ,j  ,k  ) % c
        call Compute_Metric_Tensor ( N1, N2, N3, N4, N5, N6, N7, N8, domain % blk(b) % M(i,j,k), domain % blk(b) % dl(i,j,k), domain % blk(b) % vol(i,j,k) )
        if ( inverse_metric ) domain % blk(b) % Minv(i,j,k) % c = Inverse3 ( domain % blk(b) % M(i,j,k) % c )
      enddo; enddo; enddo
      !$omp end parallel

      ! Compute Normal & Area
      call Compute_Norm_Area ( domain % blk(b) )
    enddo
    
    !$omp parallel
    ! Create the nodes for gc layers of ghost cell
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, domain % nbound
      Bm = domain % bc(i) % b
      Im = domain % bc(i) % i 
      Jm = domain % bc(i) % j 
      Km = domain % bc(i) % k 
      Fm = domain % bc(i) % f
      select case ( domain % bc(i) % type )
        case(101) ! block connection
          Bs = domain % bc(i) % bs
          Is = domain % bc(i) % is 
          Js = domain % bc(i) % js 
          Ks = domain % bc(i) % ks 
          Fs = domain % bc(i) % fs
          d11s = domain % bc(i) % d11
          d12s = domain % bc(i) % d12
          d21s = domain % bc(i) % d21
          d22s = domain % bc(i) % d22
          call BC_Connect_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                    Is, Js, Ks, Fs, domain % blk(Bs), d11s, d12s, d21s, d22s, &
                                    domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg)
        case(300) ! symmetry
          call BC_Symmetry_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                     domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg )
        case default  ! Rientra qui anche chimera
          call BC_Extrapolate_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                        domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg )
      end select
    enddo
    !$omp end parallel

  end subroutine Setup_Metrics

end module FUSS_Mod_Metrics