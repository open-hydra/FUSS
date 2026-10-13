module FUSS_Lib_BC_Fluxes_Wall_RadiationConvection
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Config_Types_m, only: obj_pyrolysis
  use FUSS_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Compute_Wall_Properties

  implicit none
  public

contains

  subroutine BC_Wall_RadiationConvection ( Im, Jm, Km, Fm, Blk, hconv, eps, Tref, Ovar )
    use FUSS_Parameters_m, only: sigma_SB
    use FUSS_Lib_Solid
    use FUSS_Lib_Diffusive
    implicit none
    integer, intent(in)  :: Im, Jm, Km, Fm
    real(R8), intent(in) :: hconv, eps, Tref
    type(FUSS_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(2), intent(out)  :: Ovar
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k
    real(R8) :: Normal(3), Area, M(3,3)
    real(R8) :: T, Gradient(3), Flux, kappa, T_wall, q_wall
    real(R8) :: K, H, qw, err_qw
    integer  :: iter_qw
    ! Convergence controls of the wall-temperature fixed-point iteration (as BC 304)
    integer,  parameter :: iter_qw_max = 100      ! hard iteration cap
    real(R8), parameter :: tol_qw      = 1.d-3    ! [W/m2] on the wall flux
    real(R8), parameter :: dT_min      = 1.d-6    ! [K] neighbourhood of Tref
    real(R8), parameter :: K_min       = 1.d-30   ! floor on the radiative link

    ! Boundary face index
    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    M = Blk % M (Im,Jm,Km) % c

    ! boundary cell variables
    T = Blk % T(Im,Jm,Km)
    if ( obj_pyrolysis%enabled ) then
      kappa = Blk % k_s(Im,Jm,Km)     ! (T, chi) cache, plan 06 section 3.2
    else
      call co_Kappa( Blk % matID(Im,Jm,Km), T, kappa )
    endif

    ! Wall temperature and heat flux: radiation (equivalent conductance K) and
    ! convection (hconv) act in parallel between Tref and the wall, in series
    ! with the half-cell conduction, so q = (K + hconv)(Tref - T_wall) =
    ! (Tref - T)/H. The wall update divides by H (K + hconv): before
    ! 2026-10-13 it divided by H K, as BC 304 does (hconv = 0 there), which put
    ! T_wall too low, K with it, and the flux into the solid 2 % short on the
    ! Zhang (2010) three-code slab (test/pyrolysis/lit/zhang2010).
    T_wall = T  ! Approximate Tw with the local cell temperature
    qw = 0.0d0
    err_qw = 1.d6
    iter_qw = 0
    do while (err_qw > tol_qw)
      iter_qw = iter_qw + 1
      ! Radiative Equivalent Conductivity. (T_wall**4 - Tref**4)/(T_wall - Tref)
      ! is a removable singularity (limit 4 Tref**3), evaluated as the limit
      ! inside a small neighbourhood of Tref.
      if ( abs(T_wall - Tref) > dT_min ) then
        K = sigma_SB*eps*(T_wall**4 - Tref**4)/(T_wall - Tref)
      else
        K = sigma_SB*eps*4.d0*Tref**3
      endif
      if ( K < K_min ) K = K_min
      ! Equivalent convective coefficient
      H = 1.d0/(2.d0*kappa*dot_product(M(Dir,:),Normal)) + 1.d0/(K + hconv)
      ! Fluxes
      Flux = Area * ( 1/H*(Tref - T) )
      ! Wall variables: the whole flux, radiative plus convective
      q_wall  = ( K + hconv ) * ( Tref - T_wall )
      ! Evaluate error and update
      err_qw = abs(qw - q_wall)
      qw = q_wall
      ! Update T_wall
      T_wall = Tref - 1/(H*(K + hconv)) * (Tref - T)
      if ( iter_qw >= iter_qw_max ) then
        write(*,'(A,I0,A,3(1X,I0))') &
          '[WARNING] BC_Wall_RadiationConvection: no convergence in ', iter_qw_max, &
          ' iterations at cell i,j,k =', Im, Jm, Km
        write(*,'(A,ES12.5,A,ES12.5)') &
          '          |dq| = ', err_qw, '  tol = ', tol_qw
        exit
      endif
    enddo

    ! Residual update
    Blk % r (Im,Jm,Km) = Blk % r (Im,Jm,Km) - Flux

    if (present(Ovar)) call Compute_Wall_Properties(Tw=T_wall, qw=q_wall, exit_array=Ovar)

  end subroutine BC_Wall_RadiationConvection

end module FUSS_Lib_BC_Fluxes_Wall_RadiationConvection