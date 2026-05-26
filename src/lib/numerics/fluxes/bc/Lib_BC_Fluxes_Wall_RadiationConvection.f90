module FUSS_Lib_BC_Fluxes_Wall_RadiationConvection
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
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

    ! Boundary face index
    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    M = Blk % M (Im,Jm,Km) % c

    ! boundary cell variables
    T = Blk % T(Im,Jm,Km)
    call co_Kappa( Blk % matID(Im,Jm,Km), T, kappa )  

    ! Loop to evaluate Wall Temperature and Heat Flux
    T_wall = T  ! Approximate Tw with the local cell temperature
    qw = 0.0d0
    err_qw = 1.d6
    do while (err_qw > 1.d-3)
      ! Radiative Equivalent Conductivity
      K = sigma_SB*eps*(T_wall**4 - Tref**4)/(T_wall - Tref)
      ! Equivalent convective coefficient
      H = 1.d0/(2.d0*kappa*dot_product(M(Dir,:),Normal)) + 1.d0/(K + hconv)
      ! Fluxes
      Flux = Area * ( 1/H*(Tref - T) )
      ! Wall variables
      q_wall  = K * ( Tref - T_wall )
      ! Evaluate error and update
      err_qw = abs(qw - q_wall)
      qw = q_wall
      ! Update T_wall
      T_wall = Tref - 1/(H*K) * (Tref - T)
    enddo

    ! Residual update
    Blk % r (Im,Jm,Km) = Blk % r (Im,Jm,Km) - Flux

    if (present(Ovar)) call Compute_Wall_Properties(Tw=T_wall, qw=q_wall, exit_array=Ovar)

  end subroutine BC_Wall_RadiationConvection

end module FUSS_Lib_BC_Fluxes_Wall_RadiationConvection