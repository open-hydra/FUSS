module FUSS_Lib_BC_Fluxes_Wall_HeatTransfer
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Compute_Wall_Properties

  implicit none
  public

contains

  subroutine BC_Wall_HeatTransfer ( Im, Jm, Km, Fm, Blk, hconv, Tref, qw, Ovar )
    use FUSS_Lib_Solid
    use FUSS_Lib_Diffusive
    implicit none
    integer, intent(in)  :: Im, Jm, Km, Fm
    real(R8), intent(in) :: hconv, Tref, qw
    type(FUSS_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(2), intent(out)  :: Ovar
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k
    real(R8) :: Normal(3), Area, M(3,3)
    real(R8) :: T, Gradient(3), Flux, kappa, T_wall, q_wall, H, Hrad

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

    ! Difference in Across-face direction: symmetry; tangential direction: null
    H = 1d0/(2d0*kappa*dot_product(M(Dir,:),Normal)) + 1d0/hconv
    Hrad = hconv/(2d0*kappa*dot_product(M(Dir,:),Normal)) + 1
    ! Fluxes
    Flux = Area * ( 1/H*(Tref - T) + 1/Hrad*qw )
    ! Wall variables
    T_wall = (hconv*Tref + T*2d0*kappa*dot_product(M(Dir,:),Normal) + qw) / ( hconv + 2d0*kappa*dot_product(M(Dir,:),Normal) )
    q_wall = hconv * ( Tref - T_wall ) + qw

    ! Residual update
    Blk % r (Im,Jm,Km) = Blk % r (Im,Jm,Km) - Flux

    if (present(Ovar)) call Compute_Wall_Properties(Tw=T_wall, qw=q_wall, exit_array=Ovar)

  end subroutine BC_Wall_HeatTransfer

end module FUSS_Lib_BC_Fluxes_Wall_HeatTransfer