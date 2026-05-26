module FUSS_Lib_BC_Fluxes_Wall_Heat
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Lib_BC_Fluxes, only: Face_Index, Compute_Wall_Properties

  implicit none
  public

contains

  subroutine BC_Wall_Heat ( Im, Jm, Km, Fm, Blk, qw, Ovar )
    use FUSS_Lib_Solid
    implicit none
    integer, intent(in)  :: Im, Jm, Km, Fm
    real(R8), intent(in) :: qw
    type(FUSS_block_type), intent(inout)  :: Blk
    real(R8), optional, dimension(2), intent(inout) :: Ovar
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k
    real(R8) :: Normal(3), Area, M(3,3)
    real(R8) :: T, Gradient(3), Flux, kappa, T_wall, q_wall

    ! Boundary face index
    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    modfm2 = 1 - 2 * mod (Fm,2)

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    M = Blk % M (Im,Jm,Km) % c

    ! boundary cell variables
    T = Blk % T(Im,Jm,Km)
    call co_kappa( Blk % matID(Im,Jm,Km), T, kappa )

    ! Difference in Across-face direction: symmetry; tangential direction: null
    Gradient = 0d0
    Gradient(Dir) = - qw / ( 2d0*kappa*dot_product( M(Dir,:), Normal ) ) * modfm2
    if (ndir == 2) Gradient(3) = 0.0
    Gradient = matmul ( Gradient, M )

    ! Fluxes
    Flux = Area * qw

    ! Residual update
    Blk % r (Im,Jm,Km) = Blk % r (Im,Jm,Km) - modfm2 * Flux

    ! Wall variables
    T_wall = T - Gradient(Dir)
    q_wall = qw
    if (present(Ovar)) call Compute_Wall_Properties(Tw=T_wall, qw=q_wall, exit_array=Ovar)

  end subroutine BC_Wall_Heat

end module FUSS_Lib_BC_Fluxes_Wall_Heat