module FUSS_Lib_BC_Fluxes_Wall_Temperature
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Compute_Wall_Properties

  implicit none
  public

contains

  subroutine BC_Wall_Temperature ( Im, Jm, Km, Fm, Blk, Tw, Ovar )
    use FUSS_Lib_Solid
    use FUSS_Lib_Diffusive
    implicit none
    integer, intent(in)  :: Im, Jm, Km, Fm
    real(R8), intent(in) :: Tw
    type(FUSS_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(2), intent(out)  :: Ovar
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k
    real(R8) :: Normal(3), Area, M(3,3)
    real(R8) :: T, Gradient(3), Flux, kappa, T_wall, q_wall

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
    Gradient = 0d0

    Gradient(Dir) = ( T - Tw ) * modfm3

    if (Fm == 1 .or. Fm == 2) then
      Gradient(2) = (Blk % T(Im,Jm-1,Km) - Blk % T(Im,Jm+1,Km))/2.0
      Gradient(3) = (Blk % T(Im,Jm,Km-1) - Blk % T(Im,Jm,Km+1))/2.0
    elseif (Fm == 3 .or. Fm == 4) then
      Gradient(1) = (Blk % T(Im-1,Jm,Km) - Blk % T(Im+1,Jm,Km))/2.0
      Gradient(3) = (Blk % T(Im,Jm,Km-1) - Blk % T(Im,Jm,Km+1))/2.0
    elseif (Fm == 5 .or. Fm == 6) then
      Gradient(1) = (Blk % T(Im-1,Jm,Km) - Blk % T(Im+1,Jm,Km))/2.0
      Gradient(2) = (Blk % T(Im,Jm-1,Km) - Blk % T(Im,Jm+1,Km))/2.0
    endif

    if (ndir == 2) Gradient(3) = 0.0
    Gradient = matmul ( Gradient, M )

    ! Fluxes
    call Compute_Diffusive_Flux ( kappa, Gradient, Area, Normal, Flux )

    ! Residual update
    Blk % r (Im,Jm,Km) = Blk % r (Im,Jm,Km) - modfm2 * Flux

    ! Wall variables
    T_wall = Tw
    q_wall = kappa * dot_product( Gradient(:), Normal )
    if (present(Ovar)) call Compute_Wall_Properties(Tw=T_wall, qw=q_wall, exit_array=Ovar)

  end subroutine BC_Wall_Temperature

end module FUSS_Lib_BC_Fluxes_Wall_Temperature