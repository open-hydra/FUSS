module FUSS_Lib_BC_Fluxes_Symmetry
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Lib_BC_Fluxes, only: Face_Index

  implicit none
  public

contains

  subroutine BC_Symmetry ( Im, Jm, Km, Fm, Blk )
    use FUSS_Lib_Solid
    use FUSS_Lib_Diffusive
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout) :: Blk
    ! Local
    integer :: modfm2, dir, Face_i, Face_j, Face_k
    real(R8) :: Normal(3), Area, M(3,3)
    real(R8) :: T_ip, T_im, T_jp, T_jm, T_kp, T_km
    real(R8) :: T, Gradient(3), Flux, kappa

    ! Boundary face index
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    modfm2 = 1 - 2 * mod (Fm,2)

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    M = Blk % M (Im,Jm,Km) % c

    ! Stencil building
    T = Blk % T(Im,Jm,Km)
    T_ip  = Blk % T (Im+1,Jm,Km)
    T_im  = Blk % T (Im-1,Jm,Km)
    T_jp  = Blk % T (Im,Jm+1,Km)
    T_jm  = Blk % T (Im,Jm-1,Km)
    T_kp  = Blk % T (Im,Jm,Km+1)
    T_km  = Blk % T (Im,Jm,Km-1)

    ! Gradient computations
    select case ( Fm )
      case ( 1 : 2 )
        Gradient (1) = 0.0d0
        Gradient (2) = 0.5d0 * ( T_jp - T_jm ) * modfm2
        Gradient (3) = 0.5d0 * ( T_kp - T_km ) * modfm2
      case ( 3 : 4 )
        Gradient (1) = 0.5d0 * ( T_ip - T_im ) * modfm2
        Gradient (2) = 0.0d0
        Gradient (3) = 0.5d0 * ( T_kp - T_km ) * modfm2
      case( 5 : 6 )
        Gradient (1) = 0.5d0 * ( T_ip - T_im ) * modfm2 
        Gradient (2) = 0.5d0 * ( T_jp - T_jm ) * modfm2
        Gradient (3) = 0.0d0
    end select

    if (ndir == 2) Gradient(3) = 0.0
    Gradient = matmul ( Gradient, M )

    call co_Kappa( Blk % matID(Im,Jm,Km), T, kappa )
    call Compute_Diffusive_Flux ( kappa, Gradient, Area, Normal, Flux )

    ! Residual update
    Blk % r(Im,Jm,Km) = Blk % r(Im,Jm,Km) - modfm2 * Flux
  
  end subroutine BC_Symmetry

end module FUSS_Lib_BC_Fluxes_Symmetry