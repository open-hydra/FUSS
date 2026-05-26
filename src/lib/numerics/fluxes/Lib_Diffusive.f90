module FUSS_Lib_Diffusive
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Diffusive_Flux, Compute_Diffusive_Flux

contains

  subroutine Diffusive_Flux ( matID, normal, area, T1, T2, T3, T4, T5, T6, T7, &
                              T8, T9, T10, M1, M2, Res1, Res2, aa, bb, cc )
    use FUSS_Global_m
    use FUSS_Lib_Solid
    implicit none
    real(R8), intent(in) :: matID
    integer, intent(in) :: aa, bb, cc
    real(R8), intent(in) :: normal(3), area
    real(R8), intent(in) :: T1, T2, T3, T4, T5
    real(R8), intent(in) :: T6, T7, T8, T9, T10
    real(R8), intent(in), dimension(3,3) :: M1, M2
    real(R8), intent(inout) :: Res1, Res2
    ! Local
    real(R8) :: Gradient(3), T, M(3,3), Flux, kappa

    ! Gradient in the same direction of the face: 1 and 2
    Gradient ( aa ) = T2 - T1

    ! Gradient in tangential directions: 3-10
    call Tangential_Gradient ( T3, T4, T5, T6,  Gradient( bb ) )
    call Tangential_Gradient ( T7, T8, T9, T10, Gradient( cc ) )

    if (ndir==2) Gradient(cc) = 0.0
    
    M = 0.5d0 * ( M1 + M2 )
    Gradient = matmul ( Gradient, M )

    T = 0.5d0 * ( T1 + T2 )
    call co_Kappa( matID, T, kappa )

    call Compute_Diffusive_Flux ( kappa, Gradient, area, normal, Flux )

    Res1 = Res1 - Flux
    Res2 = Res2 + Flux

  end subroutine Diffusive_Flux


  subroutine Tangential_Gradient ( T1, T2, T3, T4, Gradient )
    implicit none
    real(R8), intent(in)  :: T1, T2, T3, T4
    real(R8), intent(out) :: Gradient

    Gradient = ( T2 - T1 + T4 - T3 ) * 0.25d0

  end subroutine Tangential_Gradient



  subroutine Compute_Diffusive_Flux ( kappa, Gradient, area, normal, Flux )
    implicit none
    real(R8), intent(in)  :: kappa, Gradient(3), area, normal(3)
    real(R8), intent(out) :: Flux

    ! Diffusive flux
    Flux = area * kappa * dot_product ( Gradient, normal )

  end subroutine Compute_Diffusive_Flux
  
end module FUSS_Lib_Diffusive