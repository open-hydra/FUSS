module FUSS_Lib_Diffusive
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Diffusive_Flux, Compute_Diffusive_Flux, Inverse3

contains

  subroutine Diffusive_Flux ( matID, normal, area, T1, T2, T3, T4, T5, T6, T7, &
                              T8, T9, T10, M1, M2, Res1, Res2, aa, bb, cc, inverse_metric )
    use FUSS_Global_m
    use FUSS_Lib_Solid
    implicit none
    real(R8), intent(in) :: matID
    integer, intent(in) :: aa, bb, cc
    real(R8), intent(in) :: normal(3), area
    real(R8), intent(in) :: T1, T2, T3, T4, T5
    real(R8), intent(in) :: T6, T7, T8, T9, T10
    real(R8), intent(in), dimension(3,3) :: M1, M2   ! cell metrics, or their inverses (see inverse_metric)
    logical, intent(in) :: inverse_metric           ! M1, M2 are the cells' inverse metrics (diffusive-metric = inverse-mean)
    real(R8), intent(inout) :: Res1, Res2
    ! Local
    real(R8) :: Gradient(3), T, M(3,3), Flux, kappa

    ! Gradient in the same direction of the face: 1 and 2
    Gradient ( aa ) = T2 - T1

    ! Gradient in tangential directions: 3-10
    call Tangential_Gradient ( T3, T4, T5, T6,  Gradient( bb ) )
    call Tangential_Gradient ( T7, T8, T9, T10, Gradient( cc ) )

    if (ndir==2) Gradient(cc) = 0.0
    
    ! Metric tensor at the face.
    ! inverse-mean: the mean of the two inverse metrics (cell edge vectors, stored
    ! per cell in blk%Minv), inverted, so that the one-index difference across the
    ! face is divided by the centre-to-centre distance.
    ! mean: averaging M itself divides it by the mean of 1/h instead, which
    ! overstates the gradient by (1+r)^2/(4r) for a cell-size ratio r across the
    ! face (+4% at r = 1.5, +9% at 1.8), at any mesh resolution.
    if ( inverse_metric ) then
      M = Inverse3 ( 0.5d0 * ( M1 + M2 ) )
    else
      M = 0.5d0 * ( M1 + M2 )
    end if
    Gradient = matmul ( Gradient, M )

    T = 0.5d0 * ( T1 + T2 )
    call co_Kappa( matID, T, kappa )

    call Compute_Diffusive_Flux ( kappa, Gradient, area, normal, Flux )

    Res1 = Res1 - Flux
    Res2 = Res2 + Flux

  end subroutine Diffusive_Flux


  !> Inverse of a 3x3 matrix by its adjugate (transposed cofactors). Works for M
  !> as stored by Compute_Metric_Tensor (cofactors/det = transposed inverse of the
  !> edge vectors): transposition commutes with averaging and inversion.
  pure function Inverse3 ( A ) result ( B )
    implicit none
    real(R8), intent(in) :: A(3,3)
    real(R8)             :: B(3,3)
    real(R8) :: det

    B(1,1) =  ( A(2,2)*A(3,3) - A(2,3)*A(3,2) )
    B(1,2) = -( A(1,2)*A(3,3) - A(1,3)*A(3,2) )
    B(1,3) =  ( A(1,2)*A(2,3) - A(1,3)*A(2,2) )
    B(2,1) = -( A(2,1)*A(3,3) - A(2,3)*A(3,1) )
    B(2,2) =  ( A(1,1)*A(3,3) - A(1,3)*A(3,1) )
    B(2,3) = -( A(1,1)*A(2,3) - A(1,3)*A(2,1) )
    B(3,1) =  ( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
    B(3,2) = -( A(1,1)*A(3,2) - A(1,2)*A(3,1) )
    B(3,3) =  ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )
    det = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
    B = B / det

  end function Inverse3


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