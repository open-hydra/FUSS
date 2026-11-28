!> Face kernel of the Darcy pressure equation and of its gas-enthalpy flux
!> (plan 07 sections 3.1, 3.3, 3.5): the tensor analogue of Diffusive_Flux.
!>
!> With the index-space differences G (G(aa) = p2 - p1 across the face, the
!> tangential G(bb), G(cc) from the 4+4 stencil as Tangential_Gradient), the
!> face metric M_f = (M1 + M2)/2 and the physical gradient grad p = M_f^T G
!> (Diffusive_Flux's `matmul(Gradient, M)`), the flux of the symmetric face
!> tensor Gf through the face of normal n and area A is
!>     F = A n . (Gf grad p) = A sum_l G_l [ M_f(l,:) . (Gf n) ]
!> = C (p2 - p1) + E, with C = A M_f(aa,:).(Gf n) the conductance of the
!> normal difference and E the tangential part. F is the gas (or enthalpy)
!> flow INTO the lower-index cell 1, as Diffusive_Flux's Flux is the heat
!> flow into it (Res1 -= Flux, Res2 += Flux). Splitting C from E is what lets
!> the same kernel feed the line-implicit rows (C on the diagonal, E lagged),
!> the Jacobi residual and the explicit energy flux.
module FUSS_Lib_Diffusive_Tensor
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Darcy_Face, Darcy_Boundary_Face

contains

  !> Interior face. p3..p6 / p7..p10 are the tangential pairs in the order
  !> Fluxes_blk passes them to Diffusive_Flux (cell-1 side then cell-2 side
  !> for direction bb, then the same for cc). ndir == 2 drops G(cc).
  pure subroutine Darcy_Face ( Gf, normal, area, p3, p4, p5, p6, p7, p8, p9, p10, M1, M2, aa, bb, cc, C, E )
    use FUSS_Global_m, only: ndir
    implicit none
    real(R8), intent(in)  :: Gf(3,3), normal(3), area
    real(R8), intent(in)  :: p3, p4, p5, p6, p7, p8, p9, p10
    real(R8), intent(in)  :: M1(3,3), M2(3,3)
    integer,  intent(in)  :: aa, bb, cc
    real(R8), intent(out) :: C, E
    ! Local
    real(R8) :: Mf(3,3), w(3), G(3)

    Mf = 0.5d0 * ( M1 + M2 )
    w  = matmul( Gf, normal )
    G  = 0.0d0
    G(bb) = ( p4 - p3 + p6 - p5 ) * 0.25d0
    G(cc) = ( p8 - p7 + p10 - p9 ) * 0.25d0
    if ( ndir == 2 ) G(cc) = 0.0d0

    C = area * dot_product( Mf(aa,:), w )
    E = area * ( G(bb) * dot_product( Mf(bb,:), w ) + G(cc) * dot_product( Mf(cc,:), w ) )

  end subroutine Darcy_Face


  !> Boundary face of cell c (half-cell closure, plan 07 section 3.4): with
  !> the face tensor Gf, the cell metric M, the normal n of the face (pointing
  !> in the +index direction as every FUSS face normal does) and the cell-row
  !> central differences Gt (index space, the two tangential components of the
  !> face direction aa; zero in the normal one), the flow INTO the cell is
  !>     inflow = Cw (q - p_c) + s E,  Cw = 2 A M(aa,:).(Gf n),  s = -1 (low face) / +1 (high face),
  !> q the value on the other side (p_w of an open face, the ghost value of a
  !> connection). The factor 2 is the half-cell distance in index space, as
  !> BC_Wall_HeatTransfer's `2 kappa M(Dir,:).n`.
  pure subroutine Darcy_Boundary_Face ( Gf, normal, area, M, Gt, aa, Cw, E )
    use FUSS_Global_m, only: ndir
    implicit none
    real(R8), intent(in)  :: Gf(3,3), normal(3), area, M(3,3), Gt(3)
    integer,  intent(in)  :: aa
    real(R8), intent(out) :: Cw, E
    ! Local
    real(R8) :: w(3), G(3)
    integer  :: l

    w = matmul( Gf, normal )
    G = Gt
    G(aa) = 0.0d0
    if ( ndir == 2 ) G(3) = 0.0d0
    Cw = 2.0d0 * area * dot_product( M(aa,:), w )
    E  = 0.0d0
    do l = 1, 3
      if ( l /= aa ) E = E + area * G(l) * dot_product( M(l,:), w )
    enddo

  end subroutine Darcy_Boundary_Face

end module FUSS_Lib_Diffusive_Tensor
