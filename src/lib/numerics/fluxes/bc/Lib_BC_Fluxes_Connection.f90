module FUSS_Lib_BC_Fluxes_Connection
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Lib_BC_Fluxes, only: Face_Index

  implicit none
  public

contains

  subroutine BC_Connection ( Im, Jm, Km, Fm, Blk, Mg, Tg, matIDg )
    use FUSS_Lib_Solid
    use FUSS_Lib_Diffusive
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_tensor_3D_type), intent(in) :: Mg
    real(R8), intent(in)         :: Tg(6), matIDg
    type(FUSS_block_type), intent(inout) :: Blk
    ! Local
    integer :: Dir, Face_im, Face_jm, Face_km, modfm2
    real(R8) :: T_loc, T_con, T_ip_loc, T_im_loc, T_jp_loc, T_jm_loc, T_kp_loc, T_km_loc
    real(R8) :: Normal(3), Area, M_loc(3,3), M_con(3,3), Gradient(3)
    real(R8) :: kappa_loc, kappa_con, Num_loc, Num_con, Den_loc, Den_con
    real(R8) :: Tint, Flux
    
    ! Boundary face index
    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_im, Face_jm, Face_km )
    modfm2 = 1 - 2 * mod (Fm,2)

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_im,Face_jm,Face_km) % N
    Area   = Blk % dir(Dir) % f(Face_im,Face_jm,Face_km) % A
    M_loc  = Blk % M(Im,Jm,Km) % c
    M_con  = Mg % c

    ! Stencil building
    T_loc    = Blk % T (Im,Jm,Km)
    T_ip_loc = Blk % T (Im+1,Jm,Km)
    T_im_loc = Blk % T (Im-1,Jm,Km)
    T_jp_loc = Blk % T (Im,Jm+1,Km)
    T_jm_loc = Blk % T (Im,Jm-1,Km)
    if (ndir == 2) then
      T_kp_loc = huge(0.0)
      T_km_loc = huge(0.0)
    else
      T_kp_loc = Blk % T (Im,Jm,Km+1)
      T_km_loc = Blk % T (Im,Jm,Km-1)
    endif
    
    T_con = Tg(1)

    call co_Kappa( Blk % matID(Im,Jm,Km), T_loc, kappa_loc )
    call co_Kappa( matIDg, T_con, kappa_con )
          
    ! Tint computations
    select case ( Fm )
      case ( 1 : 2 )
        Den_loc = 2.0d0 * dot_product( M_loc(1,:), Normal ) * modfm2
        Num_loc = 2.0d0*T_loc * dot_product( M_loc(1,:), Normal ) * modfm2 &
                - 0.5d0*(T_jp_loc-T_jm_loc) * dot_product( M_loc(2,:), Normal ) &
                - 0.5d0*(T_kp_loc-T_km_loc) * dot_product( M_loc(3,:), Normal )
        Den_con = 2.0d0 * dot_product( M_con(1,:), Normal ) * modfm2
        Num_con = 2.0d0*T_con * dot_product( M_con(1,:), Normal ) * modfm2 &
                + 0.5d0*(Tg(4)-Tg(3)) * dot_product( M_con(2,:), Normal ) &
                + 0.5d0*(Tg(6)-Tg(5)) * dot_product( M_con(3,:), Normal )
      case ( 3 : 4 )
        Den_loc = 2.0d0 * dot_product( M_loc(2,:), Normal ) * modfm2
        Num_loc = - 0.5d0*(T_ip_loc-T_im_loc) * dot_product( M_loc(1,:), Normal ) &
                + 2.0d0*T_loc * dot_product( M_loc(2,:), Normal ) * modfm2 &
                - 0.5d0*(T_kp_loc-T_km_loc) * dot_product( M_loc(3,:), Normal )
        Den_con = 2.0d0 * dot_product( M_con(2,:), Normal ) * modfm2
        Num_con = 0.5d0*(Tg(4)-Tg(3)) * dot_product( M_con(1,:), Normal ) &
                + 2.0d0*T_con * dot_product( M_con(2,:), Normal ) * modfm2 &
                + 0.5d0*(Tg(6)-Tg(5)) * dot_product( M_con(3,:), Normal )
      case ( 5 : 6 )
        Den_loc = 2.0d0 * dot_product( M_loc(3,:), Normal ) * modfm2
        Num_loc = - 0.5d0*(T_ip_loc-T_im_loc) * dot_product( M_loc(1,:), Normal ) &
                - 0.5d0*(T_jp_loc-T_jm_loc) * dot_product( M_loc(2,:), Normal ) &
                + 2.0d0*T_loc * dot_product( M_loc(3,:), Normal ) * modfm2
        Den_con = 2.0d0 * dot_product( M_con(3,:), Normal ) * modfm2
        Num_con = 0.5d0*(Tg(4)-Tg(3)) * dot_product( M_con(1,:), Normal ) &
                + 0.5d0*(Tg(6)-Tg(5)) * dot_product( M_con(2,:), Normal ) &
                + 2.0d0*T_con * dot_product( M_con(3,:), Normal ) * modfm2
    end select

    Tint = (kappa_loc*Num_loc + kappa_con*Num_con) / (kappa_loc*Den_loc + kappa_con*Den_con)
    
    ! Gradient computations
    select case ( Fm )
    case ( 1 : 2 )
      Gradient (1) = 2.0d0 * ( Tint - T_loc ) * modfm2
      Gradient (2) = ( T_jp_loc - T_jm_loc ) * 0.5d0
      Gradient (3) = ( T_kp_loc - T_km_loc ) * 0.5d0
    case ( 3 : 4 )
      Gradient (1) = ( T_ip_loc - T_im_loc ) * 0.5d0
      Gradient (2) = 2.0d0 * ( Tint - T_loc ) * modfm2
      Gradient (3) = ( T_kp_loc - T_km_loc ) * 0.5d0
    case ( 5 : 6 )
      Gradient (1) = ( T_ip_loc - T_im_loc ) * 0.5d0
      Gradient (2) = ( T_jp_loc - T_jm_loc ) * 0.5d0
      Gradient (3) = 2.0d0 * ( Tint - T_loc ) * modfm2
    end select

    if (ndir == 2) Gradient(3) = 0.0
    Gradient = matmul ( Gradient, M_loc )

    call Compute_Diffusive_Flux ( kappa_loc, Gradient, Area, Normal, Flux )

    ! Residual update
    Blk % r(Im,Jm,Km) = Blk % r(Im,Jm,Km) - modfm2 * Flux
    
  end subroutine BC_Connection

end module FUSS_Lib_BC_Fluxes_Connection