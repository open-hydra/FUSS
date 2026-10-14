module FUSS_Lib_Multigrid
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Check_Multigrid, Coarse_Grid, Coarse_IOfield, Coarse_matID
  public :: fine2coarse_prim, coarse2fine_prim
  public :: fas_store_h0, fas_build_forcing, coarse2fine_corr

  integer, public :: fas_fallbacks = 0   ! corrections halved or dropped by the admissibility guard

contains

subroutine Check_Multigrid ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_multigrid
    use ir_precision, only: str
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, d, rap, check

    rap = 2**(obj_multigrid%MGL-1)
    do b = 1, domain % nb
      check = 0
      do d = 1, 3
        if ( Mod ( domain % Blk(b) % Dim(d), rap ) == 0 ) then
          check = check + 1
        elseif ( d == 3 .and. domain % Blk(b) % Dim(d) == 1 ) then
          check = check + 1
        endif
      enddo
      if ( check < 3 ) then
        write(*,'(A)') ' [ERROR] multigrid with '//trim(str(.true.,obj_multigrid%MGL))//' levels needs block cells '// &
                       'divisible by '//trim(str(.true.,rap))//'; block '//trim(str(.true.,b))//' is not'
        error stop 'multigrid: block cells not divisible'
      endif
      ! The coarsest level needs at least 2 cells per coarsened direction
      do d = 1, 3
        if ( domain % Blk(b) % Dim(d) > 1 .and. domain % Blk(b) % Dim(d) / rap < 2 ) then
          write(*,'(A)') ' [ERROR] multigrid with '//trim(str(.true.,obj_multigrid%MGL))//' levels leaves block '// &
                         trim(str(.true.,b))//' with '//trim(str(.true.,domain % Blk(b) % Dim(d) / rap))// &
                         ' cell in direction '//trim(str(.true.,d))//' on the coarsest level; at least 2 are needed'
          error stop 'multigrid: coarsest level thinner than 2 cells'
        endif
      enddo
    enddo

  end subroutine Check_Multigrid


  subroutine fine2coarse_prim ( fPrim, cPrim, fmatID, cmatID, fqvol, cqvol, fVol, cVol, fDim, cDim )
    use FUSS_Global_m
    use FUSS_Lib_Solid
    implicit none
    integer, intent(in), dimension(3) :: fDim, cDim
    real(R8), intent(in),  dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fmatID
    real(R8), intent(out), dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cmatID
    real(R8), intent(in),  dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fqvol
    real(R8), intent(out), dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cqvol
    real(R8), intent(in),  dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fPrim
    real(R8), intent(out), dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cPrim
    real(R8), intent(in),  dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fVol
    real(R8), intent(in),  dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cVol
    ! Local
    integer :: i, j, k, i2, j2, k2, i2d, j2d, k2d
    real(R8), dimension(1:fDim(1),1:fDim(2),1:fDim(3)) :: fCons
    real(R8), dimension(1:cDim(1),1:cDim(2),1:cDim(3)) :: cCons
    real(R8) :: tot_vol

    !$omp do collapse (2)
    do k = 1, fDim(3)
    do j = 1, fDim(2)
    do i = 1, fDim(1)
      call co_H( fmatID(i,j,k), fPrim(i,j,k), fCons(i,j,k) )
    end do ; end do ; end do

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)

      i2 = 2*i
      j2 = 2*j
      k2 = 2*k
      i2d = i2-1
      j2d = j2-1
      k2d = k2-1

      if ( fDim(3) == 1 ) then ! 2D
        cCons(i,j,k) = fCons(i2d,j2d,1) * fVol(i2d,j2d,1)  &
                     + fCons(i2, j2d,1) * fVol(i2, j2d,1)  &
                     + fCons(i2d,j2, 1) * fVol(i2d,j2, 1)  &
                     + fCons(i2, j2, 1) * fVol(i2, j2, 1)
        
        cqvol(i,j,k) = fqvol(i2d,j2d,1) * fVol(i2d,j2d,1)  &
                     + fqvol(i2, j2d,1) * fVol(i2, j2d,1)  &
                     + fqvol(i2d,j2, 1) * fVol(i2d,j2, 1)  &
                     + fqvol(i2, j2, 1) * fVol(i2, j2, 1)
        tot_vol = fVol(i2d,j2d,1) + fVol(i2, j2d,1) + fVol(i2d,j2, 1) + fVol(i2, j2, 1)
        cqvol(i,j,k) = cqvol(i,j,k) / tot_vol

      else ! 3D
        cCons(i,j,k) = fCons(i2d,j2d,k2d) * fVol(i2d,j2d,k2d)  &
                     + fCons(i2, j2d,k2d) * fVol(i2, j2d,k2d)  &
                     + fCons(i2d,j2, k2d) * fVol(i2d,j2, k2d)  &
                     + fCons(i2, j2, k2d) * fVol(i2, j2, k2d)  &
                     + fCons(i2d,j2d,k2 ) * fVol(i2d,j2d,k2 )  &
                     + fCons(i2, j2d,k2 ) * fVol(i2, j2d,k2 )  &
                     + fCons(i2d,j2, k2 ) * fVol(i2d,j2, k2 )  &
                     + fCons(i2, j2, k2 ) * fVol(i2, j2, k2 )
        
        cqvol(i,j,k) = fqvol(i2d,j2d,k2d) * fVol(i2d,j2d,k2d)  &
                     + fqvol(i2, j2d,k2d) * fVol(i2, j2d,k2d)  &
                     + fqvol(i2d,j2, k2d) * fVol(i2d,j2, k2d)  &
                     + fqvol(i2, j2, k2d) * fVol(i2, j2, k2d)  &
                     + fqvol(i2d,j2d,k2 ) * fVol(i2d,j2d,k2 )  &
                     + fqvol(i2, j2d,k2 ) * fVol(i2, j2d,k2 )  &
                     + fqvol(i2d,j2, k2 ) * fVol(i2d,j2, k2 )  &
                     + fqvol(i2, j2, k2 ) * fVol(i2, j2, k2 )
        tot_vol = fVol(i2d,j2d,k2d) + fVol(i2, j2d,k2d) + fVol(i2d,j2, k2d) + fVol(i2, j2, k2d) &
                + fVol(i2d,j2d,k2 ) + fVol(i2, j2d,k2 ) + fVol(i2d,j2, k2 ) + fVol(i2, j2, k2 )
        cqvol(i,j,k) = cqvol(i,j,k) / tot_vol

      endif
      
      cmatID(i,j,k) = fmatID(i2,j2,min(k2,fDim(3)))  ! 2D: k2 would be the ghost layer

    enddo; enddo; enddo

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)
      cCons(i,j,k) = cCons(i,j,k) / cVol(i,j,k)
      call co_T( cmatID(i,j,k), cCons(i,j,k), cPrim(i,j,k) )
    enddo; enddo; enddo

  end subroutine fine2coarse_prim


  subroutine coarse2fine_prim ( fPrim, cPrim, fDim, cDim )
    use FUSS_Global_m
    implicit none
    integer, intent(in), dimension(3) :: fDim, cDim
    real(R8), intent(out), dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fPrim
    real(R8), intent(in),  dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cPrim

    call prolong_trilinear ( cPrim, fPrim, cDim, fDim, gc, gc )

  end subroutine coarse2fine_prim


  !> Trilinear prolongation of one coarse-cell field onto the fine interior cells, stencils
  !> clamped at the block boundary; ch and fh are the halo widths of the two arrays. The
  !> worksharing directive is inert when the caller is outside a parallel region.
  subroutine prolong_trilinear ( c, f, cDim, fDim, ch, fh )
    implicit none
    integer, intent(in), dimension(3) :: cDim, fDim
    integer, intent(in) :: ch, fh
    real(R8), intent(in),  dimension(1-ch:cDim(1)+ch,1-ch:cDim(2)+ch,1-ch:cDim(3)+ch) :: c
    real(R8), intent(out), dimension(1-fh:fDim(1)+fh,1-fh:fDim(2)+fh,1-fh:fDim(3)+fh) :: f
    ! Local
    integer :: i, j, k, rap, i2, j2, k2, i2d, j2d, k2d, im, jm, km, ip, jp, kp
    integer :: ii, jj, kk, counter, mask(3), id(6)
    real(R8) :: a1, a2, a3, a4, coeffs(8), interp

    a1 = 27d0/64d0
    a2 = 9d0/64d0
    a3 = 3d0/64d0
    a4 = 1d0/64d0
    coeffs(1:8) = [ a1, a2, a2, a2, a3, a3, a3, a4 ]
    rap = 2
    
    if ( fDim(3) == 1 ) then ! 2D
      a1 = 9.d0/16.d0
      a2 = 3.d0/16.d0
      a3 = 1.d0/16.d0
      a4 = 0d0
      coeffs(1:8) = [ a1, a2, a2, a4, a3, a4, a4, a4 ]
    end if

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)            
    do i = 1, cDim(1)
            
      i2 = rap*i
      j2 = rap*j
      k2 = rap*k
      i2d = i2-(rap-1)
      j2d = j2-(rap-1)
      k2d = k2-(rap-1)
          
      im = Max (1,i-1)
      jm = Max (1,j-1)
      km = Max (1,k-1)

      ip = Min (cDim(1), i+1)
      jp = Min (cDim(2), j+1)
      kp = Min (cDim(3), k+1)

      if ( fDim(3) == 1 ) then
        k2 = 1
        k2d = 1 
      end if ! 2D

      id(1:6) = [ im, jm, km, ip, jp, kp ]
          
      counter = 1
          
      do kk = k2d, k2          
      do jj = j2d, j2            
      do ii = i2d, i2

        if ( fDim(3) == 1 ) then ! 2D
          if (counter==1) mask(1:3) = [1,2,3]
          if (counter==2) mask(1:3) = [4,2,3]
          if (counter==3) mask(1:3) = [1,5,3]
          if (counter==4) mask(1:3) = [4,5,3]
        else
          if (counter==1) mask(1:3) = [1,2,3]
          if (counter==2) mask(1:3) = [4,2,3]
          if (counter==3) mask(1:3) = [1,5,3]
          if (counter==4) mask(1:3) = [4,5,3]
          if (counter==5) mask(1:3) = [1,2,6]
          if (counter==6) mask(1:3) = [4,2,6]
          if (counter==7) mask(1:3) = [1,5,6]
          if (counter==8) mask(1:3) = [4,5,6]
        endif  

        interp = coeffs(1) * c(i,j,k) &
               + coeffs(2) * c(id(mask(1)),j,k) &
               + coeffs(3) * c(i,id(mask(2)),k)    &
               + coeffs(4) * c(i,j,id(mask(3)))    &
               + coeffs(5) * c(id(mask(1)),id(mask(2)),k)   &
               + coeffs(6) * c(id(mask(1)),j,id(mask(3)))   &
               + coeffs(7) * c(i,id(mask(2)),id(mask(3)))   &
               + coeffs(8) * c(id(mask(1)),id(mask(2)),id(mask(3)))

        f(ii,jj,kk) = interp

        counter = counter + 1

      enddo; enddo; enddo

    enddo; enddo; enddo

  end subroutine prolong_trilinear


  subroutine fas_store_h0 ( cT, cmatID, h0, cDim )
    use FUSS_Global_m
    use FUSS_Lib_Solid, only: co_H
    implicit none
    integer, intent(in), dimension(3) :: cDim
    real(R8), intent(in),  dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cT, cmatID
    real(R8), intent(out), dimension(1:cDim(1),1:cDim(2),1:cDim(3)) :: h0
    ! Local
    integer :: i, j, k

    do k = 1, cDim(3); do j = 1, cDim(2); do i = 1, cDim(1)
      call co_H ( cmatID(i,j,k), cT(i,j,k), h0(i,j,k) )
    enddo; enddo; enddo

  end subroutine fas_store_h0


  !> FAS forcing of a coarse block: the volume-integrated steady residual of the fine children,
  !> summed (conservative), minus the coarse one at the restricted state. The residual convention is
  !> N = -R + qvol*V (see Scale_Residual_Cons); R of a level under forcing already carries -fas_f,
  !> so the sum telescopes the forcing of the finer level (Jameson's Q = sum, P_k collected).
  subroutine fas_build_forcing ( fR, fqvol, fVol, cR, cqvol, cVol, fas_f, fDim, cDim )
    use FUSS_Global_m
    implicit none
    integer, intent(in), dimension(3) :: fDim, cDim
    real(R8), intent(in), dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fR, fqvol, fVol
    real(R8), intent(in), dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cR, cqvol, cVol
    real(R8), intent(out), dimension(1:cDim(1),1:cDim(2),1:cDim(3)) :: fas_f
    ! Local
    integer :: i, j, k, i2, j2, k2, i2d, j2d, k2d
    real(R8) :: s

    do k = 1, cDim(3); do j = 1, cDim(2); do i = 1, cDim(1)
      i2 = 2*i; j2 = 2*j; k2 = 2*k
      i2d = i2-1; j2d = j2-1; k2d = k2-1
      if ( fDim(3) == 1 ) then ! 2D
        s = - fR(i2d,j2d,1) + fqvol(i2d,j2d,1)*fVol(i2d,j2d,1) &
            - fR(i2, j2d,1) + fqvol(i2, j2d,1)*fVol(i2, j2d,1) &
            - fR(i2d,j2, 1) + fqvol(i2d,j2, 1)*fVol(i2d,j2, 1) &
            - fR(i2, j2, 1) + fqvol(i2, j2, 1)*fVol(i2, j2, 1)
      else
        s = - fR(i2d,j2d,k2d) + fqvol(i2d,j2d,k2d)*fVol(i2d,j2d,k2d) &
            - fR(i2, j2d,k2d) + fqvol(i2, j2d,k2d)*fVol(i2, j2d,k2d) &
            - fR(i2d,j2, k2d) + fqvol(i2d,j2, k2d)*fVol(i2d,j2, k2d) &
            - fR(i2, j2, k2d) + fqvol(i2, j2, k2d)*fVol(i2, j2, k2d) &
            - fR(i2d,j2d,k2 ) + fqvol(i2d,j2d,k2 )*fVol(i2d,j2d,k2 ) &
            - fR(i2, j2d,k2 ) + fqvol(i2, j2d,k2 )*fVol(i2, j2d,k2 ) &
            - fR(i2d,j2, k2 ) + fqvol(i2d,j2, k2 )*fVol(i2d,j2, k2 ) &
            - fR(i2, j2, k2 ) + fqvol(i2, j2, k2 )*fVol(i2, j2, k2 )
      endif
      fas_f(i,j,k) = s - ( - cR(i,j,k) + cqvol(i,j,k)*cVol(i,j,k) )
    enddo; enddo; enddo

  end subroutine fas_build_forcing


  !> FAS: prolong the coarse correction in h, delta = h(w_2h) - h0, with the trilinear weights of
  !> coarse2fine_prim, and add it to the fine state.
  subroutine coarse2fine_corr ( fT, fmatID, cT, cmatID, h0, fDim, cDim )
    use FUSS_Global_m
    use FUSS_Lib_Solid, only: co_H, co_T
    use FUSS_Lib_Newstate, only: Admissible_T, IRS_HALVINGS
    implicit none
    integer, intent(in), dimension(3) :: fDim, cDim
    real(R8), intent(inout), dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fT
    real(R8), intent(in),    dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fmatID
    real(R8), intent(in),    dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cT, cmatID
    real(R8), intent(in),    dimension(1:cDim(1),1:cDim(2),1:cDim(3)) :: h0
    ! Local
    integer :: i, j, k, n
    real(R8) :: h, Tnew, theta
    real(R8), dimension(1:cDim(1),1:cDim(2),1:cDim(3)) :: dh
    real(R8), dimension(1:fDim(1),1:fDim(2),1:fDim(3)) :: dhf

    ! Correction on the coarse cells, prolongated with the weights of coarse2fine_prim
    do k = 1, cDim(3); do j = 1, cDim(2); do i = 1, cDim(1)
      call co_H ( cmatID(i,j,k), cT(i,j,k), h )
      dh(i,j,k) = h - h0(i,j,k)
    enddo; enddo; enddo
    call prolong_trilinear ( dh, dhf, cDim, fDim, 0, 0 )

    do k = 1, fDim(3)
    do j = 1, fDim(2)
    do i = 1, fDim(1)

      call co_H ( fmatID(i,j,k), fT(i,j,k), h )
      theta = 1d0
      do n = 0, IRS_HALVINGS
        call co_T ( fmatID(i,j,k), h + theta*dhf(i,j,k), Tnew )
        if ( Admissible_T( Tnew ) ) exit
        theta = 0.5d0 * theta
        if (n == IRS_HALVINGS) Tnew = fT(i,j,k)
      enddo
      if (theta < 1d0) fas_fallbacks = fas_fallbacks + 1
      fT(i,j,k) = Tnew

    enddo; enddo; enddo

  end subroutine coarse2fine_corr


  subroutine Coarse_matID ( Fine, Coarse )
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: Fine, Coarse
    ! Local
    integer :: b, i, j, k

    do b = 1, Coarse % nb
      do k = 1, Coarse % Blk(b) % dim(3)
      do j = 1, Coarse % Blk(b) % dim(2)
      do i = 1, Coarse % Blk(b) % dim(1)
        Coarse % Blk(b) % matID(i,j,k) = Fine % Blk(b) % matID(2*i, 2*j, min(2*k, Fine % Blk(b) % dim(3)))
      enddo; enddo; enddo
    enddo

  end subroutine Coarse_matID


  subroutine Coarse_Grid ( Fine, Coarse )
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type), intent(inout) :: Fine, Coarse
    ! Local
    integer :: b

    do b = 1, Coarse % nb ! Loop over blocks
      call Coarse_Grid_Blk ( Fine % Blk(b) % node, Coarse % Blk(b) % node, &
                             Fine % Blk(b) % dim,  Coarse % Blk(b) % dim )
    enddo
  
  contains

    subroutine Coarse_Grid_Blk ( fNode, cNode, fDim, cDim )
      use FUSS_Base_Types_m
      implicit none
      integer, intent(in), dimension(3) :: fDim, cDim
      type(FUSS_vector_3D_type), intent(in),  dimension(0:fDim(1),0:fDim(2),0:fDim(3)) :: fNode
      type(FUSS_vector_3D_type), intent(out), dimension(0:cDim(1),0:cDim(2),0:cDim(3)) :: cNode
      ! Local
      integer :: i, j, k, i2, j2, k2

      !$omp do collapse (2)
      do k = 0, fDim(3), 2-mod(fDim(3),2)
      do j = 0, fDim(2), 2
      do i = 0, fDim(1), 2
        
        i2 = i / 2
        j2 = j / 2
        k2 = k / 2

        if ( fDim(3) == 1 ) then ! 2D
          k2 = k
        end if

        cNode(i2,j2,k2) % c = fNode(i,j,k) % c
      
      enddo; enddo; enddo
    
    end subroutine Coarse_Grid_Blk

  end subroutine Coarse_Grid


  subroutine Coarse_IOfield ( fIOfield, cIOfield )
    use Lib_ORION_data
    implicit none
    type(ORION_data), intent(in)    :: fIOfield
    type(ORION_data), intent(inout) :: cIOfield
    ! Local
    integer :: b, nb, Ni, Nj, Nk, Onvar

    Onvar = 3  ! T, matID, qvol
    
    nb = size( fIOfield % block )
    allocate ( cIOfield % block( nb ) )
    
    do b = 1, nb
      Ni = fIOfield % block(b)%Ni /2
      Nj = fIOfield % block(b)%Nj /2
      Nk = fIOfield % block(b)%Nk /2
      Nk = Max ( 1, Nk ) ! 2D case
      cIOfield % block(b) % Ni = Ni
      cIOfield % block(b) % Nj = Nj
      cIOfield % block(b) % Nk = Nk
      allocate( cIOfield%block(b)%mesh(1:3,0:Ni,0:Nj,0:Nk) )
      allocate( cIOfield%block(b)%vars(1:Onvar,Ni,Nj,Nk) )
    enddo

  end subroutine Coarse_IOfield

end module FUSS_Lib_Multigrid