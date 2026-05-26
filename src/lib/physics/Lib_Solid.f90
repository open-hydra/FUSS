module FUSS_Lib_Solid
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Config_Types_m, only: obj_table

  implicit none
  private
  public :: co_Kappa, co_Rho, co_Cp, co_H, co_T, co_DrhoDT, co_DcpDT

contains

  subroutine co_Kappa( matID, T, kappa )
    implicit none
    real(R8), intent(in)  :: T, matID
    real(R8), intent(out) :: kappa
    ! Local
    integer  :: Tm, Tp
    real(R8) :: km, kp
    
    Tm = int(T)
    Tp = Tm + 1
    km = obj_table%kappa(matID, Tm)
    kp = obj_table%kappa(matID, Tp)
    kappa = km + (T-Tm)/(Tp-Tm)*(kp-km)
  
  end subroutine co_Kappa


  subroutine co_Rho( matID, T, rho )
    implicit none
    real(R8), intent(in)  :: T, matID
    real(R8), intent(out) :: rho
    ! Local
    integer  :: Tm, Tp
    real(R8) :: rhom, rhop
    
    Tm = int(T)
    Tp = Tm + 1
    rhom = obj_table%rho(matID, Tm)
    rhop = obj_table%rho(matID, Tp)
    rho = rhom + (T-Tm)/(Tp-Tm)*(rhop-rhom)
  
  end subroutine co_Rho


  subroutine co_Cp( matID, T, cp )
    implicit none
    real(R8), intent(in)  :: T, matID
    real(R8), intent(out) :: cp
    ! Local
    integer  :: Tm, Tp
    real(R8) :: cpm, cpp
    
    Tm = int(T)
    Tp = Tm + 1
    cpm = obj_table%cp(matID, Tm)
    cpp = obj_table%cp(matID, Tp)
    cp = cpm + (T-Tm)/(Tp-Tm)*(cpp-cpm)
  
  end subroutine co_Cp


  subroutine co_H( matID, T, h )
    implicit none
    real(R8), intent(in)  :: T, matID
    real(R8), intent(out) :: h
    ! Local
    integer  :: Tm, Tp
    real(R8) :: hm, hp
    
    Tm = int(T)
    Tp = Tm + 1
    hm = obj_table%h(matID, Tm)
    hp = obj_table%h(matID, Tp)
    h = hm + (T-Tm)/(Tp-Tm)*(hp-hm)
  
  end subroutine co_H


  subroutine co_T( matID, h, T )
    implicit none
    real(R8), intent(in)  :: h, matID
    real(R8), intent(out) :: T
    ! Local
    integer   :: Tm, Tp, i
    real(R8)  :: hm, hp
    
    do i = 1,size(obj_table%h(matID,:)) 
      if (obj_table%h(matID,i) > h) then
        Tm = i - 1
        Tp = i
        exit
      endif
    enddo
    hm = obj_table%h(matID, Tm)
    hp = obj_table%h(matID, Tp)
    T = Tm + (h-hm)/(hp-hm)*(Tp-Tm)
  
  end subroutine co_T


  subroutine co_DrhoDT( matID, T, drhodT )
    implicit none
    real(R8), intent(in)  :: T, matID
    real(R8), intent(out) :: drhodT
    ! Local
    integer  :: Tm, Tp
    real(R8) :: rhom, rhop
    
    Tm = int(T)
    Tp = Tm + 1
    rhom = obj_table%rho(matID, Tm)
    rhop = obj_table%rho(matID, Tp)
    drhodT = (rhop-rhom)/(Tp-Tm)
  
  end subroutine co_DrhoDT


  subroutine co_DcpDT( matID, T, dcpdT )
    implicit none
    real(R8), intent(in)  :: T, matID
    real(R8), intent(out) :: dcpdT
    ! Local
    integer  :: Tm, Tp
    real(R8) :: cpm, cpp
    
    Tm = int(T)
    Tp = Tm + 1
    cpm = obj_table%cp(matID, Tm)
    cpp = obj_table%cp(matID, Tp)
    dcpdT = (cpp - cpm)/(Tp-Tm)

  end subroutine co_DcpDT

end module FUSS_Lib_Solid