module FUSS_Lib_Solid
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Config_Types_m, only: obj_table

  implicit none
  private
  public :: co_Kappa, co_Rho, co_Cp, co_H, co_T

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


  !> T(h) from the Energy column; outside the table the end interval is extrapolated (then rejected)
  subroutine co_T( matID, h, T )
    implicit none
    real(R8), intent(in)  :: h, matID
    real(R8), intent(out) :: T
    ! Local
    integer   :: Tm, Tp, i, lo, hi
    real(R8)  :: hm, hp

    ! First entry above h; none: last interval; h below the table: first interval.
    ! On a non-decreasing column the entries above h are a tail [Tp, Tmax], so bisecting [lo, hi)
    ! finds the entry of the scan (hi = Tmax+1: none above h, also for a NaN h).
    if ( obj_table%h_monotone(matID) ) then
      lo = obj_table%Tmin
      hi = obj_table%Tmax + 1
      do while ( lo < hi )
        i = lo + (hi - lo) / 2
        if (obj_table%h(matID,i) > h) then
          hi = i
        else
          lo = i + 1
        endif
      enddo
      Tp = min ( lo, obj_table%Tmax )
    else
      Tp = obj_table%Tmax
      do i = obj_table%Tmin, obj_table%Tmax
        if (obj_table%h(matID,i) > h) then
          Tp = i
          exit
        endif
      enddo
    endif
    Tp = max ( Tp, obj_table%Tmin + 1 )
    Tm = Tp - 1
    hm = obj_table%h(matID, Tm)
    hp = obj_table%h(matID, Tp)
    T = Tm + (h-hm)/(hp-hm)*(Tp-Tm)

  end subroutine co_T

end module FUSS_Lib_Solid