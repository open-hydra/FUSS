module FUSS_Lib_RK
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: RK_stage, RK_stage_time

  real(R8), parameter, private :: RKcoeff(3,3) = reshape([  1._R8,   0._R8,      0._R8,  &
                                                            1._R8, 0.50_R8,      0._R8,  &
                                                            1._R8, 0.25_R8,2._R8/3._R8 ],&
                                                            shape(RKcoeff), order=[2,1])

contains

  pure function RK_stage ( irk, n_rk, state, oldstate, residual ) result ( newstate )
    implicit none
    integer, intent(in) :: irk
    integer, intent(in) :: n_rk
    real(R8), intent(in) :: state, oldstate, residual
    ! Local
    real(R8)  :: newstate, residual_

    residual_ = ( state - oldstate + residual ) * rkcoeff(n_rk, irk)
    newstate = oldstate + residual_
  
  end function RK_stage

  !> The time, as a fraction of the step, of the state stage irk evaluates its
  !> residual at. Stage k sets u(k) = u(n) + a_k (u(k-1) - u(n) + dt L(u(k-1))), which
  !> puts u(k) at a_k (t(k-1) + 1) from u(0) at 0: Euler 0; RK2 0, 1; RK3 0, 1, 1/2.
  pure function RK_stage_time ( irk, n_rk ) result ( c )
    implicit none
    integer, intent(in) :: irk
    integer, intent(in) :: n_rk
    ! Local
    real(R8) :: c
    integer  :: k

    c = 0._R8
    do k = 1, irk - 1
      c = rkcoeff(n_rk, k) * ( c + 1._R8 )
    enddo

  end function RK_stage_time

end module FUSS_Lib_RK