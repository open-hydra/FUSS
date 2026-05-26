module FUSS_Parameters_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

  integer, parameter  :: hlen=1024
  integer, parameter  :: llen=256
  integer, parameter  :: clen=16
  character(len=clen) :: codename = 'FUSS'
  
  real(R8), parameter :: sigma_SB = 5.67d-8
  integer, dimension(6,3), parameter ::  guide  = reshape( [ 1, 0, 0, &
                                                            -1, 0, 0, &
                                                             0, 1, 0, &
                                                             0,-1, 0, &
                                                             0, 0, 1, &
                                                             0, 0,-1  ] , shape(guide), order=[2,1] )

end module FUSS_Parameters_m