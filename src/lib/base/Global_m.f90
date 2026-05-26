module FUSS_Global_m
  use FUSS_Parameters_m

  implicit none

  character(len=clen) :: FUSS_phase_prefix = ''

  integer :: ndir        ! Number of dimensions of computational frame
  integer :: gc=2        ! ghost cells
  
end module FUSS_Global_m