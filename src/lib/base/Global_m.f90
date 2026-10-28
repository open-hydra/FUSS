module FUSS_Global_m
  use FUSS_Parameters_m

  implicit none

  character(len=clen) :: FUSS_phase_prefix = ''

  integer :: ndir        ! Number of dimensions of computational frame
  integer :: gc=2        ! ghost cells

contains

  !> The string without double quotes: the variable names are stored quoted (the Tecplot
  !> convention), but quotes inside an XML attribute make the VTK files invalid
  pure function unquoted ( s ) result ( out )
    implicit none
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out
    integer :: i

    out = s
    do i = 1, len(out)
      if (out(i:i) == '"') out(i:i) = ' '
    enddo

  end function unquoted

end module FUSS_Global_m