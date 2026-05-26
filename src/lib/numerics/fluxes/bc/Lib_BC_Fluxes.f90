module FUSS_Lib_BC_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  public

contains

  subroutine Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    implicit none
    integer, intent(in) :: fm
    integer, intent(out) :: modfm, modfm1, modfm2, modfm3

    modFm  = mod(Fm,2)
    modfm1 = 1-modFm  
    modfm2 = 1-2*modFm
    modfm3 = 2*( 2*modFm-1 )
    
  end subroutine Compute_Modfm


  subroutine Face_Index ( face, dir, cell_i, cell_j, cell_k, face_i, face_j, face_k )
    implicit none
    integer, intent(in) :: face, cell_i, cell_j, cell_k
    integer, intent(out) :: dir, face_i, face_j, face_k

    select case ( face )
      case(1:2)
        dir = 1
        face_i = cell_i - mod ( face, 2 )
        face_j = cell_j
        face_k = cell_k
      case(3:4)
        dir = 2
        face_i = cell_i
        face_j = cell_j - mod ( face, 2 )
        face_k = cell_k
      case(5:6)
        dir = 3
        face_i = cell_i
        face_j = cell_j
        face_k = cell_k - mod ( face, 2 )
    end select

  end subroutine Face_Index


  subroutine Compute_Wall_Properties(Tw, qw, exit_array)
    implicit none
    real(R8), dimension(2), intent(inout) :: exit_array
    real(R8),intent(in) :: Tw, qw

    ! Wall temperature
    exit_array(1) = Tw
    ! Heat flux
    exit_array(2) = qw

  end subroutine Compute_Wall_Properties

end module FUSS_Lib_BC_Fluxes