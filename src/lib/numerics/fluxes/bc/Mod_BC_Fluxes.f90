module FUSS_Mod_BC_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: BC_Fluxes

contains

  subroutine BC_Fluxes ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Mod_MPI, only: is_local_block
    use FUSS_Lib_BC_Fluxes_Connection
    use FUSS_Lib_BC_Fluxes_Wall
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: f, lower, upper, i, b
    integer  :: Bm, Im, Jm, Km, Fm
    
    ! BC fluxes are computed in order of block and face type.
    upper = 0

    blocks: do b = 1, domain % nb
      faces: do f = 1, 6

        lower = upper + 1                   ! Update lower bound
        upper = upper + domain % n_bf(b,f)  ! Upper bound: add number of cells on face f of block b
        
        !$omp do schedule ( dynamic ) private(i, Bm, Im, Jm, Km, Fm)
        do i = lower, upper
          Bm = domain % bc(i) % b
          if (.not. is_local_block(Bm)) cycle
          Im = domain % bc(i) % i 
          Jm = domain % bc(i) % j 
          Km = domain % bc(i) % k 
          Fm = domain % bc(i) % f
          select case ( domain % bc(i) % type )
          
            case (101,102) ! connection & chimera (101=block connect, 102=chimera)
              call BC_Connection ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Mg(1), &
                                   domain % bc(i) % Tg , domain % bc(i) % matIDg )

            ! 300 (symmetry) and 0 (null): no flux across the face

            case (301:305) ! walls: heat flux (301), temperature (302), convection (303), radiation (304, 305)
              call BC_Wall ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) )

            case(103,104) ! Multi-Solver Coupling (104 = with chimera)
              domain % blk(Bm) % R(Im,Jm,Km) = domain % blk(Bm) % R(Im,Jm,Km) + domain % bc(i) % ext_flux
          
          end select

        enddo

      enddo faces
    enddo blocks

  end subroutine BC_Fluxes

end module FUSS_Mod_BC_Fluxes