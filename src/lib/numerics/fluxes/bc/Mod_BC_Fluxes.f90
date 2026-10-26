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
    integer  :: f, lower, upper, i, ii, b
    integer  :: Bm, Im, Jm, Km, Fm
    
    ! BC fluxes are computed block by block and, within a block, face by face, over the
    ! records of the blocks this rank owns only (lists built once by build_local_bc_index):
    ! no rank scans the records of the whole domain. The records of one face of one block
    ! write distinct target cells, so one worksharing loop over them is race-free; its
    ! implicit barrier keeps the face contributions of a corner cell summed in the order 1..6.
    blocks: do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle blocks
      faces: do f = 1, 6

        lower = domain % local_bf_first(b,f)
        upper = lower + domain % local_bf_count(b,f) - 1

        !$omp do schedule ( dynamic ) private(ii, i, Bm, Im, Jm, Km, Fm)
        do ii = lower, upper
          i  = domain % local_bf_idx(ii)
          Bm = domain % bc(i) % b
          Im = domain % bc(i) % i 
          Jm = domain % bc(i) % j 
          Km = domain % bc(i) % k 
          Fm = domain % bc(i) % f
          select case ( domain % bc(i) % type )
          
            case (101) ! block connection: stencil face metric, interior-face conductivity
              call BC_Connection_Stencil ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) )

            case (102) ! chimera: two-sided balance with the cell metric and the ghost metric
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