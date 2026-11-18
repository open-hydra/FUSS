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
    use FUSS_Lib_BC_Fluxes_Wall_Ablating, only: BC_Wall_Ablating
    use FUSS_Config_Types_m, only: obj_darcy, obj_ablation, obj_pyrolysis
    use FUSS_Mod_Darcy, only: Darcy_First_Step_Pending
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: f, lower, upper, i, ii, b
    integer  :: Bm, Im, Jm, Km, Fm
    real(R8) :: Ovar(2)    ! wall temperature and flux of an ablating face
    logical  :: count_clamps
    
    ! BC fluxes are computed block by block and, within a block, face by face, over the
    ! records of the blocks this rank owns only (lists built once by build_local_bc_index):
    ! no rank scans the records of the whole domain. The records of one face of one block
    ! write distinct target cells, so one worksharing loop over them is race-free; its
    ! implicit barrier keeps the face contributions of a corner cell summed in the order 1..6.
    ! A wall BC (301-307) also leaves its wall temperature in bc%Tw_face, which the
    ! type-2 gas-enthalpy flux of an open face reads later in the same stage.
    ! The 301/302 time tables are read in the ghost fill (Ghost_Wall), at the
    ! clock Mod_Explicit sets to the stage's time.
    ! The ablating wall counts its B' table clamps only once per stage: not in
    ! the first pass of a new type-2 run's first step, which is discarded.
    count_clamps = .true.
    if ( obj_darcy%enabled ) count_clamps = .not. Darcy_First_Step_Pending()
    blocks: do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle blocks
      faces: do f = 1, 6

        lower = domain % local_bf_first(b,f)
        upper = lower + domain % local_bf_count(b,f) - 1

        !$omp do schedule ( dynamic ) private(ii, i, Bm, Im, Jm, Km, Fm, Ovar)
        do ii = lower, upper
          i  = domain % local_bf_idx(ii)
          Bm = domain % bc(i) % b
          Im = domain % bc(i) % i 
          Jm = domain % bc(i) % j 
          Km = domain % bc(i) % k 
          Fm = domain % bc(i) % f
          select case ( domain % bc(i) % type )
          
            case (101,102) ! connection & chimera (101=block connect, 102=chimera)
              if ( obj_pyrolysis%enabled ) then
                ! k_s(T, chi) on both sides (plan 08 section 2.9; 101 only, type 2: guard 9);
                ! q_if is the TEST ONLY dump of interface-flux.txt (plan 08 R5)
                call BC_Connection ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Mg(1), &
                                     domain % bc(i) % Tg , domain % bc(i) % matIDg, domain % bc(i) % q_if, &
                                     chig = domain % bc(i) % chig )
              elseif ( obj_ablation%dump_interface_flux ) then
                ! TEST ONLY (plan 08 R5): keep the flux for interface-flux.txt
                call BC_Connection ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Mg(1), &
                                     domain % bc(i) % Tg , domain % bc(i) % matIDg, domain % bc(i) % q_if )
              else
                call BC_Connection ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Mg(1), &
                                     domain % bc(i) % Tg , domain % bc(i) % matIDg )
              endif

            ! 300 (symmetry) and 0 (null): no flux across the face

            case (301:305) ! walls: heat flux (301), temperature (302), convection (303), radiation (304, 305)
              call BC_Wall ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) )

            case (306, 307) ! wall: ablating (plan 08), surface energy balance with the B' table and the env
                            ! table (307: the env of the hydra exchange, lambda = 0)
              call BC_Wall_Ablating ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i), Ovar, count_clamps )
              domain % bc(i) % Tw_face = Ovar(1)
            
            case(103,104) ! Multi-Solver Coupling (104 = with chimera)
              domain % blk(Bm) % R(Im,Jm,Km) = domain % blk(Bm) % R(Im,Jm,Km) + domain % bc(i) % ext_flux
          
          end select

        enddo

      enddo faces
    enddo blocks

  end subroutine BC_Fluxes

end module FUSS_Mod_BC_Fluxes
