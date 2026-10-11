module FUSS_Lib_BC_Fluxes_Wall
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  private
  public :: BC_Wall, Wall_Temperature

contains

  !> Wall BCs 301-305: heat flux into the solid through the boundary face.
  subroutine BC_Wall ( Im, Jm, Km, Fm, Blk, bc )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout) :: Blk
    type(FUSS_bc_type), intent(in) :: bc
    ! Local
    integer  :: Dir, Face_i, Face_j, Face_k
    real(R8) :: Tw, q_in

    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    call Wall_Temperature ( Im, Jm, Km, Fm, Blk, bc, .true., Tw, q_in )
    Blk % r (Im,Jm,Km) = Blk % r (Im,Jm,Km) - Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a * q_in

  end subroutine BC_Wall


  !> Wall temperature Tw and heat flux into the solid q_in of a wall BC (301-305).
  !> The face gradient takes its normal part from Tw - T across the half cell and its tangential
  !> parts from the cells along the face, so q_in = G (Tw - T) + C with G = 2 kappa grad(xi_n).n.
  !> 302 prescribes Tw; 301, 303, 304, 305 prescribe q_in = q + hconv (Tref - Tw) + eps sigma (Tref^4 - Tw^4),
  !> with q = qw (301) or qrad (303). Tangential differences use the ghost cells (ghosts = .true.) or
  !> stay inside the block, one-sided at the face edges (ghosts = .false., used to fill the ghosts).
  subroutine Wall_Temperature ( Im, Jm, Km, Fm, Blk, bc, ghosts, Tw, q_in )
    use FUSS_Lib_Solid, only: co_Kappa
    use FUSS_Parameters_m, only: sigma_SB
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(in) :: Blk
    type(FUSS_bc_type), intent(in) :: bc
    logical, intent(in) :: ghosts
    real(R8), intent(out) :: Tw, q_in
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, d, it, lo(3), hi(3)
    real(R8) :: Normal(3), M(3,3), T, kappa, G, C, q, h, es, Tref, f, df, dT

    call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    M = Blk % M (Im,Jm,Km) % c
    T = Blk % T(Im,Jm,Km)
    call co_Kappa( Blk % matID(Im,Jm,Km), T, kappa )

    G = 2d0 * kappa * dot_product( M(Dir,:), Normal )
    C = 0d0
    do d = 1, ndir
      if (d == Dir) cycle
      lo = [Im, Jm, Km]
      hi = lo
      lo(d) = lo(d) - 1
      hi(d) = hi(d) + 1
      if (.not. ghosts) then
        lo(d) = max(lo(d), 1)
        hi(d) = min(hi(d), Blk % dim(d))
        if (hi(d) == lo(d)) cycle
      endif
      C = C + ( Blk % T(hi(1),hi(2),hi(3)) - Blk % T(lo(1),lo(2),lo(3)) ) / (hi(d) - lo(d)) &
            * dot_product( M(d,:), Normal )
    enddo
    C = modfm2 * kappa * C

    select case ( bc % type )
      case (302)
        Tw = bc % Tw
        q_in = G * (Tw - T) + C
        return
      case (301)
        q = bc % qw;   h = 0d0;        es = 0d0;                 Tref = 0d0
      case (303)
        q = bc % qrad; h = bc % hconv; es = 0d0;                 Tref = bc % Tref
      case (304)
        q = 0d0;       h = 0d0;        es = sigma_SB * bc % eps; Tref = bc % Tref
      case (305)
        q = 0d0;       h = bc % hconv; es = sigma_SB * bc % eps; Tref = bc % Tref
      case default
        Tw = T
        q_in = 0d0
        return
    end select

    ! Newton on q + h (Tref - Tw) + es (Tref^4 - Tw^4) = G (Tw - T) + C: decreasing and concave
    ! in Tw > 0, so from the f < 0 side Newton converges monotonically; a step below zero is
    ! replaced by halving (the model has no admissible root). Exact in one step without
    ! radiation; no conductance at all (k = 0, 301): Tw = T.
    Tw = max ( T, 1d0 )
    do it = 1, 50
      f  = q + h*(Tref - Tw) + es*(Tref**4 - Tw**4) - G*(Tw - T) - C
      df = h + 4d0*es*Tw**3 + G
      if (df <= 0d0) exit
      dT = f / df
      if (Tw + dT > 0d0) then
        Tw = Tw + dT
      else
        Tw = 0.5d0*Tw
        cycle
      endif
      if (es == 0d0 .or. abs(dT) <= 1d-12*abs(Tw)) exit
    enddo
    q_in = q + h*(Tref - Tw) + es*(Tref**4 - Tw**4)

  end subroutine Wall_Temperature

end module FUSS_Lib_BC_Fluxes_Wall
