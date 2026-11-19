module FUSS_Lib_BC_Fluxes_Wall
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Face_Tangential

  implicit none
  private
  public :: BC_Wall, Wall_Temperature, Wall_Conduction, Wall_Newton

  !> 304/305 wall-temperature solves of this rank that reached the 50-step Newton cap (the last iterate is
  !> used, plan 14 D17); summed over the ranks and printed once at the end of the run (Wrap_Postprocess).
  !> A real: exact to 2**53, reduced by mpi_reduce_sum_r8 as it is.
  real(R8), public :: wall_newton_capped = 0d0

contains

  !> Wall BCs 301-305: heat flux into the solid through the boundary face. The wall temperature
  !> stays in bc%Tw_face: the type-2 gas-enthalpy flux of an open face (Mod_Darcy) and the exit
  !> temperature of the type-1 gas (Mod_Pyrolysis) read it later in the same stage.
  subroutine BC_Wall ( Im, Jm, Km, Fm, Blk, bc )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout) :: Blk
    type(FUSS_bc_type), intent(inout) :: bc
    ! Local
    integer  :: Dir, Face_i, Face_j, Face_k
    real(R8) :: Tw, q_in

    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    call Wall_Temperature ( Im, Jm, Km, Fm, Blk, bc, .true., Tw, q_in )
    Blk % r (Im,Jm,Km) = Blk % r (Im,Jm,Km) - Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a * q_in
    bc % Tw_face = Tw

  end subroutine BC_Wall


  !> Wall temperature Tw and heat flux into the solid q_in of a wall BC (301-307).
  !> The face gradient takes its normal part from Tw - T across the half cell and its tangential
  !> parts from the cells along the face, so q_in = G (Tw - T) + C (Wall_Conduction).
  !> 302 prescribes Tw; 301, 303, 304, 305 prescribe q_in = q + hconv (Tref - Tw) + eps sigma (Tref^4 - Tw^4),
  !> with q = qw (301) or qrad (303). 306, 307 (the ablating walls, plan 08) give the Tw of the face's last
  !> converged surface energy balance (bc%Tw_prev; the cell's before the first one): their balance is
  !> FUSS_Lib_BC_Fluxes_Wall_Ablating's, this routine serves their ghost cells. Tangential differences use
  !> the ghost cells (ghosts = .true.) or stay inside the block, one-sided at the face edges
  !> (ghosts = .false., used to fill the ghosts).
  !> The metric M is that of the same differences of the positions (blk%bm(Fm)%wM, Lib_Stencil_Metrics, one variant per
  !> value of ghosts): D = J g for a linear field, so Tw balances the prescribed flux at T(x_f) exactly.
  subroutine Wall_Temperature ( Im, Jm, Km, Fm, Blk, bc, ghosts, Tw, q_in )
    use FUSS_Parameters_m, only: sigma_SB
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(in) :: Blk
    type(FUSS_bc_type), intent(in) :: bc
    logical, intent(in) :: ghosts
    real(R8), intent(out) :: Tw, q_in
    ! Local
    real(R8) :: T, G, C, q, h, es, Tref
    logical  :: converged

    call Wall_Conduction ( Im, Jm, Km, Fm, Blk, ghosts, T, G, C )

    select case ( bc % type )
      case (302)
        Tw = bc % Tw
        q_in = G * (Tw - T) + C
        return
      case (306, 307)
        Tw = bc % Tw_prev
        if ( .not. ( Tw > 0d0 ) ) Tw = T
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

    call Wall_Newton ( T, G, C, q, h, es, Tref, Tw, q_in, converged )

    ! Not converged in the 50 steps: the last iterate is used, as upstream (plan 14 D17, user 2026-11-11;
    ! plan 10 D12's stop retired). Every such solve is counted; the first one of the rank prints one line
    ! (this routine runs in OpenMP threads: the count and the line are one critical section; no MPI here).
    if ( .not. converged ) then
      !$omp critical (wall_newton_cap)
      wall_newton_capped = wall_newton_capped + 1d0
      if ( wall_newton_capped == 1d0 ) &
        write(*,'(A,I0,A,I0,A,3(1X,I0),A,I0,A,ES12.4,A,ES12.4,A)') ' [WARNING] wall BC ', bc % type, &
          ': the wall temperature did not converge in 50 Newton steps (block ', bc % b, ', cell', Im, Jm, Km, &
          ', face ', Fm, ', T_wall ', Tw, ' K, residual ', q_in - G*(Tw - T) - C, &
          ' W/m2): the last iterate is used; the count is printed at the end of the run'
      !$omp end critical (wall_newton_cap)
    endif

  end subroutine Wall_Temperature


  !> The conduction through a wall face, q_in = G (Tw - T) + C: G = 2 kappa grad(xi_n).n over the
  !> half cell, C the tangential part from the cells along the face (ghosts = .true.: across the face
  !> edges into the ghost cells; .false.: inside the block, one-sided at the face edges). T is the
  !> boundary cell's temperature; kappa its k_s(T, chi) under pyrolysis (plan 06 section 3.2), else co_Kappa(T).
  subroutine Wall_Conduction ( Im, Jm, Km, Fm, Blk, ghosts, T, G, C )
    use FUSS_Lib_Solid, only: co_Kappa
    use FUSS_Config_Types_m, only: obj_pyrolysis
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(in) :: Blk
    logical, intent(in) :: ghosts
    real(R8), intent(out) :: T, G, C
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, d, lo(3), hi(3), tf(2)
    real(R8) :: Normal(3), M(3,3), kappa

    call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    M = Blk % M (Im,Jm,Km) % c
    if ( allocated ( Blk % bm(Fm) % wM ) ) then
      tf = Face_Tangential ( Fm, Im, Jm, Km )
      if ( Blk % bm(Fm) % wset(tf(1),tf(2)) ) M = Blk % bm(Fm) % wM(:,:, merge(1, 2, ghosts), tf(1), tf(2))
    end if
    T = Blk % T(Im,Jm,Km)
    if ( obj_pyrolysis % enabled ) then
      kappa = Blk % k_s(Im,Jm,Km)     ! (T, chi) cache, plan 06 section 3.2
    else
      call co_Kappa( Blk % matID(Im,Jm,Km), T, kappa )
    endif

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

  end subroutine Wall_Conduction


  !> Tw and q_in of a face with q_in = q + h (Tref - Tw) + es (Tref^4 - Tw^4) = G (Tw - T) + C.
  !> converged = .false. only when the 50 steps ran out without the step falling below 1e-12 Tw.
  subroutine Wall_Newton ( T, G, C, q, h, es, Tref, Tw, q_in, converged )
    implicit none
    real(R8), intent(in)  :: T, G, C, q, h, es, Tref
    real(R8), intent(out) :: Tw, q_in
    logical,  intent(out) :: converged
    ! Local
    integer  :: it
    real(R8) :: f, df, dT

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
    converged = ( it <= 50 )
    q_in = q + h*(Tref - Tw) + es*(Tref**4 - Tw**4)

  end subroutine Wall_Newton

end module FUSS_Lib_BC_Fluxes_Wall
