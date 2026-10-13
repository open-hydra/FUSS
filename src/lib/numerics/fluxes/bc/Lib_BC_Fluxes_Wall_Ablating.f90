!> The ablating wall, BC 306 (plan 08 section 2.1): the surface energy balance
!> of thesis eqs. 3.46-3.50, PATO's BprimeBoundaryConditions, per boundary face.
!>
!>   F(T_w) = G (T_w - T_c)                                  conduction into the solid
!>          - rhoeUeCH0 beta (h_r - h_w)                     convection, blowing-corrected
!>          - mdot_pg (h_pg,w - h_w) - mdot_ch (h_ch,w - h_w) gas and char advection
!>          - alpha_s q_rad - eps_s sigma T_inf^4 + eps_s sigma T_w^4
!>          = 0
!>
!> G = 2 kappa M(Dir,:).n is the half-cell conductance BC 303 uses. h_w and
!> B'_c come from the B' table at (p_w, B'_pg, T_w) (FUSS_Lib_Bprime, PATO's
!> interpolation); beta = ln(1 + 2 lambda B')/(2 lambda B'), B' = B'_pg + B'_ch,
!> with PATO's exceptions (beta = 1 when lambda = 0, B' = 0, mdot_pg <= 1e-6 or
!> rhoeUeCH0 <= 1e-4; rhoeUeCH0 floored at 1e-6; beta updated only when the
!> chemistry is on). With the chemistry off (env column chem = 0, or the
!> 0 < chem < 1 ramp) the convection is hconv (T_edge - T_w) and the gas
!> advection is dropped, PATO's chemistryOn == 0 balance.
!>
!> Solve: a safeguarded Newton on T_w with a central finite-difference slope
!> (dT = 0.5 K; the tables are piecewise linear) and bisection whenever a step
!> leaves the bracket [max(T_c - 500, T_lo + 1), T_hi - 1] of the material's
!> tables; at every Newton iterate beta is taken to its Picard fixed point
!> (|d beta| <= seb-tol). Converged when, at the current iterate, the step
!> |dT_w| <= seb-tol T_w and |F| <= seb-tol max|term|; that iterate is
!> returned. The cap seb-max-iter (on either loop) is an error stop 3 (an
!> unconverged T_w would feed the recession).
!>
!> Chemistry off with a zero emissivity: the routine hands the face to the BC
!> 303 kernel with (hconv, T_edge, alpha_s q_rad), the same closed form, so a
!> 306 face with such an env table reproduces a 303 face bit for bit.
!>
!> On a static mesh the char is not consumed (PATO zeroes B'_c and mdot_c when
!> the mesh does not move): mdot_ch = B'_ch = sdot = 0, the pyrolysis-gas
!> blowing and h_w from the table stay.
!>
!> The flux into the residual is BC 302's at the converged T_w. The face keeps
!> its last converged state in bc (seb.txt, the wall output, the restart file).
module FUSS_Lib_BC_Fluxes_Wall_Ablating
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Config_Types_m, only: obj_pyrolysis, obj_darcy, obj_ablation, obj_mesh_motion
  use FUSS_Lib_BC_Fluxes,  only: Face_Index, Compute_Modfm

  implicit none
  private
  public :: BC_Wall_Ablating

  !> Stefan-Boltzmann constant as PATO evaluates it (OpenFOAM-7,
  !> pi^2/60 k^4/(hbar^3 c^2) with the SI values of etc/controlDict); FUSS's
  !> sigma_SB = 5.67e-8 (BC 304/305) differs by 7.2e-5.
  real(R8), parameter :: SIGMA_PATO = 5.6704085581d-8
  real(R8), parameter :: CHEM_ON = 1.0d0 - 1.0d-12      ! chem >= CHEM_ON is on (PATO tests chemistryOn == 1)
  real(R8), parameter :: DT_FD = 0.5d0                  ! [K] finite-difference half-step of dF/dT_w

contains

  subroutine BC_Wall_Ablating ( Im, Jm, Km, Fm, Blk, bc, Ovar, count_clamps )
    use FUSS_Lib_Solid,      only: co_Kappa
    use FUSS_Lib_Pyrolysis,  only: Surface_Optics, Char_Enthalpy, Gas_State
    use FUSS_Lib_Bprime,     only: Bprime_Lookup
    use FUSS_Lib_BC_Fluxes_Wall_HeatTransfer, only: BC_Wall_HeatTransfer
    use FUSS_Lib_BC_Fluxes_Wall_Temperature,  only: BC_Wall_Temperature
    implicit none
    integer,               intent(in)    :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout) :: Blk
    type(FUSS_bc_type),    intent(inout) :: bc
    real(R8), dimension(2), intent(out)  :: Ovar
    logical,               intent(in)    :: count_clamps
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Dir, Fi, Fj, Fk, mat, it, ib, nmax
    real(R8) :: Normal(3), Mmet(3,3), Tc, kappa, G, rhoe, beta, beta_new, mdot_pg, p_gas
    real(R8) :: Tw, lo, hi, F0, Fp, Fm_, dFdT, Tnew, tol, terms(8), dbeta, Bpg, Bch, hw, mdot_ch, tq(8), Bq(4)
    logical  :: on, frac, pyro, consume, converged, cl(4), clamped(4), clq(4)

    call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, Dir, Im, Jm, Km, Fi, Fj, Fk )
    Normal = Blk%dir(Dir)%f(Fi,Fj,Fk)%n
    Mmet   = Blk%M(Im,Jm,Km)%c

    pyro = obj_pyrolysis%enabled
    Tc   = Blk%T(Im,Jm,Km)
    if ( pyro ) then
      kappa = Blk%k_s(Im,Jm,Km)
      mat   = nint( Blk%matID(Im,Jm,Km) )
    else
      call co_Kappa( Blk%matID(Im,Jm,Km), Tc, kappa )
      mat = 0
    endif
    G = 2.0d0 * kappa * dot_product( Mmet(Dir,:), Normal )

    on   = ( bc%chem >= CHEM_ON )
    frac = ( .not. on ) .and. ( bc%chem > 0.0d0 )
    consume = ( trim(obj_mesh_motion%law) == 'nodeshift' ) .and. obj_mesh_motion%enabled .and. ( .not. obj_ablation%no_char )

    ! Gas mass flux leaving through the face (outward positive): type 2 owns
    ! bc%mdot_pg (Open_Face_Fluxes, once per step); type 1 reads the surface
    ! cell's line integral of the previous stage and keeps a copy for the
    ! wall output; no pyrolysis, no gas.
    mdot_pg = 0.0d0
    if ( pyro ) then
      if ( obj_darcy%enabled ) then
        mdot_pg = bc%mdot_pg
      else
        mdot_pg = Blk%mdot_pg(Im,Jm,Km)
        bc%mdot_pg = mdot_pg
      endif
    endif
    p_gas = bc%pw
    if ( pyro .and. .not. obj_darcy%enabled ) p_gas = obj_pyrolysis%p_gas

    ! ---- chemistry off with a zero emissivity: the BC 303 closed form
    if ( bc%chem == 0.0d0 .and. Optics_Zero( mat, bc%qrad_env ) ) then
      call BC_Wall_HeatTransfer ( Im, Jm, Km, Fm, Blk, bc%hconv_env, bc%Tedge, 0.0d0, Ovar )
      Tw = Ovar(1)
      bc%Tw_prev = Tw ; bc%qw_last = Ovar(2)
      bc%mdot_ch = 0.0d0 ; bc%sdot = 0.0d0 ; bc%Bprime = 0.0d0 ; bc%Bpg = 0.0d0 ; bc%Bch = 0.0d0 ; bc%hw = 0.0d0
      if ( obj_ablation%sdot_override >= 0.0d0 ) bc%sdot = obj_ablation%sdot_override
      bc%seb = 0.0d0
      bc%seb(1) = G * ( Tw - Tc )
      bc%seb(2) = bc%hconv_env * ( bc%Tedge - Tw )
      bc%seb(8) = bc%seb(1) - bc%seb(2)
      bc%seb_iter = 0
      return
    endif

    ! ---- the surface energy balance
    rhoe = bc%rhoeUeCH0
    if ( on ) rhoe = max( rhoe, 1.0d-6 )
    beta = bc%beta
    if ( .not. ( beta > 0.0d0 ) ) beta = 1.0d0
    tol  = obj_ablation%seb_tol
    nmax = obj_ablation%seb_max_iter

    if ( pyro ) then
      lo = max( Tc - 500.0d0, real(obj_pyrolysis%mat(mat)%T_brk_lo, R8) + 1.0d0 )
      hi = real( obj_pyrolysis%mat(mat)%T_brk_hi, R8 ) - 1.0d0
    else
      lo = max( Tc - 500.0d0, 1.0d0 )
      hi = Tc + 1.0d4
    endif
    Tw = bc%Tw_prev
    if ( .not. ( Tw > 0.0d0 ) ) Tw = Tc
    Tw = min( max( Tw, lo ), hi )

    ! A safeguarded Newton on T_w; at every iterate beta is first brought to
    ! its Picard fixed point at that T_w, so F is always the same function of
    ! T_w and the signs of F at the iterates keep [lo, hi] a bracket of its
    ! root (one Picard step per Newton step, the first version, moved the root
    ! between iterates and left the bracket stale: 23 iterations, halving, on
    ! ATC 2.1 at 1.43 s). F is increasing in T_w (conduction, re-radiation and
    ! h_w all grow with it). Converged when, at the current iterate, the step
    ! it would take is <= seb-tol T_w and |F| <= seb-tol max|term| (the Picard
    ! has already met |d beta| <= seb-tol); that iterate is returned, so the
    ! stored F and terms are those of the returned T_w and beta (section 3.2).
    converged = .false.
    dbeta = 0.0d0
    do it = 1, nmax
      if ( on ) then
        do ib = 1, nmax
          call Balance ( Tw, beta, F0, terms, Bpg, Bch, hw, mdot_ch, cl )
          beta_new = Blowing( bc%lambda, mdot_pg, rhoe, beta, Bch )
          dbeta = beta_new - beta
          beta  = beta_new
          if ( abs(dbeta) <= tol ) exit
        enddo
        if ( abs(dbeta) > tol ) exit
      endif
      call Balance ( Tw + DT_FD, beta, Fp,  tq, Bq(1), Bq(2), Bq(3), Bq(4), clq )
      call Balance ( Tw - DT_FD, beta, Fm_, tq, Bq(1), Bq(2), Bq(3), Bq(4), clq )
      call Balance ( Tw, beta, F0, terms, Bpg, Bch, hw, mdot_ch, cl )
      dFdT = ( Fp - Fm_ ) / ( 2.0d0 * DT_FD )
      if ( F0 == 0.0d0 ) then
        converged = .true.
      elseif ( dFdT > 0.0d0 ) then
        converged = abs(F0) <= dFdT * tol * Tw .and. abs(F0) <= tol * maxval( abs(terms(1:7)) )
      endif
      if ( converged ) exit
      if ( F0 > 0.0d0 ) then
        hi = min( hi, Tw )
      else
        lo = max( lo, Tw )
      endif
      Tnew = -1.0d0
      if ( dFdT > 0.0d0 ) Tnew = Tw - F0 / dFdT
      if ( .not. ( Tnew > lo .and. Tnew < hi ) ) Tnew = 0.5d0 * ( lo + hi )
      Tw = Tnew
    enddo
    ! the counters and the abort mode see the queries of the returned state
    ! only, not the bracketing iterates
    clamped = cl

    if ( .not. converged ) then
      write(*,'(A,I0,A)') ' [ERROR] BC ', bc%type, ': the surface energy balance did not converge (SEB did not converge)'
      write(*,'(A,4I6,A,I0)') '         face (i, j, k, f) =', Im, Jm, Km, Fm, '  iterations ', nmax
      write(*,'(A,ES14.6,A,ES14.6,A,ES14.6)') '         T_w = ', Tw, '  T_c = ', Tc, '  F = ', F0
      write(*,'(A,8ES12.4)') '         terms (cond, conv, adv_pg, adv_ch, alpha q, eps T_inf^4, -eps T_w^4, F):', terms
      if ( abs(dbeta) > tol ) write(*,'(A,ES14.6,A,ES12.4)') '         the blowing correction did not settle: beta = ', &
        beta, '  last change ', dbeta
      error stop 3
    endif

    if ( trim(obj_ablation%bprime_out_of_range) == 'abort' .and. ( clamped(1) .or. clamped(2) .or. clamped(4) ) ) then
      write(*,'(A,I0,A)') ' [ERROR] BC ', bc%type, ': the B'' table query left the table (bprime-out-of-range = abort)'
      write(*,'(A,4I6)') '         face (i, j, k, f) =', Im, Jm, Km, Fm
      write(*,'(A,ES14.6,A,ES14.6,A,ES14.6)') '         p_w = ', bc%pw, '  B''g = ', Bpg, '  T_w = ', Tw
      error stop 3
    endif
    if ( count_clamps ) then
      where ( clamped ) bc%nclamp = bc%nclamp + 1
    endif

    ! Residual flux and wall output: BC 302's kernel at the converged T_w
    call BC_Wall_Temperature ( Im, Jm, Km, Fm, Blk, Tw, Ovar )

    ! State of the converged balance
    bc%Tw_prev  = Tw
    bc%qw_last  = Ovar(2)
    bc%seb      = terms
    bc%seb(8)   = F0
    bc%seb_iter = it
    bc%hw       = hw
    bc%Bpg      = Bpg
    bc%Bch      = Bch
    bc%mdot_ch  = mdot_ch
    bc%Bprime   = Bpg + Bch
    if ( on ) then
      bc%beta    = beta
      bc%Bpg_lag = Bpg
    endif
    bc%sdot = 0.0d0
    if ( consume .and. pyro ) bc%sdot = mdot_ch / Blk%rho_s(Im,Jm,Km)
    if ( obj_ablation%sdot_override >= 0.0d0 ) bc%sdot = obj_ablation%sdot_override   ! TEST ONLY

  contains

    !> F(T_w) at a given beta, its eight terms, and the table values used.
    subroutine Balance ( Tw_, beta_, F, t, Bpg_, Bch_, hw_, mdc_, cl_ )
      real(R8), intent(in)  :: Tw_, beta_
      real(R8), intent(out) :: F, t(8), Bpg_, Bch_, hw_, mdc_
      logical,  intent(out) :: cl_(4)
      ! Local
      real(R8) :: eps_s, alpha_s, Wg, hpg, mug, hch, Bcv, Bg_look, rc

      cl_ = .false.
      eps_s = 0.0d0 ; alpha_s = 0.0d0
      if ( pyro ) call Surface_Optics ( mat, Tw_, Blk%chi(:,Im,Jm,Km), eps_s, alpha_s )
      t = 0.0d0
      Bpg_ = 0.0d0 ; Bch_ = 0.0d0 ; hw_ = 0.0d0 ; mdc_ = 0.0d0
      t(1) = G * ( Tw_ - Tc )
      if ( on ) then
        Bpg_ = mdot_pg / ( rhoe * beta_ )
        Bg_look = max( Bpg_, 0.0d0 )
        rc = rhoe * beta_
        call Bprime_Lookup ( bc%ibprime, bc%pw, Bg_look, Tw_, rc, Bcv, hw_, cl_ )
        if ( consume ) then
          Bch_ = Bcv
          mdc_ = Bcv * rhoe * beta_
        endif
        t(2) = rhoe * beta_ * ( bc%hr - hw_ )
        call Gas_State ( mat, p_gas, Tw_, Wg, hpg, mug )
        t(3) = mdot_pg * ( hpg - hw_ )
        if ( consume ) then
          hch = Char_Enthalpy ( mat, Tw_ )
          t(4) = mdc_ * ( hch - hw_ )
        endif
      else
        t(2) = bc%hconv_env * ( bc%Tedge - Tw_ )
        if ( frac .and. consume .and. bc%ibprime > 0 ) then
          ! the 60-60.1 s ramp (P-D2): PATO keeps receding with B'_c at the
          ! lagged B'_pg and the frozen beta; the balance itself is chem-off
          rc = bc%rhoeUeCH0 * bc%beta
          call Bprime_Lookup ( bc%ibprime, bc%pw, max( bc%Bpg_lag, 0.0d0 ), Tw_, rc, Bcv, hw_, cl_ )
          Bpg_ = bc%Bpg_lag
          Bch_ = Bcv
          mdc_ = Bcv * bc%rhoeUeCH0 * bc%beta
        endif
      endif
      t(5) = alpha_s * bc%qrad_env
      t(6) = eps_s * SIGMA_PATO * bc%Tinf**4
      t(7) = - eps_s * SIGMA_PATO * Tw_**4
      F = t(1) - t(2) - t(3) - t(4) - t(5) - t(6) - t(7)
    end subroutine Balance

  end subroutine BC_Wall_Ablating


  !> PATO's constant-lambda blowing correction (constantLambdaBlowingCorrectionModel.C:79-140):
  !> B' = mdot_pg/(rhoeUeCH0 beta) + B'_c, beta = ln(1 + 2 lambda B')/(2 lambda B'),
  !> 1 when lambda = 0, B' = 0, mdot_pg <= 1e-6 or rhoeUeCH0 <= 1e-4.
  pure real(R8) function Blowing ( lambda, mdot_pg, rhoe, beta, Bc ) result ( b )
    implicit none
    real(R8), intent(in) :: lambda, mdot_pg, rhoe, beta, Bc
    real(R8) :: Bp
    Bp = mdot_pg / ( rhoe * beta ) + Bc
    if ( lambda == 0.0d0 .or. Bp == 0.0d0 .or. mdot_pg <= 1.0d-6 .or. rhoe <= 1.0d-4 ) then
      b = 1.0d0
    else
      b = log( 1.0d0 + 2.0d0 * lambda * Bp ) / ( 2.0d0 * lambda * Bp )
    endif
  end function Blowing


  !> .true. when the face radiates nothing and absorbs nothing: no pyrolysis
  !> material (kinematic mode), or a material whose emissivity columns are
  !> zero everywhere and whose absorptivity is zero or meets no q_rad.
  pure logical function Optics_Zero ( m, qrad ) result ( z )
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: qrad
    z = .true.
    if ( m == 0 ) return
    associate ( mt => obj_pyrolysis%mat(m) )
      z = all( mt%em_v == 0.0d0 ) .and. all( mt%em_c == 0.0d0 )
      if ( z .and. qrad /= 0.0d0 ) z = all( mt%ab_v == 0.0d0 ) .and. all( mt%ab_c == 0.0d0 )
    end associate
  end function Optics_Zero

end module FUSS_Lib_BC_Fluxes_Wall_Ablating
