!> Pyrolysis property model and kinetics closures (plan 06 §2, §3.2, §3.3), and
!> the type-2 (Darcy) closures of plan 07 in their own section at the end.
!>
!> BIT-IDENTITY NOTE (plan 07): every routine above the type-2 section is plan
!> 06's text verbatim. The build contracts floating-point expressions into
!> fused multiply-adds across inlined calls, so restructuring a shared routine
!> (an added branch between a product and the subtraction that followed it)
!> moved type-1 results by 1 ulp (seen 2026-10-08); the type-2 closures are
!> separate routines (Energy_Of2, Temperature_Of2, ...) and the call sites
!> select them with obj_darcy%enabled.
!>
!> Everything here is cell-local and thread-safe: the only module state read
!> is obj_pyrolysis (materials, tables, the flattened reaction list), which is
!> fixed after setup. Table lookups are linear between integer rows exactly as
!> Lib_Solid's co_Cp; the per-material Newton bracket [T_brk_lo, T_brk_hi] is
!> the intersection of every table range a lookup touches (computed once by
!> the loader), so no lookup can read outside its table.
module FUSS_Lib_Pyrolysis
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Config_Types_m, only: obj_pyrolysis, obj_darcy
  use FUSS_Parameters_m,   only: R_univ

  implicit none
  private
  public :: Tab_Lookup, Tau_Of, Solid_State, Gas_State, Energy_Of, Energy_Tau, Temperature_Of
  public :: Reaction_Rate, Advance_Reaction
  public :: Gas_Lookup, Energy_Of2, Temperature_Of2, Perm_Tensor, Darcy_Cell_Coeffs
  public :: Surface_Optics, Char_Enthalpy

  real(R8), parameter :: NEWTON_DT_TOL = 1.0d-10   ! [K] a Newton step at most this long terminates
  real(R8), parameter :: NEWTON_REL_F  = 1.0d-13   ! |f(T)| <= NEWTON_REL_F*|E| terminates too
  integer,  parameter :: NEWTON_CAP    = 50

contains

  !> Linear interpolation between integer rows (the co_Cp pattern): v is a
  !> table indexed by integer temperature from `lo`; T must satisfy
  !> lo <= T and int(T)+1 <= ubound(v).
  pure function Tab_Lookup ( v, lo, T ) result ( val )
    implicit none
    integer,  intent(in) :: lo
    real(R8), intent(in) :: v(lo:)
    real(R8), intent(in) :: T
    real(R8) :: val
    ! Local
    integer :: Tm

    Tm  = int(T)
    val = v(Tm) + (T - Tm) * (v(Tm+1) - v(Tm))

  end function Tab_Lookup


  !> Virgin fraction tau = sum_r w_r (1 - chi_r) over the material's reactions
  !> (1 virgin, 0 charred); identically 1 for a material with rho_v = rho_c
  !> (no reaction, or every F = 0), where the weights are undefined.
  pure function Tau_Of ( m, chi ) result ( tau )
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: chi(:)
    real(R8) :: tau
    ! Local
    integer :: r

    tau = 1.0d0
    if ( obj_pyrolysis%mat(m)%hbar_guard ) return
    tau = 0.0d0
    do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
      tau = tau + obj_pyrolysis%reac(r)%w * ( 1.0d0 - chi(r) )
    enddo
    tau = min( max( tau, 0.0d0 ), 1.0d0 )

  end function Tau_Of


  !> Derived solid state of a cell: rho_s = sum_i eps_i0 rho_i(chi) and PATO's
  !> blends of the virgin and char tables with the weight tau rho_v/max(rho_s, rho_c)
  !> (PorousMaterialPropertiesModel.C:188-191; plan 03 §12.4), plus the gas
  !> porosity eps_pg(tau).
  subroutine Solid_State ( m, T, chi, rho_s, cp_s, k_s, tau, eps_pg )
    implicit none
    integer,  intent(in)  :: m
    real(R8), intent(in)  :: T, chi(:)
    real(R8), intent(out) :: rho_s, cp_s, k_s, tau, eps_pg
    ! Local
    integer  :: r
    real(R8) :: wv

    associate ( mt => obj_pyrolysis%mat(m) )
      tau = Tau_Of( m, chi )
      rho_s = mt%rho_v
      do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
        rho_s = rho_s - obj_pyrolysis%reac(r)%eps_rho_F * chi(r)
      enddo
      wv = tau * mt%rho_v / max( rho_s, mt%rho_c )
      cp_s = Tab_Lookup( mt%cp_v, mt%Tlo_s, T ) * wv + Tab_Lookup( mt%cp_c, mt%Tlo_s, T ) * ( 1.0d0 - wv )
      k_s  = Tab_Lookup( mt%k_v,  mt%Tlo_s, T ) * wv + Tab_Lookup( mt%k_c,  mt%Tlo_s, T ) * ( 1.0d0 - wv )
      eps_pg = mt%eps_g_c + ( mt%eps_g_v - mt%eps_g_c ) * tau
    end associate

  end subroutine Solid_State


  !> Gas properties at (p, T): linear in T between integer rows, linear in p
  !> between the two bracketing p-levels. Outside the table: error stop 3.
  subroutine Gas_State ( m, p, T, W, h_pg, mu, cell )
    implicit none
    integer,  intent(in)  :: m
    real(R8), intent(in)  :: p, T
    real(R8), intent(out) :: W, h_pg, mu
    integer,  intent(in), optional :: cell(4)
    ! Local
    integer  :: ip, ipp
    real(R8) :: wp

    associate ( mt => obj_pyrolysis%mat(m) )
      if ( T < mt%Tlo_g .or. T >= mt%Thi_g .or. p < mt%p_lev(1) .or. p > mt%p_lev(mt%Np) ) then
        write(*,'(A)') ' [ERROR] pyrolysis: gas state outside the gas table'
        write(*,'(A,A,A,ES14.6,A,ES14.6)') '         material ', trim(mt%name), ': p = ', p, '  T = ', T
        write(*,'(A,I0,A,I0,A,ES14.6,A,ES14.6)') '         table: T in [', mt%Tlo_g, ', ', mt%Thi_g, ') K, p in [', &
          mt%p_lev(1), ', ', mt%p_lev(mt%Np)
        if ( present(cell) ) write(*,'(A,4I6)') '         at b, i, j, k =', cell
        error stop 3
      endif
      if ( mt%Np == 1 ) then
        ip = 1 ; ipp = 1 ; wp = 0.0d0
      else
        ip = mt%Np - 1
        do while ( ip > 1 .and. p < mt%p_lev(ip) )
          ip = ip - 1
        enddo
        ipp = ip + 1
        wp  = ( p - mt%p_lev(ip) ) / ( mt%p_lev(ipp) - mt%p_lev(ip) )
      endif
      W    = ( 1.0d0 - wp ) * Tab_Lookup( mt%W_g(:,ip),  mt%Tlo_g, T ) + wp * Tab_Lookup( mt%W_g(:,ipp),  mt%Tlo_g, T )
      h_pg = ( 1.0d0 - wp ) * Tab_Lookup( mt%h_g(:,ip),  mt%Tlo_g, T ) + wp * Tab_Lookup( mt%h_g(:,ipp),  mt%Tlo_g, T )
      mu   = ( 1.0d0 - wp ) * Tab_Lookup( mt%mu_g(:,ip), mt%Tlo_g, T ) + wp * Tab_Lookup( mt%mu_g(:,ipp), mt%Tlo_g, T )
    end associate

  end subroutine Gas_State


  !> Stored energy per unit volume at a given virgin fraction (type 1, plan 06 §2
  !> with decision D1): E = H_s = tau rho_v h_v(T) + (1 - tau) rho_c h_c(T).
  !> No gas storage; p is accepted for interface stability (plan 07 adds the
  !> type-2 term eps_pg (rho_pg h_pg - p)).
  pure function Energy_Tau ( m, T, tau ) result ( E )
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: T, tau
    real(R8) :: E

    associate ( mt => obj_pyrolysis%mat(m) )
      E = tau * mt%rho_v * Tab_Lookup( mt%h_v, mt%Tlo_s, T ) + ( 1.0d0 - tau ) * mt%rho_c * Tab_Lookup( mt%h_c, mt%Tlo_s, T )
    end associate

  end function Energy_Tau


  function Energy_Of ( m, T, chi, p ) result ( E )
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: T, chi(:), p
    real(R8) :: E

    E = Energy_Tau( m, T, Tau_Of( m, chi ) )
    if ( .false. ) E = E + 0.0d0 * p

  end function Energy_Of


  !> Slope of E on the 1 K segment holding T, per kelvin: under type 1 E is
  !> piecewise linear in T, so this derivative is exact there and one Newton
  !> step from inside the right segment lands on the root to round-off.
  pure function Energy_Slope ( m, T, tau ) result ( fp )
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: T, tau
    real(R8) :: fp
    ! Local
    integer :: Tm

    Tm = int(T)
    associate ( mt => obj_pyrolysis%mat(m) )
      fp = tau * mt%rho_v * ( mt%h_v(Tm+1) - mt%h_v(Tm) ) + ( 1.0d0 - tau ) * mt%rho_c * ( mt%h_c(Tm+1) - mt%h_c(Tm) )
    end associate

  end function Energy_Slope


  !> Invert E -> T by a bracketed Newton iteration on f(T) = E(T; chi) - E
  !> (plan 06 §3.2). Bracket = the material's [T_brk_lo, T_brk_hi]; outside it
  !> the run stops (error stop 3) with the cell, its energy and its chi. A
  !> Newton candidate outside the current bracket is replaced by a bisection
  !> step, which does not count as converged; convergence is a Newton step of
  !> at most NEWTON_DT_TOL or |f| <= NEWTON_REL_F |E|.
  subroutine Temperature_Of ( m, E, chi, p, T_guess, T, b, i, j, k )
    implicit none
    integer,  intent(in)  :: m, b, i, j, k
    real(R8), intent(in)  :: E, chi(:), p, T_guess
    real(R8), intent(out) :: T
    ! Local
    real(R8) :: Ta, Tb, Ea, Eb, fT, fp, Tnew, tau, dT
    integer  :: it
    logical  :: converged

    tau = Tau_Of( m, chi )
    Ta  = real( obj_pyrolysis%mat(m)%T_brk_lo, R8 )
    Tb  = real( obj_pyrolysis%mat(m)%T_brk_hi, R8 )
    Ea  = Energy_Tau( m, Ta, tau )
    Eb  = Energy_Tau( m, Tb, tau )
    if ( E < Ea .or. E > Eb ) then
      call Report( 'energy outside table range' )
      error stop 3
    endif

    T = min( max( T_guess, Ta ), Tb )
    converged = .false.
    do it = 1, NEWTON_CAP
      fT = Energy_Tau( m, T, tau ) - E
      if ( abs(fT) <= NEWTON_REL_F * abs(E) ) then
        converged = .true.
        exit
      endif
      if ( fT > 0.0d0 ) then
        Tb = T
      else
        Ta = T
      endif
      fp = Energy_Slope( m, T, tau )
      if ( fp > 0.0d0 ) then
        Tnew = T - fT / fp
        dT   = abs( Tnew - T )
        if ( dT <= NEWTON_DT_TOL ) then
          ! A Newton step at most NEWTON_DT_TOL long: converged on the Newton
          ! iterate. Tested BEFORE the bracket test on purpose: near the zero
          ! of E (TACOT: ~890 K) the correction can fall below the ulp of T,
          ! so Tnew == T, which equals the bracket end just set from the sign
          ! of fT; the strict bracket test would then bisect to the middle of
          ! the bracket and the iteration would cycle until the cap (seen at
          ! ATC1 step 37, cell 98).
          if ( Tnew >= Ta .and. Tnew <= Tb ) T = Tnew
          converged = .true.
          exit
        endif
      else
        Tnew = Ta - 1.0d0      ! a non-increasing table: force the bisection branch
      endif
      if ( Tnew > Ta .and. Tnew < Tb ) then
        T = Tnew
      else
        T = 0.5d0 * ( Ta + Tb )
      endif
    enddo

    if ( .not. converged ) then
      call Report( 'Newton did not converge' )
      error stop 3
    endif
    if ( .false. ) T = T + 0.0d0 * p

  contains

    subroutine Report ( what )
      implicit none
      character(len=*), intent(in) :: what
      write(*,'(A)') ' [ERROR] pyrolysis: '//what
      write(*,'(A,4I6)') '         b, i, j, k =', b, i, j, k
      write(*,'(A,ES24.16,A,F10.3)') '         E = ', E, '  T_guess = ', T_guess
      write(*,'(A,*(ES14.6))') '         chi = ', chi( obj_pyrolysis%mat_sub_offset(m-1)+1 : obj_pyrolysis%mat_sub_offset(m) )
      write(*,'(A,I0,A,I0,A)') '         bracket [', obj_pyrolysis%mat(m)%T_brk_lo, ', ', obj_pyrolysis%mat(m)%T_brk_hi, '] K'
    end subroutine Report

  end subroutine Temperature_Of


  !> Arrhenius rate of reaction r at T: kappa = T^b A exp(-E/(R T)), zero below
  !> the reaction's threshold temperature (PATO's per-reaction switch).
  pure function Reaction_Rate ( m, r, T ) result ( kappa )
    implicit none
    integer,  intent(in) :: m, r
    real(R8), intent(in) :: T
    real(R8) :: kappa

    kappa = 0.0d0
    associate ( rc => obj_pyrolysis%reac(r) )
      if ( T < rc%T_thr .or. rc%A <= 0.0d0 ) return
      kappa = rc%A * exp( - rc%E / ( R_univ * T ) )
      if ( rc%b /= 0.0d0 ) kappa = kappa * T**rc%b
    end associate
    if ( .false. ) kappa = kappa + 0.0d0 * m

  end function Reaction_Rate


  !> Exact propagator of d chi/dt = kappa (1 - chi)^m over dt at frozen kappa
  !> (plan 06 §3.3): (1 - chi)^(1-m) = (1 - chi0)^(1-m) - (1 - m) kappa dt for
  !> m /= 1, (1 - chi) = (1 - chi0) exp(-kappa dt) for m = 1; chi <= 1 and
  !> chi >= chi0 always. Returns chi0 unchanged for kappa = 0 or chi0 >= 1
  !> (for m > 1 the bracket (1 - chi0)^(1-m) is 0^negative at chi0 = 1).
  pure function Advance_Reaction ( chi0, kappa, m_exp, dt ) result ( chi )
    implicit none
    real(R8), intent(in) :: chi0, kappa, m_exp, dt
    real(R8) :: chi
    ! Local
    real(R8) :: y0, y, br

    chi = chi0
    if ( kappa <= 0.0d0 .or. chi0 >= 1.0d0 ) return
    y0 = 1.0d0 - chi0
    if ( m_exp == 1.0d0 ) then
      y = y0 * exp( - kappa * dt )
    else
      br = y0**( 1.0d0 - m_exp ) - ( 1.0d0 - m_exp ) * kappa * dt
      if ( br <= 0.0d0 ) then
        y = 0.0d0
      else
        y = br**( 1.0d0 / ( 1.0d0 - m_exp ) )
      endif
    endif
    chi = min( 1.0d0 - y, 1.0d0 )
    chi = max( chi, chi0 )

  end function Advance_Reaction


  ! ===========================================================================
  !  TYPE-2 (DARCY) CLOSURES, plan 07
  !
  !  The stored energy under gas-closure = type2 (plan 07 section 3.5; PATO
  !  PyrolysisEnergyModel.C:56) is
  !      E = H_s(T, tau) + G h_g(p, T) - eps_pg(tau) p,
  !  with G [kg/m3] the stored gas mass per unit volume. G is a STATE of the
  !  scheme, not eps_pg rho_g(p, T) re-evaluated from the equation of state:
  !  the pressure solve fixes it (G = eta* p^(n+1), FUSS_Mod_Darcy), so the gas
  !  mass in E is exactly the mass the mass equation accounted for. With G from
  !  the equation of state, any mismatch between the stored mass and the outflow
  !  of the solve is multiplied by the ABSOLUTE gas enthalpy (TACOT: -7e6 J/kg)
  !  and the scheme is datum-dependent at first order in dt (measured on the
  !  exact 0-D gates, 2026-10-08). Under [FUSS-Darcy] quasi-steady = true the
  !  gas storage is dropped from the mass equation and, for the same
  !  consistency, from E: E = H_s.
  ! ===========================================================================

  !> The interpolation of Gas_State without its range check, plus the 1 K
  !> segment slopes dW/dT, dh/dT at (p, T). Callers stay inside the table: the
  !> Newton bracket does so in T, and every p reaches here after a checked
  !> Gas_State call on the same cell in the same step.
  pure subroutine Gas_Lookup ( m, p, T, W, h_pg, mu, dW, dh )
    implicit none
    integer,  intent(in)  :: m
    real(R8), intent(in)  :: p, T
    real(R8), intent(out) :: W, h_pg, mu, dW, dh
    ! Local
    integer  :: ip, ipp, Tm
    real(R8) :: wp

    associate ( mt => obj_pyrolysis%mat(m) )
      if ( mt%Np == 1 ) then
        ip = 1 ; ipp = 1 ; wp = 0.0d0
      else
        ip = mt%Np - 1
        do while ( ip > 1 .and. p < mt%p_lev(ip) )
          ip = ip - 1
        enddo
        ipp = ip + 1
        wp  = ( p - mt%p_lev(ip) ) / ( mt%p_lev(ipp) - mt%p_lev(ip) )
      endif
      W    = ( 1.0d0 - wp ) * Tab_Lookup( mt%W_g(:,ip),  mt%Tlo_g, T ) + wp * Tab_Lookup( mt%W_g(:,ipp),  mt%Tlo_g, T )
      h_pg = ( 1.0d0 - wp ) * Tab_Lookup( mt%h_g(:,ip),  mt%Tlo_g, T ) + wp * Tab_Lookup( mt%h_g(:,ipp),  mt%Tlo_g, T )
      mu   = ( 1.0d0 - wp ) * Tab_Lookup( mt%mu_g(:,ip), mt%Tlo_g, T ) + wp * Tab_Lookup( mt%mu_g(:,ipp), mt%Tlo_g, T )
      Tm = int(T)
      dW = ( 1.0d0 - wp ) * ( mt%W_g(Tm+1,ip) - mt%W_g(Tm,ip) ) + wp * ( mt%W_g(Tm+1,ipp) - mt%W_g(Tm,ipp) )
      dh = ( 1.0d0 - wp ) * ( mt%h_g(Tm+1,ip) - mt%h_g(Tm,ip) ) + wp * ( mt%h_g(Tm+1,ipp) - mt%h_g(Tm,ipp) )
    end associate

  end subroutine Gas_Lookup


  !> E(T; tau, p, G) = H_s + G h_g(p, T) - eps_pg(tau) p (H_s alone under quasi-steady).
  pure function Energy_Full2 ( m, T, tau, p, G ) result ( E )
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: T, tau, p, G
    real(R8) :: E
    ! Local
    real(R8) :: W, h, mu, dW, dh, eps

    E = Energy_Tau( m, T, tau )
    if ( obj_darcy%quasi_steady ) return
    call Gas_Lookup( m, p, T, W, h, mu, dW, dh )
    associate ( mt => obj_pyrolysis%mat(m) )
      eps = mt%eps_g_c + ( mt%eps_g_v - mt%eps_g_c ) * tau
    end associate
    E = E + ( G * h - eps * p )

  end function Energy_Full2


  !> dE/dT at fixed (tau, p, G): the solid segment slope plus G dh_g/dT.
  pure function Energy_Slope2 ( m, T, tau, p, G ) result ( fp )
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: T, tau, p, G
    real(R8) :: fp
    ! Local
    integer  :: Tm
    real(R8) :: W, h, mu, dW, dh

    Tm = int(T)
    associate ( mt => obj_pyrolysis%mat(m) )
      fp = tau * mt%rho_v * ( mt%h_v(Tm+1) - mt%h_v(Tm) ) + ( 1.0d0 - tau ) * mt%rho_c * ( mt%h_c(Tm+1) - mt%h_c(Tm) )
    end associate
    if ( obj_darcy%quasi_steady ) return
    call Gas_Lookup( m, p, T, W, h, mu, dW, dh )
    fp = fp + G * dh

  end function Energy_Slope2


  !> Type-2 stored energy of a cell state (T, chi, p, G).
  function Energy_Of2 ( m, T, chi, p, G ) result ( E )
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: T, chi(:), p, G
    real(R8) :: E

    E = Energy_Full2( m, T, Tau_Of( m, chi ), p, G )

  end function Energy_Of2


  !> Type-2 inversion E -> T at fixed (chi, p, G): the bracketed Newton of
  !> Temperature_Of on Energy_Full2 (E is no longer piecewise linear in T, a
  !> few more Newton steps).
  subroutine Temperature_Of2 ( m, E, chi, p, G, T_guess, T, b, i, j, k )
    implicit none
    integer,  intent(in)  :: m, b, i, j, k
    real(R8), intent(in)  :: E, chi(:), p, G, T_guess
    real(R8), intent(out) :: T
    ! Local
    real(R8) :: Ta, Tb, Ea, Eb, fT, fp, Tnew, tau, dT
    integer  :: it
    logical  :: converged

    tau = Tau_Of( m, chi )
    Ta  = real( obj_pyrolysis%mat(m)%T_brk_lo, R8 )
    Tb  = real( obj_pyrolysis%mat(m)%T_brk_hi, R8 )
    Ea  = Energy_Full2( m, Ta, tau, p, G )
    Eb  = Energy_Full2( m, Tb, tau, p, G )
    if ( E < Ea .or. E > Eb ) then
      call Report( 'energy outside table range' )
      error stop 3
    endif

    T = min( max( T_guess, Ta ), Tb )
    converged = .false.
    do it = 1, NEWTON_CAP
      fT = Energy_Full2( m, T, tau, p, G ) - E
      if ( abs(fT) <= NEWTON_REL_F * abs(E) ) then
        converged = .true.
        exit
      endif
      if ( fT > 0.0d0 ) then
        Tb = T
      else
        Ta = T
      endif
      fp = Energy_Slope2( m, T, tau, p, G )
      if ( fp > 0.0d0 ) then
        Tnew = T - fT / fp
        dT   = abs( Tnew - T )
        if ( dT <= NEWTON_DT_TOL ) then
          if ( Tnew >= Ta .and. Tnew <= Tb ) T = Tnew
          converged = .true.
          exit
        endif
      else
        Tnew = Ta - 1.0d0
      endif
      if ( Tnew > Ta .and. Tnew < Tb ) then
        T = Tnew
      else
        T = 0.5d0 * ( Ta + Tb )
      endif
    enddo

    if ( .not. converged ) then
      call Report( 'Newton did not converge' )
      error stop 3
    endif

  contains

    subroutine Report ( what )
      implicit none
      character(len=*), intent(in) :: what
      write(*,'(A)') ' [ERROR] pyrolysis (type 2): '//what
      write(*,'(A,4I6)') '         b, i, j, k =', b, i, j, k
      write(*,'(A,ES24.16,A,F10.3,A,ES14.6,A,ES14.6)') '         E = ', E, '  T_guess = ', T_guess, '  p = ', p, '  G = ', G
      write(*,'(A,*(ES14.6))') '         chi = ', chi( obj_pyrolysis%mat_sub_offset(m-1)+1 : obj_pyrolysis%mat_sub_offset(m) )
      write(*,'(A,I0,A,I0,A)') '         bracket [', obj_pyrolysis%mat(m)%T_brk_lo, ', ', obj_pyrolysis%mat(m)%T_brk_hi, '] K'
    end subroutine Report

  end subroutine Temperature_Of2


  !> Permeability tensor of material m at virgin fraction tau (thesis eqs
  !> 3.33/3.34; plan 07 section 3.1): K = K_c + (K_v - K_c) tau, plus the
  !> Klinkenberg term (mu/p) B when [FUSS-Darcy] klinkenberg = true.
  pure function Perm_Tensor ( m, tau, p, mu ) result ( K )
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: tau, p, mu
    real(R8) :: K(3,3)

    associate ( mt => obj_pyrolysis%mat(m) )
      K = mt%Kp_c + ( mt%Kp_v - mt%Kp_c ) * tau
      if ( obj_darcy%klinkenberg ) K = K + ( mu / p ) * ( mt%Bp_c + ( mt%Bp_v - mt%Bp_c ) * tau )
    end associate

  end function Perm_Tensor


  !> Lagged coefficients of the pressure equation for one cell (plan 07
  !> section 3.3; PATO DarcyLawMassModel.C:113-119): eta = eps_pg W/(R T)
  !> (eta p = eps_pg rho_g, the stored gas mass of the equation of state) and
  !> Gamma = rho_g K/mu, with the gas state range-checked (error stop 3 names
  !> the cell). Also returns eps_pg, rho_g and the gas enthalpy.
  subroutine Darcy_Cell_Coeffs ( m, p, T, chi, cell, eta, Gam, eps, rho_g, h_pg )
    implicit none
    integer,  intent(in)  :: m, cell(4)
    real(R8), intent(in)  :: p, T, chi(:)
    real(R8), intent(out) :: eta, Gam(3,3), eps, rho_g, h_pg
    ! Local
    real(R8) :: rho_s, cp_s, k_s, tau, W, mu

    call Solid_State( m, T, chi, rho_s, cp_s, k_s, tau, eps )
    call Gas_State( m, p, T, W, h_pg, mu, cell )
    rho_g = p * W / ( R_univ * T )
    eta   = eps * W / ( R_univ * T )
    Gam   = ( rho_g / mu ) * Perm_Tensor( m, tau, p, mu )

  end subroutine Darcy_Cell_Coeffs


  ! ---------------------------------------------------------------------------
  ! Ablating wall (plan 08 section 2.1), new routines only; nothing above is
  ! touched (see the bit-identity note in the header).
  ! ---------------------------------------------------------------------------

  !> Surface emissivity and absorptivity of a cell's material at the WALL
  !> temperature: the virgin and char table values blended with the weight
  !> tau rho_v/max(rho_s, rho_c) of the cell's state, exactly as PATO blends
  !> them on its boundary faces (PorousMaterialPropertiesModel.C:471-478) and
  !> as Solid_State blends cp and k.
  pure subroutine Surface_Optics ( m, Tw, chi, eps_s, alpha_s )
    implicit none
    integer,  intent(in)  :: m
    real(R8), intent(in)  :: Tw, chi(:)
    real(R8), intent(out) :: eps_s, alpha_s
    ! Local
    integer  :: r
    real(R8) :: tau, rho_s, wv

    associate ( mt => obj_pyrolysis%mat(m) )
      tau = Tau_Of( m, chi )
      rho_s = mt%rho_v
      do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
        rho_s = rho_s - obj_pyrolysis%reac(r)%eps_rho_F * chi(r)
      enddo
      wv = tau * mt%rho_v / max( rho_s, mt%rho_c )
      eps_s   = Tab_Lookup( mt%em_v, mt%Tlo_s, Tw ) * wv + Tab_Lookup( mt%em_c, mt%Tlo_s, Tw ) * ( 1.0d0 - wv )
      alpha_s = Tab_Lookup( mt%ab_v, mt%Tlo_s, Tw ) * wv + Tab_Lookup( mt%ab_c, mt%Tlo_s, Tw ) * ( 1.0d0 - wv )
    end associate

  end subroutine Surface_Optics


  !> Char enthalpy h_c(T) [J/kg] of a material (the char table), the h_ch,w of
  !> the surface energy balance's char-advection term.
  pure function Char_Enthalpy ( m, T ) result ( h )
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: T
    real(R8) :: h
    h = Tab_Lookup( obj_pyrolysis%mat(m)%h_c, obj_pyrolysis%mat(m)%Tlo_s, T )
  end function Char_Enthalpy

end module FUSS_Lib_Pyrolysis
