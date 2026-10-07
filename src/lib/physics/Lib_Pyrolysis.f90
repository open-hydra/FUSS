!> Pyrolysis property model and kinetics closures (plan 06 §2, §3.2, §3.3).
!>
!> Everything here is cell-local and thread-safe: the only module state read
!> is obj_pyrolysis (materials, tables, the flattened reaction list), which is
!> fixed after setup. Table lookups are linear between integer rows exactly as
!> Lib_Solid's co_Cp; the per-material Newton bracket [T_brk_lo, T_brk_hi] is
!> the intersection of every table range a lookup touches (computed once by
!> the loader), so no lookup can read outside its table.
module FUSS_Lib_Pyrolysis
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Config_Types_m, only: obj_pyrolysis
  use FUSS_Parameters_m,   only: R_univ

  implicit none
  private
  public :: Tab_Lookup, Tau_Of, Solid_State, Gas_State, Energy_Of, Energy_Tau, Temperature_Of
  public :: Reaction_Rate, Advance_Reaction

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

end module FUSS_Lib_Pyrolysis
