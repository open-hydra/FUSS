!> Pyrolysis stage work, energy budget and field-list glue (plan 06 §3.1,
!> §3.3, §3.4, §3.6, §4.3).
!>
!> Pyrolysis_Stage runs once per Runge-Kutta stage between BC_Fluxes and the
!> user hook, inside the parallel region Explicit_Step opens (so every loop is
!> `!$omp do`, never `!$omp parallel do`), and does four things in order:
!>   (a) the exact closed-form kinetics increment of every reaction from the
!>       start-of-step chi at the stage temperature (blk%dchi),
!>   (b) the type-1 line integral of the gas production along the outflow
!>       direction (blk%mdot_pg), serial along each line and parallel across
!>       lines, so the result does not depend on the thread count,
!>   (c) the donor-cell gas-advection flux added to blk%R, which telescopes
!>       over the interior faces of a line,
!>   (d) the energy-budget sums, single-threaded in a fixed order, folded into
!>       the (a, a_n) accumulators with the scheme's own RK_stage weights, so
!>       the discrete balances of gate 4.3 are identities of the update.
!> (d) must read blk%R AFTER (c) and BEFORE RK_Newstate (Compute_Cons_NoIRS
!> overwrites R in place with the scaled residual).
module FUSS_Mod_Pyrolysis
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Pyrolysis_Assign_Setup, Pyrolysis_Init_Cells, Pyrolysis_Stage, Write_Energy_Budget, Pyrolysis_Residual
  public :: Pyrolysis_Budget_Rewind
  public :: Convert_T_To, Convert_T_From, Convert_Chi_To, Convert_Chi_From   ! plan 08: the pyrolysis remap

  ! Energy-budget accumulators (plan 06 §4.3): Q_out = int sum_cells R dt (net
  ! outgoing flux, conduction + BC + gas advection), Q_gas_out = int sum_lines
  ! Mdot_n h_pg(T_n) dt, M_gas_out = int sum_lines Mdot_n dt, Q_src = int sum
  ! qvol V dt. Each is a pair (a, a_n): a_n = a saved at stage 1, then
  ! a = RK_stage(irk, n_RK, a, a_n, stage increment). They start at 0 on a
  ! restart as well (the budget file is appended; its Q columns count from
  ! the restart time).
  real(R8), save :: Q_out = 0.0d0,     Q_out_n = 0.0d0
  real(R8), save :: Q_gas_out = 0.0d0, Q_gas_out_n = 0.0d0
  real(R8), save :: M_gas_out = 0.0d0, M_gas_out_n = 0.0d0
  real(R8), save :: Q_src = 0.0d0,     Q_src_n = 0.0d0
  logical,  save :: budget_file_started = .false.

contains

  !> Derived settings and the field-list converters. Called from Assign_Setup
  !> after Field_List_Init, when the feature is on.
  subroutine Pyrolysis_Assign_Setup ()
    use FUSS_Config_Types_m, only: obj_pyrolysis, obj_darcy
    use FUSS_Field_List_m,   only: Field_List_Set_Converters
    implicit none

    select case ( trim(obj_pyrolysis%outflow_direction) )
      case ('+i') ; obj_pyrolysis%out_dir = 1 ; obj_pyrolysis%out_sign =  1
      case ('-i') ; obj_pyrolysis%out_dir = 1 ; obj_pyrolysis%out_sign = -1
      case ('+j') ; obj_pyrolysis%out_dir = 2 ; obj_pyrolysis%out_sign =  1
      case ('-j') ; obj_pyrolysis%out_dir = 2 ; obj_pyrolysis%out_sign = -1
      case ('+k') ; obj_pyrolysis%out_dir = 3 ; obj_pyrolysis%out_sign =  1
      case ('-k') ; obj_pyrolysis%out_dir = 3 ; obj_pyrolysis%out_sign = -1
      case default; obj_pyrolysis%out_dir = 0 ; obj_pyrolysis%out_sign =  0   ! refused by Check_Pyrolysis_Compatibility
    end select

    obj_pyrolysis%description = 'type-1 gas closure, outflow '//trim(obj_pyrolysis%outflow_direction)
    if ( obj_darcy%enabled ) obj_pyrolysis%description = 'type-2 gas closure (Darcy pressure equation, [FUSS-Darcy])'

    call Field_List_Set_Converters( 'T',   Convert_T_To,   Convert_T_From )
    call Field_List_Set_Converters( 'chi', Convert_Chi_To, Convert_Chi_From )

  end subroutine Pyrolysis_Assign_Setup


  !> Caches of every cell from (T, chi), the initial-condition bracket check,
  !> and the field-list round-trip self-check (plan 06 §3.1, O7). Called from
  !> Setup_Data_Structure once T, matID and chi are in place.
  subroutine Pyrolysis_Init_Cells ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_pyrolysis
    use FUSS_Lib_Pyrolysis,  only: Solid_State
    use FUSS_Mod_MPI,        only: is_local_block, mpi_is_root
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b, i, j, k, m
    real(R8) :: tau, eps_pg, T0, dT, dchi, max_dT, max_dchi
    real(R8), allocatable :: chi0(:), q(:)

    max_dT = 0.0d0 ; max_dchi = 0.0d0
    allocate( chi0( obj_pyrolysis%nsub ), q( max(obj_pyrolysis%nsub, 1) ) )

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        m = nint( domain%blk(b)%matID(i,j,k) )
        if ( m < 1 .or. m > obj_pyrolysis%nmat ) then
          write(*,'(A,4I6)') ' [ERROR] pyrolysis: matID outside the material list at b, i, j, k =', b, i, j, k
          error stop 1
        endif
        if ( domain%blk(b)%T(i,j,k) < obj_pyrolysis%mat(m)%T_brk_lo .or. &
             domain%blk(b)%T(i,j,k) > obj_pyrolysis%mat(m)%T_brk_hi ) then
          write(*,'(A)') ' [ERROR] pyrolysis: initial temperature outside the Newton bracket of its material'
          write(*,'(A,4I6,A,F10.3,A,I0,A,I0,A)') '         b, i, j, k =', b, i, j, k, '  T = ', domain%blk(b)%T(i,j,k), &
            '  bracket [', obj_pyrolysis%mat(m)%T_brk_lo, ', ', obj_pyrolysis%mat(m)%T_brk_hi, '] K'
          error stop 1
        endif
        call Solid_State( m, domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), &
                          domain%blk(b)%rho_s(i,j,k), domain%blk(b)%cp_s(i,j,k), domain%blk(b)%k_s(i,j,k), tau, eps_pg )

        ! Round trip through the conserved quantities: chi first (T's
        ! converter reads the cell's chi), then T; both must come back.
        T0   = domain%blk(b)%T(i,j,k)
        chi0 = domain%blk(b)%chi(:,i,j,k)
        call Convert_Chi_To  ( domain%blk(b), i, j, k, q )
        call Convert_Chi_From( domain%blk(b), i, j, k, q )
        call Convert_T_To    ( domain%blk(b), i, j, k, q(1:1) )
        call Convert_T_From  ( domain%blk(b), i, j, k, q(1:1) )
        dT   = abs( domain%blk(b)%T(i,j,k) - T0 )
        dchi = maxval( abs( domain%blk(b)%chi(:,i,j,k) - chi0 ) )
        max_dT = max( max_dT, dT ) ; max_dchi = max( max_dchi, dchi )
        if ( dT > 1.0d-9 .or. dchi > 1.0d-14 ) then
          write(*,'(A)') ' [ERROR] pyrolysis: field-list round trip does not reproduce the state'
          write(*,'(A,4I6,A,ES12.4,A,ES12.4)') '         b, i, j, k =', b, i, j, k, '  |dT| = ', dT, '  |dchi| = ', dchi
          error stop 3
        endif
        domain%blk(b)%T(i,j,k)     = T0
        domain%blk(b)%chi(:,i,j,k) = chi0
      enddo ; enddo ; enddo
    enddo

    if ( mpi_is_root ) write(*,'(A,ES10.2,A,ES10.2)') '   Field list round-trip OK, max dT = ', max_dT, ' K, max dchi = ', max_dchi

  end subroutine Pyrolysis_Init_Cells


  ! --- field-list converters (plan 06 §3.1): conserved quantity per unit volume --

  !> T -> E(T, chi) (reads the cell's chi: remap chi before T). Type 2 (plan 07):
  !> the cell's pressure and stored gas mass G = eta_st p_pg enter E.
  subroutine Convert_T_To ( blk, i, j, k, q )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    use FUSS_Config_Types_m,   only: obj_pyrolysis, obj_darcy
    use FUSS_Lib_Pyrolysis,    only: Energy_Of, Energy_Of2
    implicit none
    type(FUSS_block_type), intent(inout) :: blk
    integer,  intent(in)    :: i, j, k
    real(R8), intent(inout) :: q(:)
    if ( obj_darcy%enabled ) then
      q(1) = Energy_Of2( nint(blk%matID(i,j,k)), blk%T(i,j,k), blk%chi(:,i,j,k), blk%p_pg(i,j,k), &
                         blk%eta_st(i,j,k) * blk%p_pg(i,j,k) )
      return
    endif
    q(1) = Energy_Of( nint(blk%matID(i,j,k)), blk%T(i,j,k), blk%chi(:,i,j,k), obj_pyrolysis%p_gas )
  end subroutine Convert_T_To

  subroutine Convert_T_From ( blk, i, j, k, q )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    use FUSS_Config_Types_m,   only: obj_pyrolysis, obj_darcy
    use FUSS_Lib_Pyrolysis,    only: Temperature_Of, Temperature_Of2
    implicit none
    type(FUSS_block_type), intent(inout) :: blk
    integer,  intent(in)    :: i, j, k
    real(R8), intent(inout) :: q(:)
    real(R8) :: T
    if ( obj_darcy%enabled ) then
      call Temperature_Of2( nint(blk%matID(i,j,k)), q(1), blk%chi(:,i,j,k), blk%p_pg(i,j,k), &
                            blk%eta_st(i,j,k) * blk%p_pg(i,j,k), blk%T(i,j,k), T, 0, i, j, k )
      blk%T(i,j,k) = T
      return
    endif
    call Temperature_Of( nint(blk%matID(i,j,k)), q(1), blk%chi(:,i,j,k), obj_pyrolysis%p_gas, blk%T(i,j,k), T, 0, i, j, k )
    blk%T(i,j,k) = T
  end subroutine Convert_T_From

  !> chi -> remaining sub-phase mass eps_i0 rho_i0 F_r (1 - chi_r); entries of
  !> other materials are 0.
  subroutine Convert_Chi_To ( blk, i, j, k, q )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    use FUSS_Config_Types_m,   only: obj_pyrolysis
    implicit none
    type(FUSS_block_type), intent(inout) :: blk
    integer,  intent(in)    :: i, j, k
    real(R8), intent(inout) :: q(:)
    integer :: m, r
    m = nint( blk%matID(i,j,k) )
    q = 0.0d0
    do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
      q(r) = obj_pyrolysis%reac(r)%eps_rho_F * ( 1.0d0 - blk%chi(r,i,j,k) )
    enddo
  end subroutine Convert_Chi_To

  subroutine Convert_Chi_From ( blk, i, j, k, q )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    use FUSS_Config_Types_m,   only: obj_pyrolysis
    implicit none
    type(FUSS_block_type), intent(inout) :: blk
    integer,  intent(in)    :: i, j, k
    real(R8), intent(inout) :: q(:)
    integer :: m, r
    m = nint( blk%matID(i,j,k) )
    do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
      if ( obj_pyrolysis%reac(r)%eps_rho_F > 0.0d0 ) &
        blk%chi(r,i,j,k) = min( max( 1.0d0 - q(r) / obj_pyrolysis%reac(r)%eps_rho_F, 0.0d0 ), 1.0d0 )
    enddo
  end subroutine Convert_Chi_From


  !> Area of the face of direction d at face index (i,j,k): face l of
  !> direction d sits between cells l and l+1 along d (Mod_Allocate_Data).
  pure function Face_Area ( blk, d, i, j, k ) result ( A )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    implicit none
    type(FUSS_block_type), intent(in) :: blk
    integer, intent(in) :: d, i, j, k
    real(R8) :: A
    A = blk%dir(d)%f(i,j,k)%A
  end function Face_Area


  !> The per-stage work (a)-(d) described in the module header.
  subroutine Pyrolysis_Stage ( domain, i_rk )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_pyrolysis, obj_time_scheme, obj_darcy
    use FUSS_Lib_Pyrolysis,  only: Reaction_Rate, Advance_Reaction, Gas_State
    use FUSS_Lib_RK,         only: RK_stage
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in) :: i_rk
    ! Local
    integer  :: b, i, j, k, r, m, d, sg, t1, t2, nd, n1, n2, s, c, b1, b2, ic(3), id(3), n_rk
    real(R8) :: dt, kap, T, pi_tot, Mdot, Mdot_new, h_prev, W, h_pg, mu, A, p_gas
    real(R8) :: sumR, sumS, sumQg, sumMg

    ! Type 2 (plan 07): its own routine, so that this one stays plan 06's text
    if ( obj_darcy%enabled ) then
      call Pyrolysis_Stage2( domain, i_rk )
      return
    endif

    n_rk  = obj_time_scheme%n_RK
    p_gas = obj_pyrolysis%p_gas
    d  = obj_pyrolysis%out_dir
    sg = obj_pyrolysis%out_sign
    t1 = mod( d, 3 ) + 1
    t2 = mod( d + 1, 3 ) + 1

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle

      ! (a) kinetics: full-step increment from the start-of-step chi at the
      !     stage temperature; reactions of other materials keep dchi = 0.
      !$omp do collapse(3) schedule(static) private(i,j,k,m,r,dt,T,kap)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        m  = nint( domain%blk(b)%matID(i,j,k) )
        dt = domain%blk(b)%dtlocal(i,j,k)
        T  = domain%blk(b)%T(i,j,k)
        do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
          kap = Reaction_Rate( m, r, T )
          domain%blk(b)%dchi(r,i,j,k) = Advance_Reaction( domain%blk(b)%chiO(r,i,j,k), kap, obj_pyrolysis%reac(r)%m, dt ) &
                                      - domain%blk(b)%chiO(r,i,j,k)
        enddo
      enddo ; enddo ; enddo
      !$omp end do

      ! (b)+(c) running integral along every line of the outflow direction,
      !     cells numbered from the sealed end; mdot_pg = Mdot / A(downstream
      !     face); donor-cell advection flux into R. pi_tot is the realised
      !     stage production (increment over dt), so gate 4.3 closes to round-off.
      nd = domain%blk(b)%dim(d)
      n1 = domain%blk(b)%dim(t1)
      n2 = domain%blk(b)%dim(t2)
      !$omp do collapse(2) schedule(static) private(b1,b2,s,c,ic,id,i,j,k,m,dt,pi_tot,r,Mdot,Mdot_new,h_prev,W,h_pg,mu,A)
      do b2 = 1, n2
      do b1 = 1, n1
        Mdot   = 0.0d0
        h_prev = 0.0d0
        do s = 1, nd
          if ( sg > 0 ) then
            c = s
          else
            c = nd + 1 - s
          endif
          ic(d) = c ; ic(t1) = b1 ; ic(t2) = b2
          i = ic(1) ; j = ic(2) ; k = ic(3)
          m  = nint( domain%blk(b)%matID(i,j,k) )
          dt = domain%blk(b)%dtlocal(i,j,k)
          pi_tot = 0.0d0
          do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
            pi_tot = pi_tot + obj_pyrolysis%reac(r)%eps_rho_F * domain%blk(b)%dchi(r,i,j,k)
          enddo
          pi_tot   = pi_tot / dt
          Mdot_new = Mdot + pi_tot * domain%blk(b)%vol(i,j,k)
          call Gas_State( m, p_gas, domain%blk(b)%T(i,j,k), W, h_pg, mu, [b, i, j, k] )
          id = ic
          if ( sg < 0 ) id(d) = c - 1
          A = Face_Area( domain%blk(b), d, id(1), id(2), id(3) )
          domain%blk(b)%mdot_pg(i,j,k) = Mdot_new / A
          domain%blk(b)%R(i,j,k) = domain%blk(b)%R(i,j,k) + Mdot_new * h_pg - Mdot * h_prev
          Mdot   = Mdot_new
          h_prev = h_pg
        enddo
      enddo ; enddo
      !$omp end do
    enddo

    ! (d) budget sums: one thread, fixed order (blocks, then k, j, i; lines in
    !     the same order), read from R after the advection flux was added and
    !     before RK_Newstate rescales R in place.
    !$omp single
    sumR = 0.0d0 ; sumS = 0.0d0 ; sumQg = 0.0d0 ; sumMg = 0.0d0
    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        dt   = domain%blk(b)%dtlocal(i,j,k)
        sumR = sumR + domain%blk(b)%R(i,j,k) * dt
        sumS = sumS + domain%blk(b)%qvol(i,j,k) * domain%blk(b)%vol(i,j,k) * dt
      enddo ; enddo ; enddo
      nd = domain%blk(b)%dim(d)
      n1 = domain%blk(b)%dim(t1)
      n2 = domain%blk(b)%dim(t2)
      do b2 = 1, n2
      do b1 = 1, n1
        if ( sg > 0 ) then
          c = nd
        else
          c = 1
        endif
        ic(d) = c ; ic(t1) = b1 ; ic(t2) = b2
        i = ic(1) ; j = ic(2) ; k = ic(3)
        id = ic
        if ( sg < 0 ) id(d) = c - 1
        A  = Face_Area( domain%blk(b), d, id(1), id(2), id(3) )
        m  = nint( domain%blk(b)%matID(i,j,k) )
        dt = domain%blk(b)%dtlocal(i,j,k)
        call Gas_State( m, p_gas, domain%blk(b)%T(i,j,k), W, h_pg, mu, [b, i, j, k] )
        Mdot  = domain%blk(b)%mdot_pg(i,j,k) * A
        sumQg = sumQg + Mdot * h_pg * dt
        sumMg = sumMg + Mdot * dt
      enddo ; enddo
    enddo
    if ( i_rk == 1 ) then
      Q_out_n = Q_out ; Q_gas_out_n = Q_gas_out ; M_gas_out_n = M_gas_out ; Q_src_n = Q_src
    endif
    Q_out     = RK_stage( i_rk, n_rk, Q_out,     Q_out_n,     sumR  )
    Q_gas_out = RK_stage( i_rk, n_rk, Q_gas_out, Q_gas_out_n, sumQg )
    M_gas_out = RK_stage( i_rk, n_rk, M_gas_out, M_gas_out_n, sumMg )
    Q_src     = RK_stage( i_rk, n_rk, Q_src,     Q_src_n,     sumS  )
    !$omp end single

  end subroutine Pyrolysis_Stage


  !> Type-2 stage work (plan 07 section 3.5): (a) the kinetics increments as
  !> Pyrolysis_Stage, then the gas-enthalpy flux div(h_g Gamma grad p) of the
  !> solved pressure into R (FUSS_Mod_Darcy) in place of the type-1 line
  !> integral, then (d) the budget sums with the gas leaving through the open
  !> faces (bc%mdot_pg, bc%edot_pg; bc order is the fixed order of bc.txt).
  subroutine Pyrolysis_Stage2 ( domain, i_rk )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_pyrolysis, obj_time_scheme
    use FUSS_Lib_Pyrolysis,  only: Reaction_Rate, Advance_Reaction
    use FUSS_Lib_RK,         only: RK_stage
    use FUSS_Lib_BC_Fluxes,  only: Face_Index
    use FUSS_Mod_Darcy,      only: Darcy_Energy_Flux, Darcy_Is_Open
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer, intent(in) :: i_rk
    ! Local
    integer  :: b, i, j, k, r, m, n_rk, ii, Dir, Fi, Fj, Fk
    real(R8) :: dt, kap, T, A
    real(R8) :: sumR, sumS, sumQg, sumMg

    n_rk = obj_time_scheme%n_RK

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      !$omp do collapse(3) schedule(static) private(i,j,k,m,r,dt,T,kap)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        m  = nint( domain%blk(b)%matID(i,j,k) )
        dt = domain%blk(b)%dtlocal(i,j,k)
        T  = domain%blk(b)%T(i,j,k)
        do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
          kap = Reaction_Rate( m, r, T )
          domain%blk(b)%dchi(r,i,j,k) = Advance_Reaction( domain%blk(b)%chiO(r,i,j,k), kap, obj_pyrolysis%reac(r)%m, dt ) &
                                      - domain%blk(b)%chiO(r,i,j,k)
        enddo
      enddo ; enddo ; enddo
      !$omp end do
    enddo

    call Darcy_Energy_Flux( domain )

    !$omp single
    sumR = 0.0d0 ; sumS = 0.0d0 ; sumQg = 0.0d0 ; sumMg = 0.0d0
    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        dt   = domain%blk(b)%dtlocal(i,j,k)
        sumR = sumR + domain%blk(b)%R(i,j,k) * dt
        sumS = sumS + domain%blk(b)%qvol(i,j,k) * domain%blk(b)%vol(i,j,k) * dt
      enddo ; enddo ; enddo
    enddo
    do ii = 1, domain%nbound
      if ( .not. Darcy_Is_Open( domain%bc(ii)%type ) ) cycle
      b = domain%bc(ii)%b
      if ( .not. is_local_block(b) ) cycle
      i = domain%bc(ii)%i ; j = domain%bc(ii)%j ; k = domain%bc(ii)%k
      call Face_Index( domain%bc(ii)%f, Dir, i, j, k, Fi, Fj, Fk )
      A  = domain%blk(b)%dir(Dir)%f(Fi,Fj,Fk)%A
      dt = domain%blk(b)%dtlocal(i,j,k)
      sumQg = sumQg + domain%bc(ii)%edot_pg * A * dt
      sumMg = sumMg + domain%bc(ii)%mdot_pg * A * dt
    enddo
    if ( i_rk == 1 ) then
      Q_out_n = Q_out ; Q_gas_out_n = Q_gas_out ; M_gas_out_n = M_gas_out ; Q_src_n = Q_src
    endif
    Q_out     = RK_stage( i_rk, n_rk, Q_out,     Q_out_n,     sumR  )
    Q_gas_out = RK_stage( i_rk, n_rk, Q_gas_out, Q_gas_out_n, sumQg )
    M_gas_out = RK_stage( i_rk, n_rk, M_gas_out, M_gas_out_n, sumMg )
    Q_src     = RK_stage( i_rk, n_rk, Q_src,     Q_src_n,     sumS  )
    !$omp end single

  end subroutine Pyrolysis_Stage2


  !> Budget accumulators back to their start-of-step values (saved at stage
  !> 1), for the second pass of a redone step (plan 07: a new run's first
  !> step). Inside the parallel region.
  subroutine Pyrolysis_Budget_Rewind ()
    implicit none
    !$omp single
    Q_out = Q_out_n ; Q_gas_out = Q_gas_out_n ; M_gas_out = M_gas_out_n ; Q_src = Q_src_n
    !$omp end single
  end subroutine Pyrolysis_Budget_Rewind


  !> Convergence-residual contribution of the reaction state: the sum over
  !> cells and reactions of (chi - chiO)^2, added to the solver's total
  !> (Compute_Residual's sum of (T - TO)^2). Without it a run whose
  !> temperature is exactly constant -- the isothermal gate -- reports a zero
  !> residual at the first step and stops as converged while the reactions
  !> still advance. Called outside the parallel region, like Compute_Residual.
  subroutine Pyrolysis_Residual ( domain, total )
    use FUSS_Advanced_Types_m
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    real(R8), intent(inout) :: total
    ! Local
    integer  :: b, i, j, k, r
    real(R8) :: res

    res = 0.0d0
    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        do r = 1, size( domain%blk(b)%chi, 1 )
          res = res + ( domain%blk(b)%chi(r,i,j,k) - domain%blk(b)%chiO(r,i,j,k) )**2
        enddo
      enddo ; enddo ; enddo
    enddo
    total = total + res

  end subroutine Pyrolysis_Residual


  !> One line of OUTPUT/<prefix>energy-budget.txt: t E_tot M_s Q_out Q_gas_out
  !> M_gas_out Q_src (ES24.16). State sums in a fixed order on one thread, so
  !> the file is reproducible across thread counts like the solution. Called
  !> on root, outside any parallel region: once at the end of setup and at
  !> every solution write.
  subroutine Write_Energy_Budget ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_pyrolysis, obj_sim_param, obj_darcy, obj_ablation, obj_mesh_motion
    use FUSS_Global_m,       only: FUSS_phase_prefix
    use FUSS_Lib_Pyrolysis,  only: Energy_Of, Energy_Of2
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    integer  :: b, i, j, k, m, u, ios
    real(R8) :: E_tot, M_s, V, M_g, G
    character(len=256) :: file
    logical  :: abl

    ! A receding ablating wall (plan 08 section 2.9) adds what its faces swept
    ! out, so the balances still close on a moving surface; every other run
    ! keeps its columns byte for byte.
    abl = obj_ablation%enabled .and. trim(obj_mesh_motion%law) == 'nodeshift'

    ! Type 2 (plan 07): E with the stored gas, and an eighth column M_g = sum
    ! V G (the stored gas mass), so that the type-2 mass balance
    ! M_s + M_g + M_gas_out = const can be checked.
    if ( obj_darcy%enabled ) then
      E_tot = 0.0d0 ; M_s = 0.0d0 ; M_g = 0.0d0
      do b = 1, domain%nb
        if ( .not. is_local_block(b) ) cycle
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          m = nint( domain%blk(b)%matID(i,j,k) )
          V = domain%blk(b)%vol(i,j,k)
          G = domain%blk(b)%eta_st(i,j,k) * domain%blk(b)%p_pg(i,j,k)
          E_tot = E_tot + V * Energy_Of2( m, domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), domain%blk(b)%p_pg(i,j,k), G )
          M_s   = M_s   + V * domain%blk(b)%rho_s(i,j,k)
          M_g   = M_g   + V * G
        enddo ; enddo ; enddo
      enddo
      file = 'OUTPUT/'//trim(FUSS_phase_prefix)//'energy-budget.txt'
      if ( .not. budget_file_started ) then
        if ( obj_sim_param%newrun ) then
          open( newunit=u, file=trim(file), status='replace', action='write', iostat=ios )
          if ( ios == 0 .and. abl ) write(u,'(A)') '# t  E_tot  M_s  Q_out  Q_gas_out  M_gas_out  Q_src  M_g  E_abl  M_abl  G_abl   '// &
            '(SI; Q_* and *_abl count from the start of this run)'
          if ( ios == 0 .and. .not. abl ) write(u,'(A)') &
            '# t  E_tot  M_s  Q_out  Q_gas_out  M_gas_out  Q_src  M_g   (SI; Q_* count from the start of this run)'
        else
          open( newunit=u, file=trim(file), status='unknown', position='append', action='write', iostat=ios )
        endif
        budget_file_started = .true.
      else
        open( newunit=u, file=trim(file), status='old', position='append', action='write', iostat=ios )
      endif
      if ( ios /= 0 ) then
        write(*,'(A)') ' [ERROR] pyrolysis: cannot open '//trim(file)
        error stop 1
      endif
      if ( abl ) then
        write(u,'(11ES24.16)') domain%time, E_tot, M_s, Q_out, Q_gas_out, M_gas_out, Q_src, M_g, &
          obj_ablation%E_abl, obj_ablation%M_abl, obj_ablation%G_abl
      else
        write(u,'(8ES24.16)') domain%time, E_tot, M_s, Q_out, Q_gas_out, M_gas_out, Q_src, M_g
      endif
      close(u)
      return
    endif

    E_tot = 0.0d0 ; M_s = 0.0d0
    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        m = nint( domain%blk(b)%matID(i,j,k) )
        V = domain%blk(b)%vol(i,j,k)
        E_tot = E_tot + V * Energy_Of( m, domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), obj_pyrolysis%p_gas )
        M_s   = M_s   + V * domain%blk(b)%rho_s(i,j,k)
      enddo ; enddo ; enddo
    enddo

    file = 'OUTPUT/'//trim(FUSS_phase_prefix)//'energy-budget.txt'
    if ( .not. budget_file_started ) then
      if ( obj_sim_param%newrun ) then
        open( newunit=u, file=trim(file), status='replace', action='write', iostat=ios )
        if ( ios == 0 .and. abl ) write(u,'(A)') '# t  E_tot  M_s  Q_out  Q_gas_out  M_gas_out  Q_src  E_abl  M_abl   '// &
          '(SI; Q_* and *_abl count from the start of this run)'
        if ( ios == 0 .and. .not. abl ) write(u,'(A)') &
          '# t  E_tot  M_s  Q_out  Q_gas_out  M_gas_out  Q_src   (SI; Q_* count from the start of this run)'
      else
        open( newunit=u, file=trim(file), status='unknown', position='append', action='write', iostat=ios )
      endif
      budget_file_started = .true.
    else
      open( newunit=u, file=trim(file), status='old', position='append', action='write', iostat=ios )
    endif
    if ( ios /= 0 ) then
      write(*,'(A)') ' [ERROR] pyrolysis: cannot open '//trim(file)
      error stop 1
    endif
    if ( abl ) then
      write(u,'(9ES24.16)') domain%time, E_tot, M_s, Q_out, Q_gas_out, M_gas_out, Q_src, obj_ablation%E_abl, obj_ablation%M_abl
    else
      write(u,'(7ES24.16)') domain%time, E_tot, M_s, Q_out, Q_gas_out, M_gas_out, Q_src
    endif
    close(u)

  end subroutine Write_Energy_Budget

end module FUSS_Mod_Pyrolysis
