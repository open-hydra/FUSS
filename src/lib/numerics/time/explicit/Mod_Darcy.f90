!> The Darcy pressure equation of the type-2 gas closure (plan 07 §3.3-3.7):
!> the first equation in FUSS that is not advanced by the explicit RK loop.
!>
!> Once per material step, after Copy_State and before the RK stages,
!> Solve_Pressure solves the backward-Euler, coefficient-lagged mass balance
!>     (eta^n p - eta^(n-1) p^n) V/dt - sum_f A_f (Gamma_f grad p)_f . n_f = pi V
!> (PATO DarcyLawMassModel.C:84-90 with fvm::ddt(Eta, p) and the lagged
!> Gamma of :113-119) for p^(n+1); every RK stage then reads p^(n+1) as frozen
!> data (the gas storage in E, the gas enthalpy flux grad.(h_g Gamma grad p)).
!> Two sweep kinds behind [FUSS-Darcy] solver: line-implicit Thomas sweeps
!> along the open-face direction of a block (`tridiag`, exact in one sweep on
!> a 1-D line) and damped Jacobi (`pseudo`, the fallback of a block whose open
!> faces lie in several directions or none). Both are block-Jacobi over lines
!> or cells of the SAME iterate -- every sweep reads p and writes p_new, then
!> copies -- so the sweep count and the result do not depend on the thread
!> count; the stop test is on the cell mass residual relative to the largest
!> term of the balance (max reductions, order-independent).
!>
!> Everything here runs inside the `!$omp parallel` region of Explicit_Step
!> (orphaned `!$omp do` / `single`, never `parallel do`), except the setup
!> entries, which run serially (an orphaned `!$omp do` outside a parallel
!> region executes on the one thread). No MPI halo exchange exists for p_pg:
!> plan 11 §2.2; pyrolysis is refused with mpi_size_ > 1 (guard 8).
module FUSS_Mod_Darcy
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Darcy_Assign_Setup, Darcy_Setup_Blocks, Darcy_Init_Cells, Solve_Pressure, Fill_Ghost_p
  public :: Darcy_Energy_Flux, Darcy_Residual, Darcy_Update_Velocity, Darcy_Print_Setup, Darcy_Is_Open
  public :: Darcy_First_Step_Pending, Darcy_First_Step_History, Darcy_Restore_State
  public :: Darcy_Mass_Correction

  !> A new run's first step has no history for the mid-step extrapolations
  !> of Solve_Pressure (T_prev = T, eta^(n-1) = eta^n): Explicit_Step runs
  !> it twice, the second time with the history taken from the first pass
  !> (Darcy_First_Step_History). Cleared after that; never set on a restart,
  !> whose history comes from the solution file.
  logical, save :: first_step_pending = .false.

  !> Per-block work arrays of the sweeps: module state, hence SHARED under
  !> OpenMP (the loops that fill them are worksharing loops of one team).
  type :: darcy_work_t
    real(R8), allocatable :: res(:,:,:)    ! cell mass residual r_c [kg/s]: inflows + pi V - storage
    real(R8), allocatable :: dgl(:,:,:)    ! line-implicit diagonal: storage + conductances of the line direction + C_w
    real(R8), allocatable :: dga(:,:,:)    ! Jacobi diagonal: storage + every conductance + C_w
    real(R8), allocatable :: rhs(:,:,:)    ! right-hand side of the line rows
    real(R8), allocatable :: lo(:,:,:)     ! coefficient of p(c-1) along the line (= -C of the lower face)
    real(R8), allocatable :: hi(:,:,:)     ! coefficient of p(c+1) along the line (= -C of the upper face)
    real(R8), allocatable :: g(:,:,:,:)    ! energy pass: cell tensor h_g Gamma (9, row-major)
    real(R8), allocatable :: cq(:,:,:)     ! quasi-steady: gas released but not yet expelled at the step start [kg/m3]
    real(R8), allocatable :: dgx(:,:,:)    ! line-implicit: the part of dgl from the faces across the lines
    real(R8), allocatable :: cx(:,:,:,:)   ! line-implicit: C of the cross faces, (- first, + first, - second, + second)
  end type darcy_work_t
  type(darcy_work_t), allocatable, save :: wk(:)

  ! Sweep control, shared by the team (reductions land here).
  real(R8), save :: resid_max = 0.0d0, resid_blk = 0.0d0, rref_max = 0.0d0, fmax = 0.0d0, smax = 0.0d0
  real(R8), save :: dp_max = 0.0d0, floor_max = 0.0d0
  real(R8), allocatable, save :: resb(:)
  integer,  save :: nsweep = 0
  logical,  save :: done = .false., diag_started = .false.

contains

  ! ---------------------------------------------------------------------------
  !  setup
  ! ---------------------------------------------------------------------------

  !> Ini-derived settings (Assign_Setup, after Field_List_Init): the open
  !> code list and the description. The guard on root has already refused a
  !> bad list; this is the belt on every rank.
  subroutine Darcy_Assign_Setup ()
    use FUSS_Config_Types_m, only: obj_darcy
    use FUSS_Register_Darcy, only: Parse_Open_Codes
    implicit none
    logical :: ok
    character(len=64) :: bad

    call Parse_Open_Codes( obj_darcy%open_bc_types, obj_darcy%open_codes, obj_darcy%n_open, ok, bad )
    if ( .not. ok ) then
      write(*,'(A)') ' [ERROR] [FUSS-Darcy] open-bc-types: bad token "'//trim(bad)//'"'
      error stop 1
    endif
    obj_darcy%description = 'type-2 gas closure: Darcy pressure equation, solver '//trim(obj_darcy%solver)

  end subroutine Darcy_Assign_Setup


  !> .true. when a face of this BC code is open to the gas (p = p_w).
  pure logical function Darcy_Is_Open ( bctype ) result ( open )
    use FUSS_Config_Types_m, only: obj_darcy
    implicit none
    integer, intent(in) :: bctype
    integer :: k
    ! An ablating wall (306, 307) is open to its own pyrolysis gas whatever
    ! open-bc-types says (plan 08, decision P-D6): a face that blows the gas
    ! in its surface balance cannot be sealed in the pressure equation.
    open = ( bctype == 306 .or. bctype == 307 )
    if ( open ) return
    do k = 1, obj_darcy%n_open
      if ( obj_darcy%open_codes(k) == bctype ) open = .true.
    enddo
  end function Darcy_Is_Open


  !> End-of-step coefficient of the stored gas, eta* = 2 eta^n - eta^(n-1)
  !> (plan 07, the scheme of 2026-10-09), limited: eta ~ 1/T, so a cell whose
  !> temperature more than doubles in one step (a step-change wall on a coarse
  !> cell: 298 -> 663 K) makes the linear extrapolation negative, the storage
  !> diagonal negative and the line solve return a negative pressure (Baum &
  !> Atreya case, 2026-10-13). Below eta^n/2 the extrapolation falls back to
  !> eta^n (first order for that step). Assemble and the post-solve update use
  !> this one function, so G^(n+1) = eta* p^(n+1) is the same expression.
  pure real(R8) function Eta_Star ( eta, eta_prev ) result ( es )
    implicit none
    real(R8), intent(in) :: eta, eta_prev
    es = 2.0d0 * eta - eta_prev
    if ( es < 0.5d0 * eta ) es = eta
  end function Eta_Star


  !> Per-block solver choice (plan 07 §3.3 direction rule) and the work
  !> arrays. Called after Setup_BC on every rank. A block whose open faces lie
  !> in exactly one index direction is swept line-implicitly along it when
  !> solver = tridiag; two or three directions, or none, use the Jacobi sweep.
  subroutine Darcy_Setup_Blocks ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_darcy
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    integer :: b, i, d, nd, ni, nj, nk
    logical :: hasdir(3)

    if ( allocated(obj_darcy%blk_dir) ) deallocate( obj_darcy%blk_dir )
    allocate( obj_darcy%blk_dir( domain%nb ) )
    obj_darcy%blk_dir = 0
    if ( allocated(wk) )   deallocate( wk )
    if ( allocated(resb) ) deallocate( resb )
    allocate( wk( domain%nb ), resb( domain%nb ) )
    resb = 0.0d0

    do b = 1, domain%nb
      hasdir = .false.
      do i = 1, domain%nbound
        if ( domain%bc(i)%b /= b ) cycle
        if ( Darcy_Is_Open( domain%bc(i)%type ) ) hasdir( ( domain%bc(i)%f + 1 ) / 2 ) = .true.
      enddo
      nd = count( hasdir )
      if ( trim(obj_darcy%solver) == 'tridiag' .and. nd == 1 ) then
        do d = 1, 3
          if ( hasdir(d) ) obj_darcy%blk_dir(b) = d
        enddo
      endif
      if ( .not. is_local_block(b) ) cycle
      ni = domain%blk(b)%dim(1) ; nj = domain%blk(b)%dim(2) ; nk = domain%blk(b)%dim(3)
      allocate( wk(b)%res(ni,nj,nk), wk(b)%dgl(ni,nj,nk), wk(b)%dga(ni,nj,nk), wk(b)%rhs(ni,nj,nk), &
                wk(b)%lo(ni,nj,nk), wk(b)%hi(ni,nj,nk), wk(b)%g(9,ni,nj,nk), wk(b)%cq(ni,nj,nk), &
                wk(b)%dgx(ni,nj,nk), wk(b)%cx(4,ni,nj,nk) )
      wk(b)%res = 0.0d0 ; wk(b)%dgl = 0.0d0 ; wk(b)%dga = 0.0d0 ; wk(b)%rhs = 0.0d0
      wk(b)%lo  = 0.0d0 ; wk(b)%hi  = 0.0d0 ; wk(b)%g   = 0.0d0 ; wk(b)%cq  = 0.0d0
      wk(b)%dgx = 0.0d0 ; wk(b)%cx  = 0.0d0
    enddo

  end subroutine Darcy_Setup_Blocks


  !> Caches of the initial state (Setup_Data_Structure, after the pyrolysis
  !> caches). The history of the scheme -- eta_prev (eta^(n-1), the
  !> equation-of-state coefficient one step back), eta_st (the coefficient of
  !> the stored gas, G = eta_st p) and T_prev -- comes from the solution file
  !> on a restart (Read_vtk_tec), so that a restarted run continues the
  !> continuous one; a new run (or a file without them: negative sentinels)
  !> starts from the equation of state, eta_prev = eta_st = eta(p0, T0, chi0),
  !> PATO's Eta.oldTime() at start-up. The velocity needs ghosts and metrics
  !> and is filled later by Darcy_Update_Velocity.
  subroutine Darcy_Init_Cells ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_sim_param
    use FUSS_Lib_Pyrolysis,  only: Darcy_Cell_Coeffs
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b, i, j, k, m
    real(R8) :: eta, Gam(3,3), eps, rho_g, h

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        m = nint( domain%blk(b)%matID(i,j,k) )
        call Darcy_Cell_Coeffs( m, domain%blk(b)%p_pg(i,j,k), domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), &
                                [b, i, j, k], eta, Gam, eps, rho_g, h )
        domain%blk(b)%eta_pg(i,j,k)   = eta
        if ( domain%blk(b)%eta_prev(i,j,k) < 0.0d0 ) domain%blk(b)%eta_prev(i,j,k) = eta
        if ( domain%blk(b)%eta_st(i,j,k)   < 0.0d0 ) domain%blk(b)%eta_st(i,j,k)   = eta
        domain%blk(b)%eta_stO(i,j,k)  = domain%blk(b)%eta_st(i,j,k)
        domain%blk(b)%gam_pg(:,i,j,k) = Pack9( Gam )
        domain%blk(b)%pi_pg(i,j,k)    = 0.0d0
        domain%blk(b)%p_pgO(i,j,k)    = domain%blk(b)%p_pg(i,j,k)
      enddo ; enddo ; enddo
    enddo
    first_step_pending = obj_sim_param%newrun

  end subroutine Darcy_Init_Cells


  logical function Darcy_First_Step_Pending () result ( pending )
    implicit none
    pending = first_step_pending
  end function Darcy_First_Step_Pending


  !> History for the second pass of a new run's first step, from the first
  !> pass's end state (T, chi, p_pg) and start state (TO, chiO, p_pgO):
  !> T_prev = 2 TO - T, so that Solve_Pressure's T* = TO + (TO - T_prev)/2
  !> is the predicted mid-step temperature, and eta^(n-1) = 2 eta^n - eta(end),
  !> so that its eta* is the predicted end-of-step coefficient. Inside the
  !> parallel region.
  subroutine Darcy_First_Step_History ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Lib_Pyrolysis, only: Darcy_Cell_Coeffs
    use FUSS_Mod_MPI,       only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b, i, j, k, m
    real(R8) :: eta_end, Gam(3,3), eps, rho_g, h

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      !$omp do collapse(3) schedule(static) private(i,j,k,m,eta_end,Gam,eps,rho_g,h)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        m = nint( domain%blk(b)%matID(i,j,k) )
        call Darcy_Cell_Coeffs( m, domain%blk(b)%p_pg(i,j,k), domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), &
                                [b, i, j, k], eta_end, Gam, eps, rho_g, h )
        domain%blk(b)%T_prev(i,j,k)  = 2.0d0 * domain%blk(b)%TO(i,j,k) - domain%blk(b)%T(i,j,k)
        domain%blk(b)%eta_prev(i,j,k) = 2.0d0 * domain%blk(b)%eta_pg(i,j,k) - eta_end
      enddo ; enddo ; enddo
      !$omp end do
    enddo
    !$omp single
    first_step_pending = .false.
    !$omp end single

  end subroutine Darcy_First_Step_History


  !> Back to the start of the step (T, chi, p_pg, eta_st and the property caches)
  !> for the second pass. Inside the parallel region; ghosts are refilled by
  !> the stage loop.
  subroutine Darcy_Restore_State ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Lib_Pyrolysis, only: Solid_State
    use FUSS_Mod_MPI,       only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b, i, j, k, m
    real(R8) :: tau, eps

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      !$omp do collapse(3) schedule(static) private(i,j,k,m,tau,eps)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        m = nint( domain%blk(b)%matID(i,j,k) )
        domain%blk(b)%T(i,j,k)      = domain%blk(b)%TO(i,j,k)
        domain%blk(b)%chi(:,i,j,k)  = domain%blk(b)%chiO(:,i,j,k)
        domain%blk(b)%p_pg(i,j,k)   = domain%blk(b)%p_pgO(i,j,k)
        domain%blk(b)%eta_st(i,j,k) = domain%blk(b)%eta_stO(i,j,k)
        call Solid_State( m, domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), domain%blk(b)%rho_s(i,j,k), &
                          domain%blk(b)%cp_s(i,j,k), domain%blk(b)%k_s(i,j,k), tau, eps )
      enddo ; enddo ; enddo
      !$omp end do
    enddo

  end subroutine Darcy_Restore_State


  !> End-of-step closure of the gas mass (plan 07, the scheme of 2026-10-09).
  !> The solve produced pi* dt of gas per unit volume, from the kinetics at the
  !> extrapolated mid-step temperature; the RK stages decomposed the solid by
  !> the realised dm_s = sum_r eps_rho_F_r (chi_r^(n+1) - chi_r^n). The
  !> difference (O(dt^3) per step, large when kappa dt ~ 1) goes into the
  !> stored gas, G^(n+1) += dm_s - (pi* - pi-source) dt, and T is re-inverted
  !> at constant E: gas mass is then conserved exactly (M_s + M_g + M_out of
  !> the energy budget), the gas left over leaves in the next step's solve
  !> with its own enthalpy, and no mass is ever priced at the absolute h_g --
  !> without it the mismatch moved T by 0.6 K on the 1000 K adiabatic 0-D gate
  !> at first order in dt (2026-10-09). Not under quasi-steady (no stored gas
  !> in E). Inside the parallel region, after the last RK stage.
  subroutine Darcy_Mass_Correction ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_darcy, obj_pyrolysis
    use FUSS_Lib_Pyrolysis,  only: Energy_Of2, Temperature_Of2, Solid_State, Reaction_Rate
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b, i, j, k, m, r
    real(R8) :: dm, E, T_new, tau, eps, G, Gold, deficit, avail
    logical  :: all_r

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      !$omp do collapse(3) schedule(static) private(i,j,k,m,r,dm,E,T_new,tau,eps,G,Gold,deficit,avail,all_r)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        m  = nint( domain%blk(b)%matID(i,j,k) )
        dm = 0.0d0
        do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
          dm = dm + obj_pyrolysis%reac(r)%eps_rho_F * ( domain%blk(b)%chi(r,i,j,k) - domain%blk(b)%chiO(r,i,j,k) )
        enddo
        dm = dm - ( domain%blk(b)%pi_pg(i,j,k) - obj_darcy%pi_source ) * domain%blk(b)%dtlocal(i,j,k)
        if ( dm == 0.0d0 ) cycle
        if ( obj_darcy%quasi_steady ) then
          ! no stored gas: the difference waits in the carry, eta_st - eta_prev
          ! (Solve_Pressure set eta_st = eta_prev + carry/p), for the next solve
          domain%blk(b)%eta_st(i,j,k) = domain%blk(b)%eta_st(i,j,k) + dm / domain%blk(b)%p_pg(i,j,k)
          cycle
        endif
        E = Energy_Of2( m, domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), domain%blk(b)%p_pg(i,j,k), &
                        domain%blk(b)%eta_st(i,j,k) * domain%blk(b)%p_pg(i,j,k) )
        ! The predicted production can exceed the realised one by a whole step
        ! (threshold kinetics: a cell that drops back below T_thr during the
        ! stages). When the correction would empty the stored gas (a negative
        ! pressure for the next solve otherwise; Baum & Atreya case with a step
        ! wall, 2026-10-13), the solid supplies the deficit, up to a tenth of
        ! the start-of-step gas: the reactions that ran this step (or can run at
        ! T) advance further, in proportion to what each has left. E is held, so
        ! T comes out of the re-inversion below with the extra endotherm, and
        ! the gas mass stays exact. A run whose stored gas stays positive never
        ! enters this branch.
        Gold = domain%blk(b)%eta_stO(i,j,k) * domain%blk(b)%p_pgO(i,j,k)
        G    = domain%blk(b)%eta_st(i,j,k) * domain%blk(b)%p_pg(i,j,k)
        if ( G + dm <= 0.0d0 ) then
          deficit = 0.1d0 * Gold - ( G + dm )
          ! the reactions that ran this step (or can run at T) first; if they
          ! cannot cover it (the prediction crossed a threshold the stages did
          ! not reach: instant kinetics, a whole cell predicted), every reaction
          ! with mass left -- the cell decomposes now, as the solve assumed
          all_r = .false.
          avail = 0.0d0
          do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
            if ( domain%blk(b)%chi(r,i,j,k) > domain%blk(b)%chiO(r,i,j,k) .or. &
                 Reaction_Rate( m, r, domain%blk(b)%T(i,j,k) ) > 0.0d0 ) &
              avail = avail + obj_pyrolysis%reac(r)%eps_rho_F * ( 1.0d0 - domain%blk(b)%chi(r,i,j,k) )
          enddo
          if ( avail < deficit ) then
            all_r = .true.
            avail = 0.0d0
            do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
              avail = avail + obj_pyrolysis%reac(r)%eps_rho_F * ( 1.0d0 - domain%blk(b)%chi(r,i,j,k) )
            enddo
          endif
          ! as much as the solid has, up to the target; only a solid that
          ! cannot keep the stored gas positive stops the run
          deficit = min( deficit, avail )
          if ( G + dm + deficit <= 0.0d0 ) then
            write(*,'(A,4I6)') ' [ERROR] Darcy: the end-of-step gas correction empties the stored gas at b, i, j, k =', b, i, j, k
            write(*,'(A,ES14.6,A,ES14.6,A,ES14.6)') '         stored G = ', G, ' kg/m3, realised - predicted = ', dm, &
              ' kg/m3, decomposable solid left = ', avail
            write(*,'(A)') '         (reduce the time step, or smooth the kinetics near its threshold temperature)'
            error stop 3
          endif
          do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
            if ( all_r .or. domain%blk(b)%chi(r,i,j,k) > domain%blk(b)%chiO(r,i,j,k) .or. &
                 Reaction_Rate( m, r, domain%blk(b)%T(i,j,k) ) > 0.0d0 ) &
              domain%blk(b)%chi(r,i,j,k) = domain%blk(b)%chi(r,i,j,k) + ( deficit / avail ) * ( 1.0d0 - domain%blk(b)%chi(r,i,j,k) )
          enddo
          dm = dm + deficit
        endif
        domain%blk(b)%eta_st(i,j,k) = ( domain%blk(b)%eta_st(i,j,k) * domain%blk(b)%p_pg(i,j,k) + dm ) / domain%blk(b)%p_pg(i,j,k)
        call Temperature_Of2( m, E, domain%blk(b)%chi(:,i,j,k), domain%blk(b)%p_pg(i,j,k), &
                              domain%blk(b)%eta_st(i,j,k) * domain%blk(b)%p_pg(i,j,k), domain%blk(b)%T(i,j,k), T_new, b, i, j, k )
        domain%blk(b)%T(i,j,k) = T_new
        call Solid_State( m, T_new, domain%blk(b)%chi(:,i,j,k), domain%blk(b)%rho_s(i,j,k), &
                          domain%blk(b)%cp_s(i,j,k), domain%blk(b)%k_s(i,j,k), tau, eps )
      enddo ; enddo ; enddo
      !$omp end do
    enddo

  end subroutine Darcy_Mass_Correction


  !> Set-up banner lines (root).
  subroutine Darcy_Print_Setup ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_darcy
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    integer :: b, k
    character(len=64) :: codes
    character(len=8)  :: tok

    codes = ''
    do k = 1, obj_darcy%n_open
      write(tok,'(I0)') obj_darcy%open_codes(k)
      codes = trim(codes)//' '//trim(tok)
    enddo
    write(*,'(A)') ' Darcy pressure equation (type 2)'
    write(*,'(A,T35,A)') '   Solver', trim(obj_darcy%solver)
    do b = 1, domain%nb
      if ( obj_darcy%blk_dir(b) > 0 ) then
        write(*,'(A,I0,T35,A,I0)') '   Block ', b, 'line-implicit along direction ', obj_darcy%blk_dir(b)
      else
        write(*,'(A,I0,T35,A)') '   Block ', b, 'damped Jacobi sweep (pseudo)'
      endif
    enddo
    write(*,'(A,T35,A)') '   Open BC codes', trim(adjustl(codes))
    write(*,'(A,T35,ES12.4,A,ES12.4)') '   p-init, p-wall [Pa]', obj_darcy%p_init, '  ', obj_darcy%p_wall
    write(*,'(A,T35,ES10.2,A,I0)') '   tol, max-iter', obj_darcy%tol, '  ', obj_darcy%max_iter
    write(*,'(A,T35,L1,A,L1)') '   quasi-steady, klinkenberg', obj_darcy%quasi_steady, '  ', obj_darcy%klinkenberg
    if ( obj_darcy%pi_source > 0.0d0 ) &
      write(*,'(A,T35,ES12.4)') '   pi-source [kg/m3/s] (TEST)', obj_darcy%pi_source

  end subroutine Darcy_Print_Setup


  ! ---------------------------------------------------------------------------
  !  ghost layer
  ! ---------------------------------------------------------------------------

  !> The p instance of Fill_Ghost_Cell (plan 07 §3.4): connections copy the
  !> neighbour's cells into the ghost layer and bc%pg(1:2) (block-Jacobi
  !> across blocks, the value of the previous sweep); every other face mirrors
  !> (zero normal gradient, the impermeable default; open faces take p_w in
  !> the flux, not in the ghost). pg(3:6) is the neighbour's tangential
  !> stencil, as Tg(3:6). Chimera faces are refused under type 2.
  subroutine Fill_Ghost_p ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Global_m,     only: gc
    use FUSS_Parameters_m, only: guide
    use FUSS_Mod_MPI,      only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: ii, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, g, Ig, Jg, Kg
    integer :: ia, ja, ka, ib, jb, kb, icc, jcc, kcc, id, jd, kd, d11s, d12s, d21s, d22s

    !$omp do schedule(dynamic) private(ii,i,Bm,Im,Jm,Km,Fm,Bs,Is,Js,Ks,Fs,g,Ig,Jg,Kg)
    do ii = 1, domain%n_local_bc
      i  = domain%local_bc_idx(ii)
      Bm = domain%bc(i)%b ; Im = domain%bc(i)%i ; Jm = domain%bc(i)%j ; Km = domain%bc(i)%k ; Fm = domain%bc(i)%f
      select case ( domain%bc(i)%type )
        case (101)
          Bs = domain%bc(i)%bs
          if ( .not. is_local_block(Bs) ) cycle   ! inter-rank: plan 11
          Is = domain%bc(i)%is ; Js = domain%bc(i)%js ; Ks = domain%bc(i)%ks ; Fs = domain%bc(i)%fs
          do g = 1, gc
            Ig = Im - guide(Fm,1)*g ; Jg = Jm - guide(Fm,2)*g ; Kg = Km - guide(Fm,3)*g
            domain%blk(Bm)%p_pg(Ig,Jg,Kg) = domain%blk(Bs)%p_pg( Is + guide(Fs,1)*(g-1), &
                                                                 Js + guide(Fs,2)*(g-1), &
                                                                 Ks + guide(Fs,3)*(g-1) )
            domain%bc(i)%pg(g) = domain%blk(Bm)%p_pg(Ig,Jg,Kg)
          enddo
        case (102)
          write(*,'(A)') ' [ERROR] Darcy: chimera (102) faces are refused under gas-closure = type2'
          error stop 1
        case default
          do g = 1, gc
            select case (Fm)
              case(1); domain%blk(Bm)%p_pg(Im-g,Jm,Km) = domain%blk(Bm)%p_pg(Im+g-1,Jm,Km)
              case(2); domain%blk(Bm)%p_pg(Im+g,Jm,Km) = domain%blk(Bm)%p_pg(Im-g+1,Jm,Km)
              case(3); domain%blk(Bm)%p_pg(Im,Jm-g,Km) = domain%blk(Bm)%p_pg(Im,Jm+g-1,Km)
              case(4); domain%blk(Bm)%p_pg(Im,Jm+g,Km) = domain%blk(Bm)%p_pg(Im,Jm-g+1,Km)
              case(5); domain%blk(Bm)%p_pg(Im,Jm,Km-g) = domain%blk(Bm)%p_pg(Im,Jm,Km+g-1)
              case(6); domain%blk(Bm)%p_pg(Im,Jm,Km+g) = domain%blk(Bm)%p_pg(Im,Jm,Km-g+1)
            end select
          enddo
      end select
    enddo

    ! pg(3:6): the source block's tangential stencil (Fill_BC_Ghost_Connection)
    !$omp do schedule(dynamic) private(ii,i,Bs,Is,Js,Ks,Fs,d11s,d12s,d21s,d22s,ia,ja,ka,ib,jb,kb,icc,jcc,kcc,id,jd,kd)
    do ii = 1, domain%n_local_bs
      i  = domain%local_bs_idx(ii)
      Bs = domain%bc(i)%bs ; Is = domain%bc(i)%is ; Js = domain%bc(i)%js ; Ks = domain%bc(i)%ks ; Fs = domain%bc(i)%fs
      d11s = domain%bc(i)%d11 ; d12s = domain%bc(i)%d12 ; d21s = domain%bc(i)%d21 ; d22s = domain%bc(i)%d22
      select case (Fs)
        case (1:2)
          ia = Is ; ja = Js - d11s ; ka = Ks - d12s ; ib = Is ; jb = Js + d11s ; kb = Ks + d12s
          icc = Is ; jcc = Js - d21s ; kcc = Ks - d22s ; id = Is ; jd = Js + d21s ; kd = Ks + d22s
        case (3:4)
          ia = Is - d11s ; ja = Js ; ka = Ks - d12s ; ib = Is + d11s ; jb = Js ; kb = Ks + d12s
          icc = Is - d21s ; jcc = Js ; kcc = Ks - d22s ; id = Is + d21s ; jd = Js ; kd = Ks + d22s
        case default
          ia = Is - d11s ; ja = Js - d12s ; ka = Ks ; ib = Is + d11s ; jb = Js + d12s ; kb = Ks
          icc = Is - d21s ; jcc = Js - d22s ; kcc = Ks ; id = Is + d21s ; jd = Js + d22s ; kd = Ks
      end select
      domain%bc(i)%pg(3) = domain%blk(Bs)%p_pg(ia,ja,ka)
      domain%bc(i)%pg(4) = domain%blk(Bs)%p_pg(ib,jb,kb)
      domain%bc(i)%pg(5) = domain%blk(Bs)%p_pg(icc,jcc,kcc)
      domain%bc(i)%pg(6) = domain%blk(Bs)%p_pg(id,jd,kd)
    enddo

  end subroutine Fill_Ghost_p


  ! ---------------------------------------------------------------------------
  !  the solve
  ! ---------------------------------------------------------------------------

  !> p^(n+1) from the lagged coefficients of the start-of-step state (TO,
  !> chiO, p_pgO), then the derived fields: ghosts, cell velocity, open-face
  !> mass fluxes, eta^(n-1) <- eta^n for the next step.
  subroutine Solve_Pressure ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_darcy, obj_pyrolysis, obj_io, obj_sim_param, obj_ablation
    use FUSS_Global_m,       only: FUSS_phase_prefix
    use FUSS_Lib_Pyrolysis,  only: Darcy_Cell_Coeffs, Reaction_Rate, Advance_Reaction
    use FUSS_Mod_MPI,        only: is_local_block, mpi_is_root
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b, i, j, k, m, r, u, ios, bw, iw, jw, kw
    real(R8) :: eta, Gam(3,3), eps, rho_g, h, pi, dt, kap, rel, rw, Tstar, cq
    character(len=256) :: file

    ! 1. lagged coefficients and the gas production of this step (plan 07
    !    §3.3): pi = sum_r eps_rho_F_r [Advance_Reaction(chi^n, kappa(T*), dt) - chi^n]/dt
    !    + pi-source, eta and Gamma from (p^n, T^n, chi^n). T* is the
    !    mid-step temperature extrapolated from the previous step,
    !    T^n + (T^n - T^(n-1))/2 (T^n itself on the first step of a new run):
    !    with kappa at T^n the gas the solve expels lags the gas the RK
    !    stages release from the solid by O(dt), which the exact 0-D type-2
    !    gates measured as a FIRST-order error (1.7 K at ATC1's 2 ms step,
    !    halving with dt; 2026-10-08) -- PATO's own lag. The extrapolation
    !    makes the expelled mass second-order consistent at no extra cost.
    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      !$omp do collapse(3) schedule(static) private(i,j,k,m,r,eta,Gam,eps,rho_g,h,pi,dt,kap,Tstar,cq)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        m  = nint( domain%blk(b)%matID(i,j,k) )
        dt = domain%blk(b)%dtlocal(i,j,k)
        pi = obj_darcy%pi_source
        Tstar = domain%blk(b)%TO(i,j,k) + 0.5d0 * ( domain%blk(b)%TO(i,j,k) - domain%blk(b)%T_prev(i,j,k) )
        do r = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
          kap = Reaction_Rate( m, r, Tstar )
          pi  = pi + obj_pyrolysis%reac(r)%eps_rho_F * &
                ( Advance_Reaction( domain%blk(b)%chiO(r,i,j,k), kap, obj_pyrolysis%reac(r)%m, dt ) &
                  - domain%blk(b)%chiO(r,i,j,k) ) / dt
        enddo
        ! Quasi-steady: no stored gas absorbs the difference between the gas
        ! predicted here and the gas the RK stages release, so
        ! Darcy_Mass_Correction carries it, as (eta_st - eta_prev) p (restart
        ! columns all three), into this source: the gas mass is exact with a
        ! one-step lag. A negative carry (more expelled than released) is
        ! taken only down to a zero source; the rest waits for the next step.
        if ( obj_darcy%quasi_steady ) then
          cq = ( domain%blk(b)%eta_st(i,j,k) - domain%blk(b)%eta_prev(i,j,k) ) * domain%blk(b)%p_pgO(i,j,k)
          wk(b)%cq(i,j,k) = cq
          if ( cq /= 0.0d0 ) pi = obj_darcy%pi_source + max( 0.0d0, ( pi - obj_darcy%pi_source ) + cq / dt )
        endif
        call Darcy_Cell_Coeffs( m, domain%blk(b)%p_pgO(i,j,k), domain%blk(b)%TO(i,j,k), domain%blk(b)%chiO(:,i,j,k), &
                                [b, i, j, k], eta, Gam, eps, rho_g, h )
        domain%blk(b)%eta_pg(i,j,k)   = eta
        domain%blk(b)%gam_pg(:,i,j,k) = Pack9( Gam )
        domain%blk(b)%pi_pg(i,j,k)    = pi
      enddo ; enddo ; enddo
      !$omp end do
    enddo

    ! 2. sweeps: residual of the current iterate (ghosts refreshed first),
    !    stop test, update into p_new, copy. nsweep counts applied updates;
    !    at least one is applied per step, so a line-implicit block solves
    !    its lagged system exactly every step (one sweep on a 1-D line) and
    !    the stop test judges the UPDATED iterate.
    !$omp single
    nsweep = 0 ; done = .false. ; dp_max = 0.0d0
    !$omp end single
    do
      call Fill_Ghost_p( domain )
      !$omp single
      resid_max = 0.0d0 ; fmax = 0.0d0 ; smax = 0.0d0 ; floor_max = 0.0d0
      !$omp end single
      call Assemble( domain )
      !$omp single
      rref_max = max( fmax, smax )
      ! The residual cannot be evaluated below the round-off of its largest
      ! terms: each face flux C (p_R - p_L) carries an error of order
      ! epsilon C p, so a converged line solve on a nearly uniform field at
      ! a large permeability (gate 4.1(d): K = 1e-7 m2, fluxes vanishing as
      ! the source dies) shows |r| ~ epsilon sum_f C_f |p| whatever the sweep
      ! count. floor_max is that bound (64 epsilon times the largest
      ! diagonal-times-pressure term: 16 left the 0-D case at the cap).
      done = ( nsweep >= 1 .and. ( resid_max <= obj_darcy%tol * rref_max .or. resid_max <= floor_max ) ) &
             .or. ( nsweep >= obj_darcy%max_iter )
      !$omp end single
      if ( done ) exit
      call Update( domain )
      !$omp single
      nsweep = nsweep + 1
      !$omp end single
    enddo

    ! 3. failure, diagnostics (one thread; the arrays are shared)
    !$omp single
    if ( resid_max > obj_darcy%tol * rref_max .and. resid_max > floor_max ) then
      rw = -1.0d0 ; bw = 0 ; iw = 0 ; jw = 0 ; kw = 0
      do b = 1, domain%nb
        if ( .not. is_local_block(b) ) cycle
        do k = 1, domain%blk(b)%dim(3) ; do j = 1, domain%blk(b)%dim(2) ; do i = 1, domain%blk(b)%dim(1)
          if ( abs( wk(b)%res(i,j,k) ) > rw ) then
            rw = abs( wk(b)%res(i,j,k) ) ; bw = b ; iw = i ; jw = j ; kw = k
          endif
        enddo ; enddo ; enddo
      enddo
      write(*,'(A)') ' [ERROR] Darcy: the pressure sweeps did not converge within max-iter'
      write(*,'(A,I0,A,ES12.4,A,ES12.4)') '         sweeps ', nsweep, ', max |r|/R_ref = ', rw / max( rref_max, tiny(1.0d0) ), &
        ', tol = ', obj_darcy%tol
      write(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A)') '         max |r| = ', rw, ' kg/s, R_ref = ', rref_max, &
        ' kg/s, round-off floor = ', floor_max, ' kg/s'
      write(*,'(A,4I6)') '         worst cell b, i, j, k =', bw, iw, jw, kw
      error stop 3
    endif
    if ( rref_max > 0.0d0 ) then
      rel = resid_max / rref_max
    else
      rel = resid_max
    endif
    ! (the first pass of a new run's doubled first step is not logged: its
    ! pressure is discarded, the second pass writes the step's row)
    if ( mpi_is_root .and. .not. first_step_pending ) then
      file = 'OUTPUT/'//trim(FUSS_phase_prefix)//'darcy.txt'
      if ( .not. diag_started ) then
        if ( obj_sim_param%newrun ) then
          open( newunit=u, file=trim(file), status='replace', action='write', iostat=ios )
          if ( ios == 0 ) write(u,'(A)') '# iter  time  block  sweeps  resid  dp_last   (resid = max|r|/R_ref of the converged iterate; dp_last = max|dp|/p-ref of the last sweep)'
        else
          open( newunit=u, file=trim(file), status='unknown', position='append', action='write', iostat=ios )
        endif
        diag_started = .true.
      else
        open( newunit=u, file=trim(file), status='old', position='append', action='write', iostat=ios )
      endif
      if ( ios == 0 ) then
        do b = 1, domain%nb
          if ( .not. is_local_block(b) ) cycle
          if ( rref_max > 0.0d0 ) then
            write(u,'(I8,ES24.16,I6,I10,2ES16.8)') domain%iter, domain%time, b, nsweep, resb(b) / rref_max, dp_max
          else
            write(u,'(I8,ES24.16,I6,I10,2ES16.8)') domain%iter, domain%time, b, nsweep, resb(b), dp_max
          endif
        enddo
        close(u)
      endif
      if ( mod( domain%iter, obj_io%shell_diter ) == 0 ) &
        write(*,'(A,I0,A,ES10.3,A,ES10.3)') '   [Darcy] sweeps ', nsweep, '  resid ', rel, '  dp_last ', dp_max
    endif
    !$omp end single

    ! 4. derived fields of p^(n+1): ghosts (the loop left them current),
    !    cell velocity and flux, open-face mass fluxes; the stored gas of the
    !    step, G^(n+1) = eta* p^(n+1) (eta_st = eta*, the same expression as
    !    Assemble's; under quasi-steady E carries no gas, eta_st = eta^n), and
    !    eta^(n-1) for the next step
    call Darcy_Update_Velocity( domain )
    call Open_Face_Fluxes( domain )
    if ( obj_ablation%dump_interface_flux ) call Connection_Mass_Fluxes( domain )
    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      !$omp do collapse(3) schedule(static) private(i,j,k)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        if ( obj_darcy%quasi_steady ) then
          domain%blk(b)%eta_st(i,j,k) = domain%blk(b)%eta_pg(i,j,k) + wk(b)%cq(i,j,k) / domain%blk(b)%p_pg(i,j,k)
        else
          domain%blk(b)%eta_st(i,j,k) = Eta_Star( domain%blk(b)%eta_pg(i,j,k), domain%blk(b)%eta_prev(i,j,k) )
        endif
        domain%blk(b)%eta_prev(i,j,k) = domain%blk(b)%eta_pg(i,j,k)
      enddo ; enddo ; enddo
      !$omp end do
    enddo

  end subroutine Solve_Pressure


  !> One pass over the current iterate: the cell residual, the row
  !> coefficients of the line solve and the Jacobi diagonal, the reference
  !> magnitudes (max |face flux|, max(|pi V|, |storage|)).
  subroutine Assemble ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_darcy
    use FUSS_Lib_Pyrolysis,  only: Darcy_Cell_Coeffs
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b, i, j, k, ni, nj, nk, bd, fc, lower, upper, ii, Im, Jm, Km, Fm, Dir, m
    real(R8) :: C, E, F, Cw, Ew, s, A, q, inflow, Gw(3,3), eta_w, eps_w, rho_w, h_w, V, dt, ast, st, r_c, etas

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      ni = domain%blk(b)%dim(1) ; nj = domain%blk(b)%dim(2) ; nk = domain%blk(b)%dim(3)
      bd = obj_darcy%blk_dir(b)

      !$omp do collapse(3) schedule(static) private(i,j,k)
      do k = 1, nk ; do j = 1, nj ; do i = 1, ni
        wk(b)%res(i,j,k) = 0.0d0 ; wk(b)%dgl(i,j,k) = 0.0d0 ; wk(b)%dga(i,j,k) = 0.0d0
        wk(b)%rhs(i,j,k) = 0.0d0 ; wk(b)%lo(i,j,k)  = 0.0d0 ; wk(b)%hi(i,j,k)  = 0.0d0
        wk(b)%dgx(i,j,k) = 0.0d0 ; wk(b)%cx(:,i,j,k) = 0.0d0
      enddo ; enddo ; enddo
      !$omp end do

      ! interior faces: each thread owns whole lines of the face direction,
      ! so the two-cell scatter cannot race (Mod_Fluxes' collapse(2) pattern)
      !$omp do collapse(2) schedule(static) private(i,j,k,C,E,F) reduction(max: fmax)
      do k = 1, nk ; do j = 1, nj
        do i = 1, ni - 1
          call Face_Terms( domain%blk(b), domain%blk(b)%gam_pg, 1, i, j, k, C, E )
          F = C * ( domain%blk(b)%p_pg(i+1,j,k) - domain%blk(b)%p_pg(i,j,k) ) + E
          call Scatter( wk(b), bd, 1, i, j, k, i+1, j, k, C, E, F, &
                        domain%blk(b)%p_pg(i,j,k), domain%blk(b)%p_pg(i+1,j,k) )
          fmax = max( fmax, abs(F) )
        enddo
      enddo ; enddo
      !$omp end do

      !$omp do collapse(2) schedule(static) private(i,j,k,C,E,F) reduction(max: fmax)
      do k = 1, nk ; do i = 1, ni
        do j = 1, nj - 1
          call Face_Terms( domain%blk(b), domain%blk(b)%gam_pg, 2, i, j, k, C, E )
          F = C * ( domain%blk(b)%p_pg(i,j+1,k) - domain%blk(b)%p_pg(i,j,k) ) + E
          call Scatter( wk(b), bd, 2, i, j, k, i, j+1, k, C, E, F, &
                        domain%blk(b)%p_pg(i,j,k), domain%blk(b)%p_pg(i,j+1,k) )
          fmax = max( fmax, abs(F) )
        enddo
      enddo ; enddo
      !$omp end do

      !$omp do collapse(2) schedule(static) private(i,j,k,C,E,F) reduction(max: fmax)
      do j = 1, nj ; do i = 1, ni
        do k = 1, nk - 1
          call Face_Terms( domain%blk(b), domain%blk(b)%gam_pg, 3, i, j, k, C, E )
          F = C * ( domain%blk(b)%p_pg(i,j,k+1) - domain%blk(b)%p_pg(i,j,k) ) + E
          call Scatter( wk(b), bd, 3, i, j, k, i, j, k+1, C, E, F, &
                        domain%blk(b)%p_pg(i,j,k), domain%blk(b)%p_pg(i,j,k+1) )
          fmax = max( fmax, abs(F) )
        enddo
      enddo ; enddo
      !$omp end do
    enddo

    ! boundary faces, block by block and face by face as BC_Fluxes does: a
    ! corner cell has entries on several faces, never twice on one
    upper = 0
    do b = 1, domain%nb
      do fc = 1, 6
        lower = upper + 1
        upper = upper + domain%n_bf(b,fc)
        if ( .not. is_local_block(b) ) cycle
        bd = obj_darcy%blk_dir(b)
        !$omp do schedule(dynamic) private(ii,Im,Jm,Km,Fm,Dir,m,Cw,Ew,s,A,q,inflow,Gw,eta_w,eps_w,rho_w,h_w) reduction(max: fmax)
        do ii = lower, upper
          Im = domain%bc(ii)%i ; Jm = domain%bc(ii)%j ; Km = domain%bc(ii)%k ; Fm = domain%bc(ii)%f
          Dir = ( Fm + 1 ) / 2
          s = real( 1 - 2 * mod(Fm,2), R8 )
          if ( Darcy_Is_Open( domain%bc(ii)%type ) ) then
            q = domain%bc(ii)%pw
            m = nint( domain%blk(b)%matID(Im,Jm,Km) )
            call Darcy_Cell_Coeffs( m, q, domain%blk(b)%TO(Im,Jm,Km), domain%blk(b)%chiO(:,Im,Jm,Km), &
                                    [b, Im, Jm, Km], eta_w, Gw, eps_w, rho_w, h_w )
            call Wall_Terms( domain%blk(b), 0.5d0 * ( Tensor9( domain%blk(b)%gam_pg(:,Im,Jm,Km) ) + Gw ), &
                             Im, Jm, Km, Fm, .true., Cw, Ew, A )
          elseif ( domain%bc(ii)%type == 101 ) then
            ! block-Jacobi connection (plan 08 section 2.9 (3)): the ghost
            ! value of the previous sweep and the interior face kernel on both
            ! blocks' data -- face metric (M + Mg)/2, face tensor from the two
            ! cells' own lagged coefficients, both tangential stencils -- so
            ! the two sides evaluate one flux and the gas mass is conserved
            q = domain%bc(ii)%pg(1)
            call Connection_Terms( domain%blk(b), domain%bc(ii), domain%blk(b)%gam_pg(:,Im,Jm,Km), &
                                   domain%blk(domain%bc(ii)%bs)%gam_pg(:,domain%bc(ii)%is,domain%bc(ii)%js, &
                                   domain%bc(ii)%ks), Cw, Ew, A )
          else
            cycle
          endif
          inflow = Cw * ( q - domain%blk(b)%p_pg(Im,Jm,Km) ) + s * Ew
          wk(b)%res(Im,Jm,Km) = wk(b)%res(Im,Jm,Km) + inflow
          wk(b)%dga(Im,Jm,Km) = wk(b)%dga(Im,Jm,Km) + Cw
          ! the face's conductance on the line diagonal whatever its direction
          ! (a face across the lines too, as Scatter's; the value beyond it
          ! from the previous sweep): line Jacobi with the full diagonal
          wk(b)%dgl(Im,Jm,Km) = wk(b)%dgl(Im,Jm,Km) + Cw
          wk(b)%rhs(Im,Jm,Km) = wk(b)%rhs(Im,Jm,Km) + Cw * q + s * Ew
          fmax = max( fmax, abs(inflow) )
        enddo
        !$omp end do
      enddo
    enddo

    ! cells: source and storage, the residual, the diagonals
    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      !$omp single
      resid_blk = 0.0d0
      !$omp end single
      !$omp do collapse(3) schedule(static) private(i,j,k,V,dt,ast,st,r_c,etas) reduction(max: resid_max, resid_blk, smax, floor_max)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        V  = domain%blk(b)%vol(i,j,k)
        dt = domain%blk(b)%dtlocal(i,j,k)
        if ( obj_darcy%quasi_steady ) then
          ast = 0.0d0 ; st = 0.0d0
        else
          ! Storage of the stored gas mass G (plan 07, the scheme of 2026-10-09):
          ! G^(n+1) - G^n with G^n = eta_stO p^n the mass the previous solve left
          ! stored and G^(n+1) = eta* p the mass this one stores, eta* = 2 eta^n -
          ! eta^(n-1) the end-of-step coefficient extrapolated from the two
          ! lagged equation-of-state values (eta^n on the first step). The sum
          ! over steps telescopes, so the gas mass is conserved exactly, and E
          ! carries this same G (FUSS_Lib_Pyrolysis, type-2 section), so no gas
          ! mass is ever priced at the absolute h_g; G departs from the equation
          ! of state only by the O(dt^2) error of the extrapolation, re-anchored
          ! every step. Written as eta* (p - p^n) + (eta* - eta_stO) p^n: on a
          ! stationary sealed field both parts are exactly zero, so the residual
          ! is exactly zero as the stop test demands when no flux and no source
          ! exist; the literal eta* p - G^n is not (the compiler contracts it into
          ! a fused multiply-add and the rounding error of one product survives,
          ! seen: 1.4e-21 kg/s on gate 4.1(a)). The reference magnitude counts
          ! both terms whole.
          etas = Eta_Star( domain%blk(b)%eta_pg(i,j,k), domain%blk(b)%eta_prev(i,j,k) )
          ast  = etas * ( V / dt )
          st   = ( etas * ( domain%blk(b)%p_pg(i,j,k) - domain%blk(b)%p_pgO(i,j,k) ) &
                 + ( etas - domain%blk(b)%eta_stO(i,j,k) ) * domain%blk(b)%p_pgO(i,j,k) ) * ( V / dt )
          wk(b)%rhs(i,j,k) = wk(b)%rhs(i,j,k) + domain%blk(b)%eta_stO(i,j,k) * domain%blk(b)%p_pgO(i,j,k) * ( V / dt )
          smax = max( smax, abs( etas * domain%blk(b)%p_pg(i,j,k) * ( V / dt ) ), &
                            abs( domain%blk(b)%eta_stO(i,j,k) * domain%blk(b)%p_pgO(i,j,k) * ( V / dt ) ) )
        endif
        r_c = wk(b)%res(i,j,k) + domain%blk(b)%pi_pg(i,j,k) * V - st
        wk(b)%res(i,j,k) = r_c
        wk(b)%dgl(i,j,k) = wk(b)%dgl(i,j,k) + ast
        wk(b)%dga(i,j,k) = wk(b)%dga(i,j,k) + ast
        wk(b)%rhs(i,j,k) = wk(b)%rhs(i,j,k) + domain%blk(b)%pi_pg(i,j,k) * V
        resid_max = max( resid_max, abs(r_c) )
        resid_blk = max( resid_blk, abs(r_c) )
        smax = max( smax, abs( domain%blk(b)%pi_pg(i,j,k) * V ) )
        floor_max = max( floor_max, 64.0d0 * epsilon(1.0d0) * wk(b)%dga(i,j,k) * abs( domain%blk(b)%p_pg(i,j,k) ) )
      enddo ; enddo ; enddo
      !$omp end do
      !$omp single
      resb(b) = resid_blk
      !$omp end single
    enddo

  end subroutine Assemble


  !> One update of every local block from the assembled rows: Thomas along
  !> the lines of a line-implicit block, damped Jacobi elsewhere; then the
  !> copy p_new -> p_pg with the per-sweep change.
  subroutine Update ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_darcy
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b, i, j, k, ni, nj, nk, bd
    real(R8) :: om

    om = obj_darcy%pseudo_vnn
    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      ni = domain%blk(b)%dim(1) ; nj = domain%blk(b)%dim(2) ; nk = domain%blk(b)%dim(3)
      bd = obj_darcy%blk_dir(b)
      select case ( bd )
        case (1)
          !$omp do collapse(2) schedule(static) private(j,k)
          do k = 1, nk ; do j = 1, nj
            call Thomas_Line( domain%blk(b), wk(b), 1, j, k )
          enddo ; enddo
          !$omp end do
          !$omp single
          call Line_Uniform_Correction( domain%blk(b), wk(b), 1 )
          !$omp end single
        case (2)
          !$omp do collapse(2) schedule(static) private(i,k)
          do k = 1, nk ; do i = 1, ni
            call Thomas_Line( domain%blk(b), wk(b), 2, i, k )
          enddo ; enddo
          !$omp end do
          !$omp single
          call Line_Uniform_Correction( domain%blk(b), wk(b), 2 )
          !$omp end single
        case (3)
          !$omp do collapse(2) schedule(static) private(i,j)
          do j = 1, nj ; do i = 1, ni
            call Thomas_Line( domain%blk(b), wk(b), 3, i, j )
          enddo ; enddo
          !$omp end do
          !$omp single
          call Line_Uniform_Correction( domain%blk(b), wk(b), 3 )
          !$omp end single
        case default
          !$omp do collapse(3) schedule(static) private(i,j,k)
          do k = 1, nk ; do j = 1, nj ; do i = 1, ni
            if ( wk(b)%dga(i,j,k) > 0.0d0 ) then
              domain%blk(b)%p_new(i,j,k) = domain%blk(b)%p_pg(i,j,k) + om * wk(b)%res(i,j,k) / wk(b)%dga(i,j,k)
            else
              domain%blk(b)%p_new(i,j,k) = domain%blk(b)%p_pg(i,j,k)
            endif
          enddo ; enddo ; enddo
          !$omp end do
      end select
      !$omp do collapse(3) schedule(static) private(i,j,k) reduction(max: dp_max)
      do k = 1, nk ; do j = 1, nj ; do i = 1, ni
        dp_max = max( dp_max, abs( domain%blk(b)%p_new(i,j,k) - domain%blk(b)%p_pg(i,j,k) ) / obj_darcy%p_ref )
        domain%blk(b)%p_pg(i,j,k) = domain%blk(b)%p_new(i,j,k)
      enddo ; enddo ; enddo
      !$omp end do
    enddo

  end subroutine Update


  !> Thomas solve of one line of direction d at transverse indices (t1, t2):
  !> lo p(c-1) + dgl p(c) + hi p(c+1) = rhs, result into p_new.
  subroutine Thomas_Line ( blk, w, d, t1, t2 )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    implicit none
    type(FUSS_block_type), intent(inout) :: blk
    type(darcy_work_t),    intent(in)    :: w
    integer, intent(in) :: d, t1, t2
    ! Local
    integer  :: n, s, ic(3)
    real(R8) :: cp( blk%dim(d) ), rp( blk%dim(d) ), den, xnext

    n = blk%dim(d)
    call Line_Index( d, 1, t1, t2, ic )
    den = w%dgl( ic(1), ic(2), ic(3) )
    if ( den == 0.0d0 ) call Singular( ic )
    cp(1) = w%hi( ic(1), ic(2), ic(3) ) / den
    rp(1) = w%rhs( ic(1), ic(2), ic(3) ) / den
    do s = 2, n
      call Line_Index( d, s, t1, t2, ic )
      den = w%dgl( ic(1), ic(2), ic(3) ) - w%lo( ic(1), ic(2), ic(3) ) * cp(s-1)
      if ( den == 0.0d0 ) call Singular( ic )
      cp(s) = w%hi( ic(1), ic(2), ic(3) ) / den
      rp(s) = ( w%rhs( ic(1), ic(2), ic(3) ) - w%lo( ic(1), ic(2), ic(3) ) * rp(s-1) ) / den
    enddo
    call Line_Index( d, n, t1, t2, ic )
    blk%p_new( ic(1), ic(2), ic(3) ) = rp(n)
    do s = n - 1, 1, -1
      call Line_Index( d, s + 1, t1, t2, ic )
      xnext = blk%p_new( ic(1), ic(2), ic(3) )
      call Line_Index( d, s, t1, t2, ic )
      blk%p_new( ic(1), ic(2), ic(3) ) = rp(s) - cp(s) * xnext
    enddo

  contains

    subroutine Singular ( ic )
      implicit none
      integer, intent(in) :: ic(3)
      write(*,'(A)') ' [ERROR] Darcy: singular line system (a sealed, quasi-steady line has no unique pressure)'
      write(*,'(A,3I6)') '         at i, j, k =', ic
      error stop 3
    end subroutine Singular

  end subroutine Thomas_Line


  !> The additive correction of the error uniform across the lines of a
  !> line-implicit block (Hutchinson and Raithby's additive correction, in the
  !> line direction). After the Thomas solves, the residual of the full system
  !> at p_new is the cross-line coupling of the change, r = sum_cross C (p_new -
  !> p_pg)(neighbour). The correction delta(s), one value per position s
  !> along the lines, solves the line-summed system, in which the cross terms
  !> cancel: sum_t [lo delta(s-1) + (dgl - dgx) delta(s) + hi delta(s+1)] =
  !> sum_t r(s,t). On lines that carry the same data (a one-dimensional case
  !> cut into identical lines) the corrected iterate is the exact line
  !> solution, so one sweep solves the step to round-off, as the lagged-flux
  !> splitting did before the cross conductance moved onto the diagonal; on a
  !> two-dimensional block it removes line Jacobi's slowest mode. A singular
  !> summed system (no line with a fixed pressure under quasi-steady) skips
  !> the correction. Serial (inside an omp single).
  subroutine Line_Uniform_Correction ( blk, w, d )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    implicit none
    type(FUSS_block_type), intent(inout) :: blk
    type(darcy_work_t),    intent(in)    :: w
    integer, intent(in) :: d
    ! Local
    integer  :: n, s, t1, t2, n1, n2, c1, c2, ic(3), jc(3), q
    real(R8) :: Ls(blk%dim(d)), Ds(blk%dim(d)), Hs(blk%dim(d)), Rs(blk%dim(d)), cp(blk%dim(d)), rp(blk%dim(d))
    real(R8) :: den, r, dl

    n = blk%dim(d)
    c1 = merge( 2, 1, d == 1 ) ; c2 = merge( 2, 3, d == 3 )
    n1 = blk%dim(c1) ; n2 = blk%dim(c2)
    if ( n1 * n2 == 1 ) return                      ! one line: the Thomas solve is exact
    Ls = 0.0d0 ; Ds = 0.0d0 ; Hs = 0.0d0 ; Rs = 0.0d0
    do s = 1, n
      do t2 = 1, n2 ; do t1 = 1, n1
        call Line_Index( d, s, t1, t2, ic )
        r = 0.0d0
        do q = 1, 4
          if ( w%cx(q,ic(1),ic(2),ic(3)) == 0.0d0 ) cycle
          jc = ic
          select case ( q )
            case (1) ; jc(c1) = ic(c1) - 1
            case (2) ; jc(c1) = ic(c1) + 1
            case (3) ; jc(c2) = ic(c2) - 1
            case (4) ; jc(c2) = ic(c2) + 1
          end select
          r = r + w%cx(q,ic(1),ic(2),ic(3)) * ( blk%p_new(jc(1),jc(2),jc(3)) - blk%p_pg(jc(1),jc(2),jc(3)) )
        enddo
        Rs(s) = Rs(s) + r
        Ds(s) = Ds(s) + ( w%dgl(ic(1),ic(2),ic(3)) - w%dgx(ic(1),ic(2),ic(3)) )
        Ls(s) = Ls(s) + w%lo(ic(1),ic(2),ic(3))
        Hs(s) = Hs(s) + w%hi(ic(1),ic(2),ic(3))
      enddo ; enddo
    enddo
    ! Thomas on the summed line
    if ( Ds(1) == 0.0d0 ) return
    cp(1) = Hs(1) / Ds(1) ; rp(1) = Rs(1) / Ds(1)
    do s = 2, n
      den = Ds(s) - Ls(s) * cp(s-1)
      if ( den == 0.0d0 ) return
      cp(s) = Hs(s) / den
      rp(s) = ( Rs(s) - Ls(s) * rp(s-1) ) / den
    enddo
    do s = n - 1, 1, -1
      rp(s) = rp(s) - cp(s) * rp(s+1)
    enddo
    do s = 1, n
      dl = rp(s)
      if ( dl == 0.0d0 ) cycle
      do t2 = 1, n2 ; do t1 = 1, n1
        call Line_Index( d, s, t1, t2, ic )
        blk%p_new(ic(1),ic(2),ic(3)) = blk%p_new(ic(1),ic(2),ic(3)) + dl
      enddo ; enddo
    enddo

  end subroutine Line_Uniform_Correction


  !> Cell (i,j,k) of position s along direction d at transverse (t1, t2):
  !> d = 1: (s, t1, t2); d = 2: (t1, s, t2); d = 3: (t1, t2, s).
  pure subroutine Line_Index ( d, s, t1, t2, ic )
    implicit none
    integer, intent(in)  :: d, s, t1, t2
    integer, intent(out) :: ic(3)
    select case ( d )
      case (1) ; ic = [ s, t1, t2 ]
      case (2) ; ic = [ t1, s, t2 ]
      case default ; ic = [ t1, t2, s ]
    end select
  end subroutine Line_Index


  ! ---------------------------------------------------------------------------
  !  face helpers
  ! ---------------------------------------------------------------------------

  !> (9) row-major <-> 3x3
  pure function Tensor9 ( v ) result ( T )
    implicit none
    real(R8), intent(in) :: v(9)
    real(R8) :: T(3,3)
    integer :: a, c
    do a = 1, 3
      do c = 1, 3
        T(a,c) = v( 3*(a-1) + c )
      enddo
    enddo
  end function Tensor9

  pure function Pack9 ( T ) result ( v )
    implicit none
    real(R8), intent(in) :: T(3,3)
    real(R8) :: v(9)
    integer :: a, c
    do a = 1, 3
      do c = 1, 3
        v( 3*(a-1) + c ) = T(a,c)
      enddo
    enddo
  end function Pack9


  !> C and E of the interior face of direction d between cell (i,j,k) and
  !> its +d neighbour, with the face tensor (gsrc(c) + gsrc(c+d))/2 and the
  !> 10-point stencil of Fluxes_blk.
  subroutine Face_Terms ( blk, gsrc, d, i, j, k, C, E )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    use FUSS_Lib_Diffusive_Tensor, only: Darcy_Face
    implicit none
    type(FUSS_block_type), intent(in) :: blk
    real(R8), intent(in)  :: gsrc(:,:,:,:)
    integer,  intent(in)  :: d, i, j, k
    real(R8), intent(out) :: C, E
    ! Local
    real(R8) :: Gf(3,3)

    select case ( d )
      case (1)
        Gf = 0.5d0 * ( Tensor9( gsrc(:,i,j,k) ) + Tensor9( gsrc(:,i+1,j,k) ) )
        call Darcy_Face( Gf, blk%dir(1)%f(i,j,k)%n, blk%dir(1)%f(i,j,k)%A, &
                         blk%p_pg(i,j-1,k), blk%p_pg(i,j+1,k), blk%p_pg(i+1,j-1,k), blk%p_pg(i+1,j+1,k), &
                         blk%p_pg(i,j,k-1), blk%p_pg(i,j,k+1), blk%p_pg(i+1,j,k-1), blk%p_pg(i+1,j,k+1), &
                         blk%M(i,j,k)%c, blk%M(i+1,j,k)%c, 1, 2, 3, C, E )
      case (2)
        Gf = 0.5d0 * ( Tensor9( gsrc(:,i,j,k) ) + Tensor9( gsrc(:,i,j+1,k) ) )
        call Darcy_Face( Gf, blk%dir(2)%f(i,j,k)%n, blk%dir(2)%f(i,j,k)%A, &
                         blk%p_pg(i-1,j,k), blk%p_pg(i+1,j,k), blk%p_pg(i-1,j+1,k), blk%p_pg(i+1,j+1,k), &
                         blk%p_pg(i,j,k-1), blk%p_pg(i,j,k+1), blk%p_pg(i,j+1,k-1), blk%p_pg(i,j+1,k+1), &
                         blk%M(i,j,k)%c, blk%M(i,j+1,k)%c, 2, 1, 3, C, E )
      case default
        Gf = 0.5d0 * ( Tensor9( gsrc(:,i,j,k) ) + Tensor9( gsrc(:,i,j,k+1) ) )
        call Darcy_Face( Gf, blk%dir(3)%f(i,j,k)%n, blk%dir(3)%f(i,j,k)%A, &
                         blk%p_pg(i-1,j,k), blk%p_pg(i+1,j,k), blk%p_pg(i-1,j,k+1), blk%p_pg(i+1,j,k+1), &
                         blk%p_pg(i,j-1,k), blk%p_pg(i,j+1,k), blk%p_pg(i,j-1,k+1), blk%p_pg(i,j+1,k+1), &
                         blk%M(i,j,k)%c, blk%M(i,j,k+1)%c, 3, 1, 2, C, E )
    end select

  end subroutine Face_Terms


  !> Cw, E and the area of the boundary face Fm of cell (Im,Jm,Km) with the
  !> face tensor Gf: half-cell distance for a wall (the Dirichlet closure of
  !> plan 07 §3.4), full-cell distance for a connection ghost. The tangential
  !> differences are the cell row's central differences (positive sense).
  subroutine Wall_Terms ( blk, Gf, Im, Jm, Km, Fm, half, Cw, E, A )
    use FUSS_Advanced_Types_m, only: FUSS_block_type
    use FUSS_Lib_BC_Fluxes,        only: Face_Index
    use FUSS_Lib_Diffusive_Tensor, only: Darcy_Boundary_Face
    implicit none
    type(FUSS_block_type), intent(in) :: blk
    real(R8), intent(in)  :: Gf(3,3)
    integer,  intent(in)  :: Im, Jm, Km, Fm
    logical,  intent(in)  :: half
    real(R8), intent(out) :: Cw, E, A
    ! Local
    integer  :: Dir, Fi, Fj, Fk
    real(R8) :: Gt(3), n(3)

    call Face_Index( Fm, Dir, Im, Jm, Km, Fi, Fj, Fk )
    n = blk%dir(Dir)%f(Fi,Fj,Fk)%n
    A = blk%dir(Dir)%f(Fi,Fj,Fk)%A
    Gt(1) = 0.5d0 * ( blk%p_pg(Im+1,Jm,Km) - blk%p_pg(Im-1,Jm,Km) )
    Gt(2) = 0.5d0 * ( blk%p_pg(Im,Jm+1,Km) - blk%p_pg(Im,Jm-1,Km) )
    Gt(3) = 0.5d0 * ( blk%p_pg(Im,Jm,Km+1) - blk%p_pg(Im,Jm,Km-1) )
    call Darcy_Boundary_Face( Gf, n, A, blk%M(Im,Jm,Km)%c, Gt, Dir, Cw, E )
    if ( .not. half ) Cw = 0.5d0 * Cw

  end subroutine Wall_Terms


  !> C and E of the block-connection (101) face of bc entry bcf (plan 08
  !> section 2.9 (3)): the interior kernel Darcy_Face with the face metric
  !> (M + Mg)/2 -- Mg the neighbour cell's, bcf%Mg(1) --, the face tensor
  !> (gloc + gnb)/2 of the two cells' own coefficients and the tangential pairs
  !> of both sides (this cell's row and the source block's, bcf%pg(3:6), as
  !> BC_Connection reads Tg(3:6)); the inflow into the cell is C (pg(1) - p) +
  !> s E, s = +1 on a high face. Both blocks evaluate the same expression, so
  !> the face conserves what it carries. Before plan 08 the row used this
  !> cell's metric alone at a full-cell distance and the neighbour's tensor
  !> re-evaluated at its ghost state: on unequal spacing the two sides'
  !> conductances differed (by the spacing ratio at R5's dx2/dx1 = 4).
  subroutine Connection_Terms ( blk, bcf, gloc, gnb, C, E, A )
    use FUSS_Advanced_Types_m,     only: FUSS_block_type, FUSS_bc_type
    use FUSS_Lib_BC_Fluxes,        only: Face_Index
    use FUSS_Lib_Diffusive_Tensor, only: Darcy_Face
    implicit none
    type(FUSS_block_type), intent(in) :: blk
    type(FUSS_bc_type),    intent(in) :: bcf
    real(R8), intent(in)  :: gloc(9), gnb(9)
    real(R8), intent(out) :: C, E, A
    ! Local
    integer  :: Dir, Fi, Fj, Fk, Im, Jm, Km, bb, cc
    real(R8) :: n(3), Gf(3,3), lb(2), lc(2)

    Im = bcf%i ; Jm = bcf%j ; Km = bcf%k
    call Face_Index( bcf%f, Dir, Im, Jm, Km, Fi, Fj, Fk )
    n = blk%dir(Dir)%f(Fi,Fj,Fk)%n
    A = blk%dir(Dir)%f(Fi,Fj,Fk)%A
    select case ( Dir )
      case (1)
        bb = 2 ; cc = 3
        lb = [ blk%p_pg(Im,Jm-1,Km), blk%p_pg(Im,Jm+1,Km) ]
        lc = [ blk%p_pg(Im,Jm,Km-1), blk%p_pg(Im,Jm,Km+1) ]
      case (2)
        bb = 1 ; cc = 3
        lb = [ blk%p_pg(Im-1,Jm,Km), blk%p_pg(Im+1,Jm,Km) ]
        lc = [ blk%p_pg(Im,Jm,Km-1), blk%p_pg(Im,Jm,Km+1) ]
      case default
        bb = 1 ; cc = 2
        lb = [ blk%p_pg(Im-1,Jm,Km), blk%p_pg(Im+1,Jm,Km) ]
        lc = [ blk%p_pg(Im,Jm-1,Km), blk%p_pg(Im,Jm+1,Km) ]
    end select
    Gf = 0.5d0 * ( Tensor9( gloc ) + Tensor9( gnb ) )
    call Darcy_Face( Gf, n, A, lb(1), lb(2), bcf%pg(3), bcf%pg(4), lc(1), lc(2), bcf%pg(5), bcf%pg(6), &
                     blk%M(Im,Jm,Km)%c, bcf%Mg(1)%c, Dir, bb, cc, C, E )

  end subroutine Connection_Terms


  !> TEST ONLY (dump-interface-flux, plan 08 R5): the gas mass flow [kg/s]
  !> into the cell through every 101 face from the converged p^(n+1) and its
  !> ghosts, the expression of Assemble's connection row.
  subroutine Connection_Mass_Fluxes ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: ii, b, Im, Jm, Km
    real(R8) :: Cw, Ew, A, s

    !$omp do schedule(dynamic) private(ii,b,Im,Jm,Km,Cw,Ew,A,s)
    do ii = 1, domain%nbound
      if ( domain%bc(ii)%type /= 101 ) cycle
      b = domain%bc(ii)%b
      if ( .not. is_local_block(b) ) cycle
      Im = domain%bc(ii)%i ; Jm = domain%bc(ii)%j ; Km = domain%bc(ii)%k
      s = real( 1 - 2 * mod(domain%bc(ii)%f,2), R8 )
      call Connection_Terms( domain%blk(b), domain%bc(ii), domain%blk(b)%gam_pg(:,Im,Jm,Km), &
                             domain%blk(domain%bc(ii)%bs)%gam_pg(:,domain%bc(ii)%is,domain%bc(ii)%js, &
                             domain%bc(ii)%ks), Cw, Ew, A )
      domain%bc(ii)%qm_if = Cw * ( domain%bc(ii)%pg(1) - domain%blk(b)%p_pg(Im,Jm,Km) ) + s * Ew
    enddo
    !$omp end do

  end subroutine Connection_Mass_Fluxes


  !> Scatter of an interior face (lower cell L, upper cell R; F = C (pR - pL)
  !> + E) into the rows: the residual and the Jacobi diagonal always; the line
  !> direction gets its conductance on the diagonal and off-diagonals with the
  !> tangential part E lagged. A face ACROSS the lines puts its conductance on
  !> the diagonal too and only the neighbour's part, C p_neighbour +- E of the
  !> current iterate, on the right-hand side: line Jacobi with the full
  !> diagonal, whose fixed point is the same system. The whole flux on the
  !> right-hand side (as first built) left the diagonal without the cross
  !> conductance, and the sweeps diverged once the storage V/dt fell below it
  !> -- a 2-D slab with its open face on top, at any time step (2026-10-13).
  subroutine Scatter ( w, bd, dd, iL, jL, kL, iR, jR, kR, C, E, F, pL, pR )
    implicit none
    type(darcy_work_t), intent(inout) :: w
    integer,  intent(in) :: bd, dd          ! the block's line direction (0: none) and the face's
    integer,  intent(in) :: iL, jL, kL, iR, jR, kR
    real(R8), intent(in) :: C, E, F, pL, pR
    ! Local
    logical :: implicit_line
    integer :: q

    implicit_line = ( bd == dd )

    w%res(iL,jL,kL) = w%res(iL,jL,kL) + F
    w%res(iR,jR,kR) = w%res(iR,jR,kR) - F
    w%dga(iL,jL,kL) = w%dga(iL,jL,kL) + C
    w%dga(iR,jR,kR) = w%dga(iR,jR,kR) + C
    if ( implicit_line ) then
      w%dgl(iL,jL,kL) = w%dgl(iL,jL,kL) + C
      w%dgl(iR,jR,kR) = w%dgl(iR,jR,kR) + C
      w%hi(iL,jL,kL)  = - C
      w%lo(iR,jR,kR)  = - C
      w%rhs(iL,jL,kL) = w%rhs(iL,jL,kL) + E
      w%rhs(iR,jR,kR) = w%rhs(iR,jR,kR) - E
    else
      w%dgl(iL,jL,kL) = w%dgl(iL,jL,kL) + C
      w%dgl(iR,jR,kR) = w%dgl(iR,jR,kR) + C
      w%rhs(iL,jL,kL) = w%rhs(iL,jL,kL) + ( C * pR + E )
      w%rhs(iR,jR,kR) = w%rhs(iR,jR,kR) + ( C * pL - E )
      if ( bd >= 1 .and. bd <= 3 ) then
        ! for the line-uniform correction (Line_Uniform_Correction): slot 1/2
        ! the first cross direction (lower, upper neighbour), 3/4 the second
        q = 3
        if ( dd == merge( 2, 1, bd == 1 ) ) q = 1
        w%dgx(iL,jL,kL) = w%dgx(iL,jL,kL) + C
        w%dgx(iR,jR,kR) = w%dgx(iR,jR,kR) + C
        w%cx(q+1,iL,jL,kL) = C          ! cell L's upper neighbour along dd is R
        w%cx(q,iR,jR,kR)   = C          ! cell R's lower neighbour along dd is L
      endif
    endif

  end subroutine Scatter


  ! ---------------------------------------------------------------------------
  !  derived fields
  ! ---------------------------------------------------------------------------

  !> Cell gas velocity u = -(K/(mu eps)) grad p = -Gamma grad p / (eta p)
  !> (eps rho_g = eta p) and the mass-flux magnitude |Gamma grad p| into
  !> mdot_pg, from central differences of p (ghosts current) and the cell
  !> metric. Inside or outside a parallel region.
  subroutine Darcy_Update_Velocity ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Global_m, only: ndir
    use FUSS_Mod_MPI,  only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b, i, j, k
    real(R8) :: G(3), grad(3), flux(3), den

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      !$omp do collapse(3) schedule(static) private(i,j,k,G,grad,flux,den)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        G(1) = 0.5d0 * ( domain%blk(b)%p_pg(i+1,j,k) - domain%blk(b)%p_pg(i-1,j,k) )
        G(2) = 0.5d0 * ( domain%blk(b)%p_pg(i,j+1,k) - domain%blk(b)%p_pg(i,j-1,k) )
        G(3) = 0.5d0 * ( domain%blk(b)%p_pg(i,j,k+1) - domain%blk(b)%p_pg(i,j,k-1) )
        if ( ndir == 2 ) G(3) = 0.0d0
        grad = matmul( G, domain%blk(b)%M(i,j,k)%c )
        flux = - matmul( Tensor9( domain%blk(b)%gam_pg(:,i,j,k) ), grad )
        domain%blk(b)%mdot_pg(i,j,k) = norm2( flux )
        den = domain%blk(b)%eta_pg(i,j,k) * domain%blk(b)%p_pg(i,j,k)
        if ( den > 0.0d0 ) then
          domain%blk(b)%u_pg(:,i,j,k) = flux / den
        else
          domain%blk(b)%u_pg(:,i,j,k) = 0.0d0
        endif
      enddo ; enddo ; enddo
      !$omp end do
    enddo

  end subroutine Darcy_Update_Velocity


  !> Mass flux leaving through every open face (bc%mdot_pg, kg/m2/s, outward
  !> positive) from the converged p: the quantity plan 08 §2.1 consumes.
  subroutine Open_Face_Fluxes ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Lib_Pyrolysis, only: Darcy_Cell_Coeffs
    use FUSS_Mod_MPI,       only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: ii, i, b, Im, Jm, Km, Fm, m
    real(R8) :: Cw, Ew, s, A, q, inflow, Gw(3,3), eta_w, eps_w, rho_w, h_w

    !$omp do schedule(dynamic) private(ii,i,b,Im,Jm,Km,Fm,m,Cw,Ew,s,A,q,inflow,Gw,eta_w,eps_w,rho_w,h_w)
    do ii = 1, domain%n_local_bc
      i = domain%local_bc_idx(ii)
      if ( .not. Darcy_Is_Open( domain%bc(i)%type ) ) cycle
      b = domain%bc(i)%b
      if ( .not. is_local_block(b) ) cycle
      Im = domain%bc(i)%i ; Jm = domain%bc(i)%j ; Km = domain%bc(i)%k ; Fm = domain%bc(i)%f
      s = real( 1 - 2 * mod(Fm,2), R8 )
      q = domain%bc(i)%pw
      m = nint( domain%blk(b)%matID(Im,Jm,Km) )
      call Darcy_Cell_Coeffs( m, q, domain%blk(b)%TO(Im,Jm,Km), domain%blk(b)%chiO(:,Im,Jm,Km), &
                              [b, Im, Jm, Km], eta_w, Gw, eps_w, rho_w, h_w )
      call Wall_Terms( domain%blk(b), 0.5d0 * ( Tensor9( domain%blk(b)%gam_pg(:,Im,Jm,Km) ) + Gw ), &
                       Im, Jm, Km, Fm, .true., Cw, Ew, A )
      inflow = Cw * ( q - domain%blk(b)%p_pg(Im,Jm,Km) ) + s * Ew
      domain%bc(i)%mdot_pg = - inflow / A
    enddo

  end subroutine Open_Face_Fluxes


  !> The gas-enthalpy flux of the energy equation under type 2 (plan 07 §3.5;
  !> PATO's fvc::laplacian(GammaHg, p)): + div(h_g Gamma grad p) added to
  !> blk%R with the convention of Diffusive_Flux (Res1 -= Flux, Res2 += Flux),
  !> the face tensor (g_L + g_R)/2 with g = h_g(p^(n+1), T_stage) Gamma^n; on
  !> an open face h_g(p_w, T_wall) times the pressure solve's face mass flux
  !> (T_wall from the face's BC routine, bc%Tw_face). Called by Pyrolysis_Stage in place of the type-1 line
  !> integral, after BC_Fluxes, inside the parallel region. bc%edot_pg records
  !> the enthalpy flux leaving each open face for the energy budget.
  subroutine Darcy_Energy_Flux ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Lib_Pyrolysis, only: Darcy_Cell_Coeffs, Gas_State
    use FUSS_Mod_MPI,       only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b, i, j, k, ni, nj, nk, fc, lower, upper, ii, Im, Jm, Km, Fm, m
    real(R8) :: C, E, F, Cw, Ew, s, A, q, inflow, Gw(3,3), eta_w, eps_w, rho_w, h_w, W, h, mu, Gf(3,3), Tf

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      ni = domain%blk(b)%dim(1) ; nj = domain%blk(b)%dim(2) ; nk = domain%blk(b)%dim(3)

      !$omp do collapse(3) schedule(static) private(i,j,k,m,W,h,mu)
      do k = 1, nk ; do j = 1, nj ; do i = 1, ni
        m = nint( domain%blk(b)%matID(i,j,k) )
        call Gas_State( m, domain%blk(b)%p_pg(i,j,k), domain%blk(b)%T(i,j,k), W, h, mu, [b, i, j, k] )
        wk(b)%g(:,i,j,k) = h * domain%blk(b)%gam_pg(:,i,j,k)
      enddo ; enddo ; enddo
      !$omp end do

      !$omp do collapse(2) schedule(static) private(i,j,k,C,E,F)
      do k = 1, nk ; do j = 1, nj
        do i = 1, ni - 1
          call Face_Terms( domain%blk(b), wk(b)%g, 1, i, j, k, C, E )
          F = C * ( domain%blk(b)%p_pg(i+1,j,k) - domain%blk(b)%p_pg(i,j,k) ) + E
          domain%blk(b)%R(i,j,k)   = domain%blk(b)%R(i,j,k)   - F
          domain%blk(b)%R(i+1,j,k) = domain%blk(b)%R(i+1,j,k) + F
        enddo
      enddo ; enddo
      !$omp end do

      !$omp do collapse(2) schedule(static) private(i,j,k,C,E,F)
      do k = 1, nk ; do i = 1, ni
        do j = 1, nj - 1
          call Face_Terms( domain%blk(b), wk(b)%g, 2, i, j, k, C, E )
          F = C * ( domain%blk(b)%p_pg(i,j+1,k) - domain%blk(b)%p_pg(i,j,k) ) + E
          domain%blk(b)%R(i,j,k)   = domain%blk(b)%R(i,j,k)   - F
          domain%blk(b)%R(i,j+1,k) = domain%blk(b)%R(i,j+1,k) + F
        enddo
      enddo ; enddo
      !$omp end do

      !$omp do collapse(2) schedule(static) private(i,j,k,C,E,F)
      do j = 1, nj ; do i = 1, ni
        do k = 1, nk - 1
          call Face_Terms( domain%blk(b), wk(b)%g, 3, i, j, k, C, E )
          F = C * ( domain%blk(b)%p_pg(i,j,k+1) - domain%blk(b)%p_pg(i,j,k) ) + E
          domain%blk(b)%R(i,j,k)   = domain%blk(b)%R(i,j,k)   - F
          domain%blk(b)%R(i,j,k+1) = domain%blk(b)%R(i,j,k+1) + F
        enddo
      enddo ; enddo
      !$omp end do
    enddo

    upper = 0
    do b = 1, domain%nb
      do fc = 1, 6
        lower = upper + 1
        upper = upper + domain%n_bf(b,fc)
        if ( .not. is_local_block(b) ) cycle
        !$omp do schedule(dynamic) private(ii,Im,Jm,Km,Fm,m,Cw,Ew,s,A,q,inflow,Gw,eta_w,eps_w,rho_w,h_w,W,h,mu,Gf,Tf)
        do ii = lower, upper
          Im = domain%bc(ii)%i ; Jm = domain%bc(ii)%j ; Km = domain%bc(ii)%k ; Fm = domain%bc(ii)%f
          s = real( 1 - 2 * mod(Fm,2), R8 )
          if ( Darcy_Is_Open( domain%bc(ii)%type ) ) then
            q = domain%bc(ii)%pw
            m = nint( domain%blk(b)%matID(Im,Jm,Km) )
            call Darcy_Cell_Coeffs( m, q, domain%blk(b)%TO(Im,Jm,Km), domain%blk(b)%chiO(:,Im,Jm,Km), &
                                    [b, Im, Jm, Km], eta_w, Gw, eps_w, rho_w, h_w )
            ! The gas crosses the face at the wall temperature its BC routine
            ! computed this stage (bc%Tw_face; the cell's where no wall BC set
            ! one), carried by the face mass flux of the pressure solve itself
            ! (Open_Face_Fluxes' expression): edot_pg = h_g(p_w, T_w) mdot_pg.
            ! Before 2026-10-13 the face value averaged h Gamma of the cell with
            ! h_g(p_w, T_cell) Gamma_w: the enthalpy at the cell temperature, an
            ! O(dx) flux error that made open faces first order (exact steady
            ! transpiration, test/pyrolysis/exact/transpiration).
            Tf = domain%blk(b)%T(Im,Jm,Km)
            if ( domain%bc(ii)%Tw_face > 0.0d0 ) Tf = domain%bc(ii)%Tw_face
            call Gas_State( m, q, Tf, W, h, mu, [b, Im, Jm, Km] )
            Gf = 0.5d0 * ( Tensor9( domain%blk(b)%gam_pg(:,Im,Jm,Km) ) + Gw )
            call Wall_Terms( domain%blk(b), Gf, Im, Jm, Km, Fm, .true., Cw, Ew, A )
            inflow = h * ( Cw * ( q - domain%blk(b)%p_pg(Im,Jm,Km) ) + s * Ew )
            domain%blk(b)%R(Im,Jm,Km) = domain%blk(b)%R(Im,Jm,Km) - inflow
            domain%bc(ii)%edot_pg = - inflow / A
          elseif ( domain%bc(ii)%type == 101 ) then
            ! the gas enthalpy flux of the connection: Connection_Terms on the
            ! two cells' h Gamma of this stage (plan 08 section 2.9 (3))
            q = domain%bc(ii)%pg(1)
            call Connection_Terms( domain%blk(b), domain%bc(ii), wk(b)%g(:,Im,Jm,Km), &
                                   wk(domain%bc(ii)%bs)%g(:,domain%bc(ii)%is,domain%bc(ii)%js,domain%bc(ii)%ks), &
                                   Cw, Ew, A )
            inflow = Cw * ( q - domain%blk(b)%p_pg(Im,Jm,Km) ) + s * Ew
            domain%blk(b)%R(Im,Jm,Km) = domain%blk(b)%R(Im,Jm,Km) - inflow
            domain%bc(ii)%qh_if = inflow
          endif
        enddo
        !$omp end do
      enddo
    enddo

  end subroutine Darcy_Energy_Flux


  !> Convergence-residual contribution of the pressure: sum over cells of
  !> ((p - p^n)/p-ref)^2, so that a run whose temperature is stationary while
  !> p still moves (the quasi-steady Picard iterations of gates 4.2/4.3, a
  !> solid with cp -> infinity) does not stop as converged. Outside the
  !> parallel region, like Compute_Residual.
  subroutine Darcy_Residual ( domain, total )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_darcy
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    real(R8), intent(inout) :: total
    ! Local
    integer  :: b, i, j, k
    real(R8) :: res

    res = 0.0d0
    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        res = res + ( ( domain%blk(b)%p_pg(i,j,k) - domain%blk(b)%p_pgO(i,j,k) ) / obj_darcy%p_ref )**2
      enddo ; enddo ; enddo
    enddo
    total = total + res

  end subroutine Darcy_Residual

end module FUSS_Mod_Darcy
