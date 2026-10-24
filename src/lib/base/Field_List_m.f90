!> The single list of per-cell fields every list walker consumes (plan 06 §3.1).
!>
!> Rows with managed = .true. are allocated, zeroed, freed, written and read
!> back by the walkers (Allocate_Block, deallocate_remote_computation_data,
!> Setup_Output_Solution / Write_vtk_tec / Read_vtk_tec); an unknown managed
!> name aborts at setup, so a row added without its handler fails loudly. The
!> legacy fields (T, TO, R, matID, qvol, dtlocal, RS1, RS2) are registered with
!> managed = .false. for documentation and for plan 11's exchange list: their
!> explicit allocation and I/O code is untouched, which is what keeps the
!> feature-off path byte-identical.
!>
!> The to_conserved / from_conserved pair of a row with remap = 'conserved' is
!> what plan 08's ALE remap and plan 11's exchange call. This module sits in
!> base/ and must not depend on the physics modules, so it only declares the
!> procedure-pointer components; FUSS_Mod_Pyrolysis fills them through
!> Field_List_Set_Converters right after Field_List_Init.
module FUSS_Field_List_m
  use FUSS_Mod_MPI, only: fuss_abort
  use iso_fortran_env, only: R8 => real64
  use FUSS_Advanced_Types_m, only: FUSS_block_type
  use FUSS_Global_m, only: gc

  implicit none
  private
  public :: field_entry_t, field_convert_if
  public :: Field_List_Init, Field_List_Count, Field_List_Get, Field_List_Lookup
  public :: Field_List_Set_Converters
  public :: Field_Cell_Get, Field_Cell_Set, Field_Allocate
  public :: FLD_T, FLD_CHI, FLD_P_PG, FLD_TO, FLD_RHO_S, FLD_MDOT_PG, FLD_U_PG, FLD_ETA_PREV, FLD_ETA_ST, FLD_T_PREV

  !> The arrays plan 11 moves between ranks: the halo exchange (the rows with exchange and
  !> has_ghost: T, chi, p_pg) and the gathers to root for the solution file (the output rows,
  !> and TO for the type-2 T_prev column). A row's code names its array for the one accessor
  !> pair Field_Cell_Get / Field_Cell_Set (0: no array plan 11 knows; the exchange and the
  !> gather refuse such a row).
  integer, parameter :: FLD_T = 1, FLD_CHI = 2, FLD_P_PG = 3, FLD_TO = 4, FLD_RHO_S = 5, FLD_MDOT_PG = 6, &
                        FLD_U_PG = 7, FLD_ETA_PREV = 8, FLD_ETA_ST = 9, FLD_T_PREV = 10

  abstract interface
    !> q(1:ncomp) = the conserved quantity per unit volume of cell (i,j,k)
    subroutine field_convert_if ( blk, i, j, k, q )
      import :: FUSS_block_type, R8
      type(FUSS_block_type), intent(inout) :: blk
      integer,  intent(in)    :: i, j, k
      real(R8), intent(inout) :: q(:)
    end subroutine field_convert_if
  end interface

  type :: field_entry_t
    character(len=32) :: name      = ''      ! 'chi', 'chiO', 'dchi', 'rho_s', 'k_s', 'cp_s', 'mdot_pg'; (07) 'p_pg', 'p_pgO', 'u_pg', 'p_new', 'eta_pg', 'eta_prev', 'eta_st', 'eta_stO', 'gam_pg', 'pi_pg', 'T_prev'
    character(len=4)  :: kind      = 'cell'  ! 'cell' | 'node'
    integer           :: ncomp     = 1       ! 1, or nsub for chi/chiO/dchi
    logical           :: has_ghost = .false. ! bounds 1-gc:n+gc when .true., 1:n otherwise
    character(len=9)  :: remap     = 'none'  ! 'none' | 'conserved'
    logical           :: restart   = .false. ! written AND read back on restart
    logical           :: output    = .false. ! appended to Ovarnames (ncomp columns, name_1..name_ncomp)
    logical           :: exchange  = .false. ! plan 11: halo exchange / root gather needed
    logical           :: managed   = .false. ! .true. = allocated/zeroed/freed/written by the list walkers
    integer           :: code      = 0       ! plan 11: FLD_* of the array the halo exchange moves, 0 = none
    procedure(field_convert_if), pointer, nopass :: to_conserved   => null()
    procedure(field_convert_if), pointer, nopass :: from_conserved => null()
  end type field_entry_t

  type(field_entry_t), allocatable, save :: list(:)
  integer, save :: nlist = 0

contains

  !> Build the list. Legacy rows always; the pyrolysis rows only when the
  !> feature is on (nsub = number of reactions over all materials); the Darcy
  !> rows (plan 07 section 3.1) only under gas-closure = type2.
  subroutine Field_List_Init ( nsub, with_pyrolysis, with_darcy )
    implicit none
    integer, intent(in) :: nsub
    logical, intent(in) :: with_pyrolysis, with_darcy

    if ( allocated(list) ) deallocate(list)
    allocate( list(32) )
    nlist = 0

    ! Legacy rows: explicit code paths, documented here (managed = .false.).
    call add( 'T',       1, .true.,  'conserved', .true.,  .true.,  .true.,  .false. )
    call add( 'TO',      1, .true.,  'none',      .false., .false., .false., .false. )
    call add( 'R',       1, .true.,  'none',      .false., .false., .false., .false. )
    call add( 'matID',   1, .true.,  'none',      .true.,  .true.,  .false., .false. )
    call add( 'qvol',    1, .true.,  'none',      .false., .true.,  .false., .false. )
    call add( 'dtlocal', 1, .false., 'none',      .false., .false., .false., .false. )
    call add( 'RS1',     1, .true.,  'none',      .false., .false., .false., .false. )
    call add( 'RS2',     1, .true.,  'none',      .false., .false., .false., .false. )

    if ( with_pyrolysis ) then
      ! name     ncomp  ghost    remap        restart  output   exchange managed
      call add( 'chi',     nsub, .true.,  'conserved', .true.,  .true.,  .true.,  .true. )
      call add( 'chiO',    nsub, .false., 'none',      .false., .false., .false., .true. )
      call add( 'dchi',    nsub, .false., 'none',      .false., .false., .false., .true. )
      call add( 'rho_s',   1,    .false., 'none',      .false., .true.,  .false., .true. )
      call add( 'k_s',     1,    .false., 'none',      .false., .false., .false., .true. )  ! written into the legacy "k" slot
      call add( 'cp_s',    1,    .false., 'none',      .false., .false., .false., .true. )  ! written into the legacy "cs" slot
      call add( 'mdot_pg', 1,    .false., 'none',      .false., .true.,  .false., .true. )
    endif

    if ( with_darcy ) then
      ! name     ncomp  ghost    remap        restart  output   exchange managed
      call add( 'p_pg',    1,    .true.,  'conserved', .true.,  .true.,  .true.,  .true. )   ! gas pressure [Pa]
      call add( 'p_pgO',   1,    .false., 'none',      .false., .false., .false., .true. )   ! p at step n (restart: p^n = p_pg)
      call add( 'u_pg',    3,    .false., 'none',      .false., .true.,  .false., .true. )   ! cell gas velocity [m/s]
      call add( 'p_new',   1,    .false., 'none',      .false., .false., .false., .true. )   ! sweep scratch
      call add( 'eta_pg',  1,    .false., 'none',      .false., .false., .false., .true. )   ! eta^n = eps_pg W/(R T) [kg/m3/Pa]
      call add( 'eta_prev',1,    .false., 'none',      .true.,  .true.,  .false., .true. )   ! eta^(n-1): the extrapolation's history
      call add( 'eta_st',  1,    .false., 'conserved', .true.,  .true.,  .true.,  .true. )   ! stored gas G = eta_st p_pg (a state); quasi-steady: (eta_st - eta_prev) p_pg = carried gas
      call add( 'eta_stO', 1,    .false., 'none',      .false., .false., .false., .true. )   ! eta_st at step n
      call add( 'gam_pg',  9,    .false., 'none',      .false., .false., .false., .true. )   ! Gamma^n = rho_g K/mu, row-major 3x3
      call add( 'pi_pg',   1,    .false., 'none',      .false., .false., .false., .true. )   ! pi_tot^n of the solve [kg/m3/s]
      call add( 'T_prev',  1,    .false., 'none',      .true.,  .true.,  .false., .true. )   ! T^(n-1) of the solve's mid-step extrapolation; the file column is T one step before T (TO)
    endif

  contains

    subroutine add ( name, ncomp, has_ghost, remap, restart, output, exchange, managed )
      implicit none
      character(len=*), intent(in) :: name, remap
      integer,          intent(in) :: ncomp
      logical,          intent(in) :: has_ghost, restart, output, exchange, managed
      nlist = nlist + 1
      list(nlist)%name      = name
      list(nlist)%kind      = 'cell'
      list(nlist)%ncomp     = ncomp
      list(nlist)%has_ghost = has_ghost
      list(nlist)%remap     = remap
      list(nlist)%restart   = restart
      list(nlist)%output    = output
      list(nlist)%exchange  = exchange
      list(nlist)%managed   = managed
      select case ( name )
        case ('T')        ; list(nlist)%code = FLD_T
        case ('chi')      ; list(nlist)%code = FLD_CHI
        case ('p_pg')     ; list(nlist)%code = FLD_P_PG
        case ('TO')       ; list(nlist)%code = FLD_TO
        case ('rho_s')    ; list(nlist)%code = FLD_RHO_S
        case ('mdot_pg')  ; list(nlist)%code = FLD_MDOT_PG
        case ('u_pg')     ; list(nlist)%code = FLD_U_PG
        case ('eta_prev') ; list(nlist)%code = FLD_ETA_PREV
        case ('eta_st')   ; list(nlist)%code = FLD_ETA_ST
        case ('T_prev')   ; list(nlist)%code = FLD_T_PREV
        case default      ; list(nlist)%code = 0
      end select
    end subroutine add

  end subroutine Field_List_Init


  integer function Field_List_Count () result ( n )
    implicit none
    n = nlist
  end function Field_List_Count


  function Field_List_Get ( i ) result ( e )
    implicit none
    integer, intent(in) :: i
    type(field_entry_t) :: e
    e = list(i)
  end function Field_List_Get


  !> Index of the row called `name`, 0 when absent.
  integer function Field_List_Lookup ( name ) result ( idx )
    implicit none
    character(len=*), intent(in) :: name
    integer :: i
    idx = 0
    do i = 1, nlist
      if ( trim(list(i)%name) == trim(name) ) then
        idx = i
        return
      endif
    enddo
  end function Field_List_Lookup


  !> Attach the conserved-quantity converters of a row (physics lives elsewhere).
  subroutine Field_List_Set_Converters ( name, to, from )
    implicit none
    character(len=*), intent(in) :: name
    procedure(field_convert_if)  :: to, from
    integer :: i
    i = Field_List_Lookup(name)
    if ( i == 0 ) then
      write(*,'(A)') ' [ERROR] Field_List_Set_Converters: no field called '//trim(name)
      call fuss_abort ( 1 )
    endif
    list(i)%to_conserved   => to
    list(i)%from_conserved => from
  end subroutine Field_List_Set_Converters


  !> v(1:ncomp) = the values in cell (i,j,k) of blk of the array with code `code` (plan 11:
  !> the pack of the halo exchange and of the gathers; one accessor for every field).
  subroutine Field_Cell_Get ( blk, code, i, j, k, v )
    implicit none
    type(FUSS_block_type), intent(in) :: blk
    integer,  intent(in)  :: code, i, j, k
    real(R8), intent(out) :: v(:)
    select case ( code )
      case ( FLD_T )        ; v(1) = blk%T(i,j,k)
      case ( FLD_CHI )      ; v    = blk%chi(:,i,j,k)
      case ( FLD_P_PG )     ; v(1) = blk%p_pg(i,j,k)
      case ( FLD_TO )       ; v(1) = blk%TO(i,j,k)
      case ( FLD_RHO_S )    ; v(1) = blk%rho_s(i,j,k)
      case ( FLD_MDOT_PG )  ; v(1) = blk%mdot_pg(i,j,k)
      case ( FLD_U_PG )     ; v    = blk%u_pg(:,i,j,k)
      case ( FLD_ETA_PREV ) ; v(1) = blk%eta_prev(i,j,k)
      case ( FLD_ETA_ST )   ; v(1) = blk%eta_st(i,j,k)
      case ( FLD_T_PREV )   ; v(1) = blk%T_prev(i,j,k)
    end select
  end subroutine Field_Cell_Get

  !> Cell (i,j,k) of blk of the array with code `code` = v(1:ncomp) (the unpack).
  subroutine Field_Cell_Set ( blk, code, i, j, k, v )
    implicit none
    type(FUSS_block_type), intent(inout) :: blk
    integer,  intent(in) :: code, i, j, k
    real(R8), intent(in) :: v(:)
    select case ( code )
      case ( FLD_T )        ; blk%T(i,j,k)        = v(1)
      case ( FLD_CHI )      ; blk%chi(:,i,j,k)    = v
      case ( FLD_P_PG )     ; blk%p_pg(i,j,k)     = v(1)
      case ( FLD_TO )       ; blk%TO(i,j,k)       = v(1)
      case ( FLD_RHO_S )    ; blk%rho_s(i,j,k)    = v(1)
      case ( FLD_MDOT_PG )  ; blk%mdot_pg(i,j,k)  = v(1)
      case ( FLD_U_PG )     ; blk%u_pg(:,i,j,k)   = v
      case ( FLD_ETA_PREV ) ; blk%eta_prev(i,j,k) = v(1)
      case ( FLD_ETA_ST )   ; blk%eta_st(i,j,k)   = v(1)
      case ( FLD_T_PREV )   ; blk%T_prev(i,j,k)   = v(1)
    end select
  end subroutine Field_Cell_Set

  !> The array with code `code` of blk, allocated with Allocate_Block's bounds when it is
  !> not (a remote block's copy on root, freed at setup by deallocate_remote_computation_data,
  !> before a gather fills it); zeroed when allocated here.
  subroutine Field_Allocate ( blk, code, ncomp )
    implicit none
    type(FUSS_block_type), intent(inout) :: blk
    integer, intent(in) :: code, ncomp
    integer :: ni, nj, nk
    ni = blk%dim(1) ; nj = blk%dim(2) ; nk = blk%dim(3)
    select case ( code )
      case ( FLD_T )
        if ( .not. allocated(blk%T) )        allocate( blk%T(1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc), source = 0.0d0 )
      case ( FLD_CHI )
        if ( .not. allocated(blk%chi) )      allocate( blk%chi(ncomp, 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc), source = 0.0d0 )
      case ( FLD_P_PG )
        if ( .not. allocated(blk%p_pg) )     allocate( blk%p_pg(1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc), source = 0.0d0 )
      case ( FLD_TO )
        if ( .not. allocated(blk%TO) )       allocate( blk%TO(1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc), source = 0.0d0 )
      case ( FLD_RHO_S )
        if ( .not. allocated(blk%rho_s) )    allocate( blk%rho_s(ni, nj, nk), source = 0.0d0 )
      case ( FLD_MDOT_PG )
        if ( .not. allocated(blk%mdot_pg) )  allocate( blk%mdot_pg(ni, nj, nk), source = 0.0d0 )
      case ( FLD_U_PG )
        if ( .not. allocated(blk%u_pg) )     allocate( blk%u_pg(ncomp, ni, nj, nk), source = 0.0d0 )
      case ( FLD_ETA_PREV )
        if ( .not. allocated(blk%eta_prev) ) allocate( blk%eta_prev(ni, nj, nk), source = 0.0d0 )
      case ( FLD_ETA_ST )
        if ( .not. allocated(blk%eta_st) )   allocate( blk%eta_st(ni, nj, nk), source = 0.0d0 )
      case ( FLD_T_PREV )
        if ( .not. allocated(blk%T_prev) )   allocate( blk%T_prev(ni, nj, nk), source = 0.0d0 )
    end select
  end subroutine Field_Allocate

end module FUSS_Field_List_m
