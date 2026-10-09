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
  use iso_fortran_env, only: R8 => real64
  use FUSS_Advanced_Types_m, only: FUSS_block_type

  implicit none
  private
  public :: field_entry_t, field_convert_if
  public :: Field_List_Init, Field_List_Count, Field_List_Get, Field_List_Lookup
  public :: Field_List_Set_Converters

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
      call add( 'eta_st',  1,    .false., 'conserved', .true.,  .true.,  .true.,  .true. )   ! stored gas G = eta_st p_pg (a state)
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
      error stop 1
    endif
    list(i)%to_conserved   => to
    list(i)%from_conserved => from
  end subroutine Field_List_Set_Converters

end module FUSS_Field_List_m
