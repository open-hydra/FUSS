!> Registry entries of the [FUSS-Darcy] section (plan 07 §3.8) and the parser
!> of its one list-valued key.
!>
!> The section configures the gas momentum closure selected by
!> [FUSS-Pyrolysis] gas-closure = type2: the Darcy pressure equation solved once
!> per material step (FUSS_Mod_Darcy). The keys are registered unconditionally,
!> as Register_Pyrolysis does, so an absent section leaves every one of them at
!> its default; a [FUSS-Darcy] section with another closure is refused by
!> Check_Pyrolysis_Compatibility, because its keys would be read by nothing.
!> Every key is immutable mid-run (none is in RUNTIME_MUTABLE).
!>
!> open-bc-types is a STRING on purpose: add_int_array reads ONE scalar from the
!> default text and copies it to every element (Registry.f90), so it cannot hold
!> the default '303 305 306'. Parse_Open_Codes turns the string into the code
!> list and is called from the ini-level guard (which refuses an unparsable
!> token or a code outside the allowed set) and from Assign_Setup.
module FUSS_Register_Darcy
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Config_Types_m
  use FUSS_Input_Registry

  implicit none
  private
  public :: Register_Darcy, Parse_Open_Codes

  !> BC codes a face may be opened with: the wall codes 301-306 and the
  !> coupled wall 103 (symmetry 300, connection 101 and chimera 102 never).
  integer, parameter, public :: OPENABLE_CODES(7) = [ 301, 302, 303, 304, 305, 306, 103 ]

contains

  subroutine Register_Darcy ()
    use FUSS_Parameters_m
    implicit none
    character(len=llen) :: section

    section = trim(codename)//'-Darcy'

    call reg%add( trim(section), 'solver', obj_darcy%solver, 'tridiag', &
                  'Pressure solver: tridiag = line-implicit sweeps along the open-face direction of each block '// &
                  '(pseudo on a block whose open faces lie in several directions, or none); pseudo = damped Jacobi everywhere', &
                  'tridiag ,  pseudo', .false. )

    call reg%add( trim(section), 'tol', obj_darcy%tol, '1.0e-8', &
                  'Stop test of the pressure sweeps: max cell mass residual <= tol * max(|face flux|, |source|, |storage|)', &
                  '> 0', .false. )

    call reg%add( trim(section), 'max-iter', obj_darcy%max_iter, '100000', &
                  'Sweep cap per material step; reached -> error stop 3', &
                  '> 0', .false. )

    call reg%add( trim(section), 'p-ref', obj_darcy%p_ref, '101325.0', &
                  'Pressure scale [Pa] of the reported per-sweep change (not the stop test)', &
                  '> 0', .false. )

    call reg%add( trim(section), 'p-init', obj_darcy%p_init, '101325.0', &
                  'Uniform gas pressure [Pa] on newrun; ignored on restart', &
                  '> 0', .false. )

    call reg%add( trim(section), 'p-wall', obj_darcy%p_wall, '101325.0', &
                  'Dirichlet gas pressure [Pa] on an open face whose bc.txt line carries no p_w of its own', &
                  '> 0', .false. )

    call reg%add( trim(section), 'open-bc-types', obj_darcy%open_bc_types, '303 305 306', &
                  'Blank-separated BC codes whose faces are open to the gas (p = p_w); from {301,302,303,304,305,306,103}', &
                  '', .false. )

    call reg%add( trim(section), 'quasi-steady', obj_darcy%quasi_steady, '.false.', &
                  'Drop the gas storage: d(eta p)/dt from the pressure equation and eps (rho_g h_g - p) from the energy E', &
                  'logical', .false. )

    call reg%add( trim(section), 'klinkenberg', obj_darcy%klinkenberg, '.false.', &
                  'K_eff = K + (mu/p) B with the tensors B_v, B_c of pyrolysis.dat (thesis eq. 3.15)', &
                  'logical', .false. )

    call reg%add( trim(section), 'pseudo-vnn', obj_darcy%pseudo_vnn, '0.5', &
                  'Damping of the Jacobi (pseudo) sweep: p <- p + pseudo-vnn * residual / diagonal; (0, 1]', &
                  '> 0', .false. )

    call reg%add( trim(section), 'axes-follow-grid', obj_darcy%axes_follow_grid, '.false.', &
                  'Rotate the permeability tensor with the grid lines (not built: false only)', &
                  'logical', .false. )

    call reg%add( trim(section), 'pi-source', obj_darcy%pi_source, '0.0', &
                  'TEST ONLY: constant gas production [kg/m3/s] added to the kinetics in every cell', &
                  '>= 0', .false. )

  end subroutine Register_Darcy


  !> Split open-bc-types into integer codes. ok = .false. (and bad = the
  !> offending token) when a token is not an integer or not an openable code.
  subroutine Parse_Open_Codes ( text, codes, n, ok, bad )
    implicit none
    character(len=*), intent(in)  :: text
    integer,          intent(out) :: codes(:)
    integer,          intent(out) :: n
    logical,          intent(out) :: ok
    character(len=*), intent(out) :: bad
    ! Local
    integer :: i, j, ios, c, k
    character(len=len(text)) :: tok

    codes = 0 ; n = 0 ; ok = .true. ; bad = ''
    i = 1
    do while ( i <= len_trim(text) )
      if ( text(i:i) == ' ' ) then
        i = i + 1
        cycle
      endif
      j = i
      do while ( j <= len_trim(text) )
        if ( text(j:j) == ' ' ) exit
        j = j + 1
      enddo
      tok = text(i:j-1)
      read( tok, *, iostat=ios ) c
      if ( ios /= 0 ) then
        ok = .false. ; bad = trim(tok) ; return
      endif
      if ( .not. any( OPENABLE_CODES == c ) ) then
        ok = .false. ; bad = trim(tok) ; return
      endif
      if ( n >= size(codes) ) then
        ok = .false. ; bad = trim(tok) ; return
      endif
      do k = 1, n
        if ( codes(k) == c ) exit
      enddo
      if ( k > n ) then
        n = n + 1
        codes(n) = c
      endif
      i = j
    enddo

  end subroutine Parse_Open_Codes

end module FUSS_Register_Darcy
