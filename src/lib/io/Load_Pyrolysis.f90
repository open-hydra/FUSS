!> Loader and validator of the pyrolysis material files (plan 06 §3.2):
!>   INPUT/<prefix>pyrolysis.dat   phase constants and kinetics, PATO key names
!>   INPUT/<prefix>virgin.dat      1 K rows: cp h ki kj kk emissivity absorptivity, one zone per material
!>   INPUT/<prefix>char.dat        same layout and the same T range as virgin.dat
!>   INPUT/<prefix>gas.dat         1 K rows: p M h mu, one zone per (material, p-level), p ascending
!> The Tecplot files go through ORION's tec_read_points_multivars exactly as
!> INPUT/<prefix>properties.dat does (Load_Table.f90); zone i of a solid file
!> is material i of phase.txt. Structural failures stop here on every rank
!> (error stop 1); value checks and the banner are Check_Pyrolysis_Tables (root).
module FUSS_Load_Pyrolysis
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Load_Pyrolysis_Tables, Check_Pyrolysis_Tables

  character(len=*), parameter :: ERR = ' [ERROR] pyrolysis tables: '

contains

  subroutine Load_Pyrolysis_Tables ()
    use FUSS_Config_Types_m, only: obj_table, obj_pyrolysis
    use FUSS_Global_m,       only: FUSS_phase_prefix
    use FUSS_Parameters_m
    implicit none
    ! Local
    character(len=llen) :: pfile, vfile, cfile, gfile
    integer :: m

    pfile = 'INPUT/'//trim(FUSS_phase_prefix)//'pyrolysis.dat'
    vfile = 'INPUT/'//trim(FUSS_phase_prefix)//'virgin.dat'
    cfile = 'INPUT/'//trim(FUSS_phase_prefix)//'char.dat'
    gfile = 'INPUT/'//trim(FUSS_phase_prefix)//'gas.dat'

    if ( obj_table%n < 1 ) then
      write(*,'(A)') ERR//'no material in phase.txt'
      error stop 1
    endif
    obj_pyrolysis%nmat = obj_table%n
    if ( allocated(obj_pyrolysis%mat) ) deallocate( obj_pyrolysis%mat )
    allocate( obj_pyrolysis%mat( obj_pyrolysis%nmat ) )
    if ( allocated(obj_pyrolysis%mat_nsub) ) deallocate( obj_pyrolysis%mat_nsub )
    if ( allocated(obj_pyrolysis%mat_sub_offset) ) deallocate( obj_pyrolysis%mat_sub_offset )
    allocate( obj_pyrolysis%mat_nsub( obj_pyrolysis%nmat ), obj_pyrolysis%mat_sub_offset( 0:obj_pyrolysis%nmat ) )

    call Read_Pyrolysis_Dat ( pfile )
    call Read_Solid_Table   ( vfile, .true. )
    call Read_Solid_Table   ( cfile, .false. )
    call Read_Gas_Table     ( gfile )

    ! The Newton bracket: the intersection of the ranges every Energy_Of /
    ! Gas_State lookup touches, one row in from each edge because a lookup at
    ! T reads rows int(T) and int(T)+1.
    do m = 1, obj_pyrolysis%nmat
      associate ( mt => obj_pyrolysis%mat(m) )
        mt%T_brk_lo = max( mt%Tlo_s, mt%Tlo_g ) + 1
        mt%T_brk_hi = min( mt%Thi_s, mt%Thi_g ) - 1
        if ( mt%T_brk_lo > mt%T_brk_hi ) then
          write(*,'(A)') ERR//'empty Newton bracket for material '//trim(mt%name)
          write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '         virgin/char T range [', mt%Tlo_s, ', ', mt%Thi_s, &
            '] K, gas T range [', mt%Tlo_g, ', ', mt%Thi_g, '] K'
          error stop 1
        endif
      end associate
    enddo

  end subroutine Load_Pyrolysis_Tables


  !> pyrolysis.dat: `//` comments, blank lines ignored, one `material <name>`
  !> block per material in phase.txt order, `key value` lines in any order,
  !> each key exactly once per block, no key the grammar does not name.
  subroutine Read_Pyrolysis_Dat ( file )
    use FUSS_Config_Types_m, only: obj_table, obj_pyrolysis, pyro_reaction_t
    use FUSS_Parameters_m,   only: llen, clen
    implicit none
    character(len=*), intent(in) :: file
    ! Local
    integer, parameter :: MAXENT = 8192
    character(len=64), allocatable :: keys(:), vals(:)
    integer, allocatable :: lines(:)
    logical, allocatable :: used(:)
    integer, allocatable :: bstart(:)
    character(len=llen) :: line
    integer :: nent, nb, u, ios, ln, p, m, i, j, r, b0, b1, nsub, nreac, k
    real(R8) :: sumF
    type(pyro_reaction_t) :: rc

    allocate( keys(MAXENT), vals(MAXENT), lines(MAXENT), used(MAXENT) )
    open( newunit=u, file=trim(file), status='old', action='read', iostat=ios )
    if ( ios /= 0 ) then
      write(*,'(A)') ERR//'cannot open '//trim(file)
      error stop 1
    endif

    ! Tokenise: key = first blank-delimited token, value = the rest.
    nent = 0 ; ln = 0
    do
      read( u, '(A)', iostat=ios ) line
      if ( ios /= 0 ) exit
      ln = ln + 1
      p = index( line, '//' )
      if ( p > 0 ) line = line(1:p-1)
      do k = 1, len_trim(line)
        if ( line(k:k) == achar(9) ) line(k:k) = ' '
      enddo
      line = adjustl( line )
      if ( len_trim(line) == 0 ) cycle
      if ( nent >= MAXENT ) call Fail( 'more than '//Itoa(MAXENT)//' entries', ln )
      p = index( trim(line), ' ' )
      if ( p == 0 ) call Fail( 'expected "key value", found "'//trim(line)//'"', ln )
      nent = nent + 1
      keys(nent)  = line(1:p-1)
      vals(nent)  = adjustl( line(p+1:) )
      lines(nent) = ln
    enddo
    close( u )
    used = .false.

    ! Material blocks.
    nb = count( keys(1:nent) == 'material' )
    if ( nb /= obj_table%n ) call Fail( 'expected '//Itoa(obj_table%n)//' material block(s) (phase.txt), found '//Itoa(nb), 0 )
    allocate( bstart(nb+1) )
    j = 0
    do i = 1, nent
      if ( keys(i) == 'material' ) then
        j = j + 1
        bstart(j) = i
      endif
    enddo
    bstart(nb+1) = nent + 1
    if ( nent > 0 ) then
      if ( bstart(1) /= 1 ) call Fail( 'the first entry must be "material <name>"', lines(1) )
    endif

    ! Pass 1: counts, so the flattened reaction list can be allocated.
    nsub = 0
    obj_pyrolysis%mat_sub_offset(0) = 0
    do m = 1, nb
      b0 = bstart(m) ; b1 = bstart(m+1)
      used(b0) = .true.
      associate ( mt => obj_pyrolysis%mat(m) )
        mt%name = vals(b0)(1:clen)
        if ( trim(vals(b0)) /= trim(obj_table%name(m)) ) &
          call Fail( 'material '//Itoa(m)//' is "'//trim(vals(b0))//'", phase.txt says "'//trim(obj_table%name(m))//'"', lines(b0) )
        mt%nphase = Rd_Int( 'nSolidPhases', b0, b1 )
        if ( mt%nphase < 1 ) call Fail( 'nSolidPhases must be >= 1', lines(b0) )
        allocate( mt%rho_I(mt%nphase), mt%eps_I(mt%nphase), mt%nreac_phase(mt%nphase) )
        nreac = 0
        do i = 1, mt%nphase
          mt%nreac_phase(i) = Rd_Int( 'nPyroReac['//Itoa(i)//']', b0, b1 )
          if ( mt%nreac_phase(i) < 0 ) call Fail( 'nPyroReac['//Itoa(i)//'] must be >= 0', lines(b0) )
          nreac = nreac + mt%nreac_phase(i)
        enddo
        mt%nsub = nreac
        obj_pyrolysis%mat_nsub(m) = nreac
        obj_pyrolysis%mat_sub_offset(m) = obj_pyrolysis%mat_sub_offset(m-1) + nreac
        nsub = nsub + nreac
      end associate
    enddo
    obj_pyrolysis%nsub = nsub
    if ( allocated(obj_pyrolysis%reac) ) deallocate( obj_pyrolysis%reac )
    allocate( obj_pyrolysis%reac( max(nsub,1) ) )

    ! Pass 2: everything, with the once-only and no-unknown-key rules.
    used = .false.
    r = 0
    do m = 1, nb
      b0 = bstart(m) ; b1 = bstart(m+1)
      used(b0) = .true.
      associate ( mt => obj_pyrolysis%mat(m) )
        i = Rd_Int( 'nSolidPhases', b0, b1 )
        mt%rho_v = 0.0d0
        do i = 1, mt%nphase
          mt%rho_I(i) = Rd_Real( 'rhoI['//Itoa(i)//']', b0, b1 )
          mt%eps_I(i) = Rd_Real( 'epsI['//Itoa(i)//']', b0, b1 )
          j = Rd_Int( 'nPyroReac['//Itoa(i)//']', b0, b1 )
          mt%rho_v = mt%rho_v + mt%eps_I(i) * mt%rho_I(i)
        enddo
        mt%eps_g_v = Rd_Real( 'eps_g_v', b0, b1 )
        mt%eps_g_c = Rd_Real( 'eps_g_c', b0, b1 )
        sumF = 0.0d0
        do i = 1, mt%nphase
          do j = 1, mt%nreac_phase(i)
            r = r + 1
            rc%phase  = i
            rc%F      = Rd_Real( 'F['//Itoa(i)//']['//Itoa(j)//']', b0, b1 )
            rc%A      = Rd_Real( 'A['//Itoa(i)//']['//Itoa(j)//']', b0, b1 )
            rc%E      = Rd_Real( 'E['//Itoa(i)//']['//Itoa(j)//']', b0, b1 )
            rc%m      = Rd_Real( 'm['//Itoa(i)//']['//Itoa(j)//']', b0, b1 )
            rc%b      = Rd_Real( 'n['//Itoa(i)//']['//Itoa(j)//']', b0, b1 )
            rc%T_thr  = Rd_Real( 'T['//Itoa(i)//']['//Itoa(j)//']', b0, b1 )
            rc%h_reac = Rd_Real( 'h['//Itoa(i)//']['//Itoa(j)//']', b0, b1 )
            rc%eps_rho_F = mt%eps_I(i) * mt%rho_I(i) * rc%F
            rc%w = 0.0d0
            obj_pyrolysis%reac(r) = rc
            sumF = sumF + rc%eps_rho_F
          enddo
        enddo
        mt%rho_c = mt%rho_v - sumF
        mt%hbar_guard = ( mt%nsub == 0 .or. mt%rho_c == mt%rho_v )
        if ( .not. mt%hbar_guard ) then
          do k = obj_pyrolysis%mat_sub_offset(m-1) + 1, obj_pyrolysis%mat_sub_offset(m)
            obj_pyrolysis%reac(k)%w = obj_pyrolysis%reac(k)%eps_rho_F / ( mt%rho_v - mt%rho_c )
          enddo
        endif
      end associate
      do i = b0 + 1, b1 - 1
        if ( .not. used(i) ) call Fail( 'unknown key "'//trim(keys(i))//'" in material block '//Itoa(m), lines(i) )
      enddo
    enddo

  contains

    subroutine Fail ( reason, ln )
      implicit none
      character(len=*), intent(in) :: reason
      integer,          intent(in) :: ln
      if ( ln > 0 ) then
        write(*,'(A)') ERR//trim(file)//' line '//Itoa(ln)//': '//reason
      else
        write(*,'(A)') ERR//trim(file)//': '//reason
      endif
      error stop 1
    end subroutine Fail

    function Find_Key ( key, b0, b1 ) result ( idx )
      implicit none
      character(len=*), intent(in) :: key
      integer,          intent(in) :: b0, b1
      integer :: idx, e
      idx = 0
      do e = b0 + 1, b1 - 1
        if ( trim(keys(e)) == trim(key) ) then
          if ( idx /= 0 ) call Fail( 'key "'//trim(key)//'" given twice (first at line '//Itoa(lines(idx))//')', lines(e) )
          idx = e
        endif
      enddo
      if ( idx == 0 ) call Fail( 'missing key "'//trim(key)//'" in material block starting at line '//Itoa(lines(b0)), 0 )
      used(idx) = .true.
    end function Find_Key

    function Rd_Real ( key, b0, b1 ) result ( v )
      implicit none
      character(len=*), intent(in) :: key
      integer,          intent(in) :: b0, b1
      real(R8) :: v
      integer  :: idx, ios
      idx = Find_Key( key, b0, b1 )
      read( vals(idx), *, iostat=ios ) v
      if ( ios /= 0 ) call Fail( 'non-numeric value "'//trim(vals(idx))//'" for key "'//trim(key)//'"', lines(idx) )
      if ( v /= v ) call Fail( 'value of "'//trim(key)//'" is not finite', lines(idx) )
    end function Rd_Real

    function Rd_Int ( key, b0, b1 ) result ( v )
      implicit none
      character(len=*), intent(in) :: key
      integer,          intent(in) :: b0, b1
      integer :: v, idx, ios
      idx = Find_Key( key, b0, b1 )
      read( vals(idx), *, iostat=ios ) v
      if ( ios /= 0 ) call Fail( 'non-integer value "'//trim(vals(idx))//'" for key "'//trim(key)//'"', lines(idx) )
    end function Rd_Int

  end subroutine Read_Pyrolysis_Dat


  !> virgin.dat / char.dat: one zone per material, integer T rows 1 K apart,
  !> columns cp h ki kj kk emissivity absorptivity.
  subroutine Read_Solid_Table ( file, virgin )
    use FUSS_Config_Types_m, only: obj_pyrolysis
    use Lib_Tecplot,         only: tec_read_points_multivars
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in) :: file
    logical,          intent(in) :: virgin
    ! Local
    type(ORION_data) :: orion
    integer  :: ios, m, i, n, Tlo, Thi, v
    real(R8) :: T1, val

    ios = tec_read_points_multivars( orion, 7, trim(file) )
    if ( ios /= 0 .or. .not. allocated(orion%block) ) then
      write(*,'(A)') ERR//'cannot read '//trim(file)
      error stop 1
    endif
    if ( size(orion%block) /= obj_pyrolysis%nmat ) then
      write(*,'(A)') ERR//trim(file)//': '//Itoa(size(orion%block))//' zone(s), expected one per material ('// &
        Itoa(obj_pyrolysis%nmat)//')'
      error stop 1
    endif

    do m = 1, obj_pyrolysis%nmat
      n  = orion%block(m)%Ni
      if ( n < 2 ) call Fail( m, 1, 'fewer than two rows' )
      T1 = orion%block(m)%mesh(1,1,1,1)
      if ( abs( T1 - nint(T1) ) > 1.0d-9 ) call Fail( m, 1, 'T must be an integer' )
      Tlo = nint(T1) ; Thi = Tlo + n - 1
      do i = 2, n
        if ( abs( orion%block(m)%mesh(1,i,1,1) - T1 - (i-1) ) > 1.0d-9 ) call Fail( m, i, 'T must step by exactly 1 K' )
      enddo
      do i = 1, n
        do v = 1, 7
          val = orion%block(m)%vars(v,i,1,1)
          if ( val /= val .or. abs(val) > huge(val) ) call Fail( m, i, 'value is not finite' )
        enddo
        if ( orion%block(m)%vars(3,i,1,1) /= orion%block(m)%vars(4,i,1,1) .or. &
             orion%block(m)%vars(3,i,1,1) /= orion%block(m)%vars(5,i,1,1) ) &
          call Fail( m, i, 'ki, kj, kk differ: orthotropic conductivity is not supported' )
      enddo
      associate ( mt => obj_pyrolysis%mat(m) )
        if ( virgin ) then
          mt%Tlo_s = Tlo ; mt%Thi_s = Thi
          allocate( mt%cp_v(Tlo:Thi), mt%h_v(Tlo:Thi), mt%k_v(Tlo:Thi), mt%em_v(Tlo:Thi), mt%ab_v(Tlo:Thi) )
          mt%cp_v(Tlo:Thi) = orion%block(m)%vars(1,1:n,1,1)
          mt%h_v (Tlo:Thi) = orion%block(m)%vars(2,1:n,1,1)
          mt%k_v (Tlo:Thi) = orion%block(m)%vars(3,1:n,1,1)
          mt%em_v(Tlo:Thi) = orion%block(m)%vars(6,1:n,1,1)
          mt%ab_v(Tlo:Thi) = orion%block(m)%vars(7,1:n,1,1)
        else
          if ( Tlo /= mt%Tlo_s .or. Thi /= mt%Thi_s ) then
            write(*,'(A)') ERR//trim(file)//' zone '//Itoa(m)//': T range ['//Itoa(Tlo)//', '//Itoa(Thi)// &
              '] K differs from virgin.dat ['//Itoa(mt%Tlo_s)//', '//Itoa(mt%Thi_s)//'] K'
            error stop 1
          endif
          allocate( mt%cp_c(Tlo:Thi), mt%h_c(Tlo:Thi), mt%k_c(Tlo:Thi), mt%em_c(Tlo:Thi), mt%ab_c(Tlo:Thi) )
          mt%cp_c(Tlo:Thi) = orion%block(m)%vars(1,1:n,1,1)
          mt%h_c (Tlo:Thi) = orion%block(m)%vars(2,1:n,1,1)
          mt%k_c (Tlo:Thi) = orion%block(m)%vars(3,1:n,1,1)
          mt%em_c(Tlo:Thi) = orion%block(m)%vars(6,1:n,1,1)
          mt%ab_c(Tlo:Thi) = orion%block(m)%vars(7,1:n,1,1)
        endif
      end associate
    enddo

  contains

    subroutine Fail ( zone, row, reason )
      implicit none
      integer,          intent(in) :: zone, row
      character(len=*), intent(in) :: reason
      write(*,'(A)') ERR//trim(file)//' zone '//Itoa(zone)//' row '//Itoa(row)//': '//reason
      error stop 1
    end subroutine Fail

  end subroutine Read_Solid_Table


  !> gas.dat: Np zones per material (material-major), p constant within a zone
  !> and strictly increasing across a material's zones, identical integer T
  !> rows 1 K apart in every zone of a material; columns p M h mu.
  subroutine Read_Gas_Table ( file )
    use FUSS_Config_Types_m, only: obj_pyrolysis
    use Lib_Tecplot,         only: tec_read_points_multivars
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in) :: file
    ! Local
    type(ORION_data) :: orion
    integer  :: ios, m, i, n, Tlo, Thi, v, nz, Np, ip, z
    real(R8) :: T1, p, val

    ios = tec_read_points_multivars( orion, 4, trim(file) )
    if ( ios /= 0 .or. .not. allocated(orion%block) ) then
      write(*,'(A)') ERR//'cannot read '//trim(file)
      error stop 1
    endif
    nz = size(orion%block)
    if ( nz < 1 .or. mod( nz, obj_pyrolysis%nmat ) /= 0 ) then
      write(*,'(A)') ERR//trim(file)//': '//Itoa(nz)//' zone(s) is not a multiple of the number of materials ('// &
        Itoa(obj_pyrolysis%nmat)//')'
      error stop 1
    endif
    Np = nz / obj_pyrolysis%nmat

    do m = 1, obj_pyrolysis%nmat
      associate ( mt => obj_pyrolysis%mat(m) )
        mt%Np = Np
        allocate( mt%p_lev(Np) )
        do ip = 1, Np
          z = (m-1)*Np + ip
          n = orion%block(z)%Ni
          if ( n < 2 ) call Fail( z, 1, 'fewer than two rows' )
          T1 = orion%block(z)%mesh(1,1,1,1)
          if ( abs( T1 - nint(T1) ) > 1.0d-9 ) call Fail( z, 1, 'T must be an integer' )
          Tlo = nint(T1) ; Thi = Tlo + n - 1
          do i = 2, n
            if ( abs( orion%block(z)%mesh(1,i,1,1) - T1 - (i-1) ) > 1.0d-9 ) call Fail( z, i, 'T must step by exactly 1 K' )
          enddo
          p = orion%block(z)%vars(1,1,1,1)
          do i = 1, n
            do v = 1, 4
              val = orion%block(z)%vars(v,i,1,1)
              if ( val /= val .or. abs(val) > huge(val) ) call Fail( z, i, 'value is not finite' )
            enddo
            if ( abs( orion%block(z)%vars(1,i,1,1) - p ) > 1.0d-12 * abs(p) ) call Fail( z, i, 'p is not constant within the zone' )
          enddo
          if ( ip == 1 ) then
            mt%Tlo_g = Tlo ; mt%Thi_g = Thi
            allocate( mt%W_g(Tlo:Thi,Np), mt%h_g(Tlo:Thi,Np), mt%mu_g(Tlo:Thi,Np) )
          else
            if ( Tlo /= mt%Tlo_g .or. Thi /= mt%Thi_g ) call Fail( z, 1, 'T range differs from the material''s first gas zone' )
            if ( p <= mt%p_lev(ip-1) ) call Fail( z, 1, 'p-levels must be strictly increasing across a material''s zones' )
          endif
          mt%p_lev(ip) = p
          mt%W_g (Tlo:Thi,ip) = orion%block(z)%vars(2,1:n,1,1)
          mt%h_g (Tlo:Thi,ip) = orion%block(z)%vars(3,1:n,1,1)
          mt%mu_g(Tlo:Thi,ip) = orion%block(z)%vars(4,1:n,1,1)
        enddo
      end associate
    enddo

  contains

    subroutine Fail ( zone, row, reason )
      implicit none
      integer,          intent(in) :: zone, row
      character(len=*), intent(in) :: reason
      write(*,'(A)') ERR//trim(file)//' zone '//Itoa(zone)//' row '//Itoa(row)//': '//reason
      error stop 1
    end subroutine Fail

  end subroutine Read_Gas_Table


  !> Value checks of pyrolysis.dat (root, after the loader) and the banner.
  subroutine Check_Pyrolysis_Tables ()
    use FUSS_Config_Types_m, only: obj_pyrolysis
    use FUSS_Lib_Pyrolysis,  only: Gas_State, Tab_Lookup
    implicit none
    ! Local
    integer  :: m, i, j, r
    real(R8) :: sumF, sumeps, hv, hc, hg, W, mu
    logical  :: bad, has_hg
    character(len=256) :: msg

    bad = .false.
    do m = 1, obj_pyrolysis%nmat
      associate ( mt => obj_pyrolysis%mat(m) )
        if ( mt%nphase < 1 ) call Flag( mt%name, 'nSolidPhases', 'must be >= 1' )
        sumeps = 0.0d0
        do i = 1, mt%nphase
          if ( mt%rho_I(i) <= 0.0d0 ) call Flag( mt%name, 'rhoI['//Itoa(i)//']', 'must be > 0' )
          if ( mt%eps_I(i) <= 0.0d0 .or. mt%eps_I(i) > 1.0d0 ) call Flag( mt%name, 'epsI['//Itoa(i)//']', 'must be in (0, 1]' )
          sumeps = sumeps + mt%eps_I(i)
        enddo
        if ( mt%eps_g_v <= 0.0d0 .or. mt%eps_g_v > 1.0d0 ) call Flag( mt%name, 'eps_g_v', 'must be in (0, 1]' )
        if ( mt%eps_g_c <= 0.0d0 .or. mt%eps_g_c > 1.0d0 ) call Flag( mt%name, 'eps_g_c', 'must be in (0, 1]' )
        if ( sumeps + mt%eps_g_v > 1.0d0 + 1.0d-12 ) call Flag( mt%name, 'epsI + eps_g_v', 'sum exceeds 1' )
        r = obj_pyrolysis%mat_sub_offset(m-1)
        do i = 1, mt%nphase
          sumF = 0.0d0
          do j = 1, mt%nreac_phase(i)
            r = r + 1
            associate ( rc => obj_pyrolysis%reac(r) )
              if ( rc%F < 0.0d0 .or. rc%F > 1.0d0 ) call Flag( mt%name, 'F['//Itoa(i)//']['//Itoa(j)//']', 'must be in [0, 1]' )
              if ( rc%A < 0.0d0 )     call Flag( mt%name, 'A['//Itoa(i)//']['//Itoa(j)//']', 'must be >= 0' )
              if ( rc%E < 0.0d0 )     call Flag( mt%name, 'E['//Itoa(i)//']['//Itoa(j)//']', 'must be >= 0' )
              if ( rc%m < 0.0d0 )     call Flag( mt%name, 'm['//Itoa(i)//']['//Itoa(j)//']', 'must be >= 0' )
              if ( rc%T_thr < 0.0d0 ) call Flag( mt%name, 'T['//Itoa(i)//']['//Itoa(j)//']', 'must be >= 0' )
              sumF = sumF + rc%F
            end associate
          enddo
          if ( sumF > 1.0d0 + 1.0d-12 ) call Flag( mt%name, 'F['//Itoa(i)//'][*]', 'sum over the reactions of a phase exceeds 1' )
        enddo
      end associate
    enddo

    if ( bad ) then
      write(*,'(A,T35,A)') '   Pyrolysis tables', 'FAIL'
      error stop 1
    endif

    write(*,'(A,T35,A)') '   Pyrolysis tables', 'OK'
    do m = 1, obj_pyrolysis%nmat
      associate ( mt => obj_pyrolysis%mat(m) )
        write(*,'(A,A,A,I0,A,I0,A,F10.3,A,F10.3,A)') '     material ', trim(mt%name), ': phases ', mt%nphase, &
          ', reactions ', mt%nsub, ', rho_v = ', mt%rho_v, ', rho_c = ', mt%rho_c, ' kg/m3'
        hv = 0.0d0 ; hc = 0.0d0
        if ( mt%Tlo_s <= 298 .and. mt%Thi_s >= 299 ) then
          hv = Tab_Lookup( mt%h_v, mt%Tlo_s, 298.0d0 )
          hc = Tab_Lookup( mt%h_c, mt%Tlo_s, 298.0d0 )
        endif
        has_hg = ( mt%Tlo_g <= 298 .and. mt%Thi_g > 298 .and. obj_pyrolysis%p_gas >= mt%p_lev(1) .and. &
                   obj_pyrolysis%p_gas <= mt%p_lev(mt%Np) )
        if ( has_hg ) then
          call Gas_State( m, obj_pyrolysis%p_gas, 298.0d0, W, hg, mu )
          write(msg,'(A,ES12.4,A,ES12.4,A,ES12.4,A)') '       h_v(298) = ', hv, '  h_c(298) = ', hc, &
            '  h_g(298, p_gas) = ', hg, ' J/kg'
        else
          write(msg,'(A,ES12.4,A,ES12.4,A)') '       h_v(298) = ', hv, '  h_c(298) = ', hc, '  h_g(298, p_gas) = n/a'
        endif
        write(*,'(A)') trim(msg)
        write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') '       tables: solid ', mt%Tlo_s, '-', mt%Thi_s, ' K, gas ', &
          mt%Tlo_g, '-', mt%Thi_g, ' K x ', mt%Np, ' p-level(s), Newton bracket [', mt%T_brk_lo, ', ', mt%T_brk_hi, '] K'
        if ( mt%nsub > 0 ) write(*,'(A)') '       h[i][j] (heat of pyrolysis per reaction) parsed and ignored: the decomposition heat is in h_v, h_c'
        if ( obj_pyrolysis%hbar_guard_report .and. mt%hbar_guard ) &
          write(*,'(A)') '       hbar guard fired: rho_v = rho_c, tau = 1 (no net decomposition)'
      end associate
    enddo

  contains

    subroutine Flag ( name, key, reason )
      implicit none
      character(len=*), intent(in) :: name, key, reason
      write(*,'(A)') ERR//'material '//trim(name)//', key '//trim(key)//': '//reason
      bad = .true.
    end subroutine Flag

  end subroutine Check_Pyrolysis_Tables


  function Itoa ( i ) result ( s )
    implicit none
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=16) :: buf
    write( buf, '(I0)' ) i
    s = trim(buf)
  end function Itoa

end module FUSS_Load_Pyrolysis
