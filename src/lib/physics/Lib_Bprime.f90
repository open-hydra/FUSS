!> B' tables of the ablating wall (plan 08 section 2.1): PATO's five-column
!> file, read and interpolated as PATO 3.1 does (BprimeTable::BprimeC).
!>
!> FILE. `//` lines are comments; blank lines are skipped; every other line
!> holds exactly 5 reals `p B'g B'c T h_w` (6 with `// axes: rhoeueCM p
!> Bprime-g T`: rhoe ue C_M first). `// kind: equilibrium | finite-rate |
!> forced-equilibrium` is provenance (absent: unknown, with a warning).
!>
!> ORDER (PATO's searches assume it and mis-interpolate silently otherwise, so
!> it is checked here): p-blocks = runs of equal p, strictly decreasing; in a
!> block, B'g sub-blocks = runs of equal B'g, strictly decreasing; in a
!> sub-block at least two rows with T increasing or strictly decreasing
!> (stored ascending). An increasing sub-block may repeat a T: createBprimeTable
!> prints T with 6 significant digits, and with the 2000 B'c points a 1 K
!> table needs, 3633 of 59586 consecutive rows tie (2026-10-13). PATO takes
!> the first row at or above the query, and so does this lookup. Each
!> sub-block has its own T grid and row count; the table need not be
!> rectangular. A 4-D file stacks such tables, one per rhoe ue C_M level,
!> levels ascending.
!>
!> LOOKUP. p clamped to the table and interpolated linearly in ln p between
!> the two bracketing blocks; in each block B'g clamped to that block and
!> interpolated linearly between its two bracketing sub-blocks; in each
!> sub-block T clamped and B'c, h_w interpolated linearly in T; every step as
!> PATO writes it, f*v_upper + (1 - f)*v_lower, and a query on a grid value
!> returns the table value. A 4-D table evaluates the two bracketing levels
!> and interpolates linearly in ln(rhoe ue C_M). Clamping is reported per axis
!> (p, B'g, T, rhoe ue C_M) to the caller, which counts it.
!>
!> The registry is filled at setup and only read afterwards, so the lookup is
!> thread-safe.
module FUSS_Lib_Bprime
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Load_Bprime_Table, Bprime_Lookup, Bprime_Describe, Bprime_Free, Bprime_Ndim, Bprime_Count

  character(len=*), parameter :: ERR = ' [ERROR] B'' table '

  type :: bprime_3d_t
    integer :: np = 0
    real(R8), allocatable :: p(:)          ! (np) strictly decreasing
    integer,  allocatable :: g_off(:)      ! (np+1) sub-blocks of block ip: g_off(ip)+1 .. g_off(ip+1)
    real(R8), allocatable :: bg(:)         ! one per sub-block, strictly decreasing inside a block
    integer,  allocatable :: t_off(:)      ! (nsub+1) rows of sub-block s: t_off(s)+1 .. t_off(s+1)
    real(R8), allocatable :: T(:), bc(:), hw(:)   ! rows, T non-decreasing inside a sub-block
  end type bprime_3d_t

  type :: bprime_table_t
    character(len=256) :: file = ''
    character(len=20)  :: kind = 'unknown'
    integer :: ndim = 3
    integer :: nr   = 1
    real(R8), allocatable :: rc(:)         ! (nr) rhoe ue C_M levels, strictly increasing (4-D)
    type(bprime_3d_t), allocatable :: tab(:)
    real(R8) :: lim(2,4) = 0.0d0           ! (min,max) of p, B'g, T, rhoe ue C_M over the file
  end type bprime_table_t

  type(bprime_table_t), allocatable, save :: reg(:)
  integer, save :: nreg = 0

contains

  !> Index of the table read from `path`, loading it on first use.
  integer function Load_Bprime_Table ( path ) result ( it )
    implicit none
    character(len=*), intent(in) :: path
    ! Local
    type(bprime_table_t), allocatable :: tmp(:)
    integer :: i

    do i = 1, nreg
      if ( trim(reg(i)%file) == trim(path) ) then
        it = i
        return
      endif
    enddo
    if ( .not. allocated(reg) ) allocate( reg(4) )
    if ( nreg == size(reg) ) then
      allocate( tmp(2*nreg) )
      tmp(1:nreg) = reg(1:nreg)
      call move_alloc( tmp, reg )
    endif
    nreg = nreg + 1
    it = nreg
    call Read_Table( path, reg(it) )

  end function Load_Bprime_Table


  !> Number of tables loaded so far.
  integer function Bprime_Count () result ( n )
    implicit none
    n = nreg
  end function Bprime_Count


  integer function Bprime_Ndim ( it ) result ( n )
    implicit none
    integer, intent(in) :: it
    n = reg(it)%ndim
  end function Bprime_Ndim


  !> One banner line per table: file, kind, layout and ranges.
  function Bprime_Describe ( it ) result ( s )
    implicit none
    integer, intent(in) :: it
    character(len=512) :: s
    associate ( t => reg(it) )
      if ( t%ndim == 4 ) then
        write(s,'(A,4(A,ES10.3,A,ES10.3),A)') trim(t%file)//': kind '//trim(t%kind)//', 4-D', &
          ', p ', t%lim(1,1), '..', t%lim(2,1), ' Pa, B''g ', t%lim(1,2), '..', t%lim(2,2), &
          ', T ', t%lim(1,3), '..', t%lim(2,3), ' K, rhoeueCM ', t%lim(1,4), '..', t%lim(2,4), ' kg/m2/s'
      else
        write(s,'(A,3(A,ES10.3,A,ES10.3),A)') trim(t%file)//': kind '//trim(t%kind)//', 3-D', &
          ', p ', t%lim(1,1), '..', t%lim(2,1), ' Pa, B''g ', t%lim(1,2), '..', t%lim(2,2), &
          ', T ', t%lim(1,3), '..', t%lim(2,3), ' K'
      endif
    end associate
  end function Bprime_Describe


  subroutine Bprime_Free ()
    implicit none
    if ( allocated(reg) ) deallocate( reg )
    nreg = 0
  end subroutine Bprime_Free


  !> B'c and h_w at (p, B'g, T) [and rhoe ue C_M for a 4-D table]; clamped(1:4)
  !> = the query left the table range in p, B'g, T, rhoe ue C_M.
  pure subroutine Bprime_Lookup ( it, p, bg, T, rc, Bc, hw, clamped )
    implicit none
    integer,  intent(in)  :: it
    real(R8), intent(in)  :: p, bg, T, rc
    real(R8), intent(out) :: Bc, hw
    logical,  intent(out) :: clamped(4)
    ! Local
    integer  :: lo, hi
    real(R8) :: x, f, b0, h0, b1, h1

    clamped = .false.
    associate ( tb => reg(it) )
      if ( tb%ndim == 3 ) then
        call Lookup_3d( tb%tab(1), p, bg, T, Bc, hw, clamped )
        return
      endif
      x = rc
      if ( x < tb%rc(1) .or. x > tb%rc(tb%nr) ) clamped(4) = .true.
      x = max( tb%rc(1), min( tb%rc(tb%nr), x ) )
      hi = 1
      do while ( tb%rc(hi) < x )
        hi = hi + 1
      enddo
      if ( tb%rc(hi) == x ) then
        call Lookup_3d( tb%tab(hi), p, bg, T, Bc, hw, clamped )
        return
      endif
      lo = hi - 1
      f  = ( log(x) - log(tb%rc(lo)) ) / ( log(tb%rc(hi)) - log(tb%rc(lo)) )
      call Lookup_3d( tb%tab(lo), p, bg, T, b0, h0, clamped )
      call Lookup_3d( tb%tab(hi), p, bg, T, b1, h1, clamped )
      Bc = f * b1 + (1.0d0 - f) * b0
      hw = f * h1 + (1.0d0 - f) * h0
    end associate

  end subroutine Bprime_Lookup


  !> PATO's BprimeTable::BprimeC on one 3-D table.
  pure subroutine Lookup_3d ( tab3, p, bg, Tq, Bc, hw, clamped )
    implicit none
    type(bprime_3d_t), intent(in)    :: tab3
    real(R8),          intent(in)    :: p, bg, Tq
    real(R8),          intent(out)   :: Bc, hw
    logical,           intent(inout) :: clamped(4)
    ! Local
    integer  :: lo, hi
    real(R8) :: x, f, b0, h0, b1, h1

    x = p
    if ( x > tab3%p(1) .or. x < tab3%p(tab3%np) ) clamped(1) = .true.
    x = min( tab3%p(1), max( tab3%p(tab3%np), x ) )
    lo = 1
    do while ( tab3%p(lo) > x )
      lo = lo + 1
    enddo
    if ( tab3%p(lo) == x ) then
      call Lookup_Block( tab3, lo, bg, Tq, Bc, hw, clamped )
      return
    endif
    hi = lo - 1
    f  = ( log(x) - log(tab3%p(lo)) ) / ( log(tab3%p(hi)) - log(tab3%p(lo)) )
    call Lookup_Block( tab3, lo, bg, Tq, b0, h0, clamped )
    call Lookup_Block( tab3, hi, bg, Tq, b1, h1, clamped )
    hw = f * h1 + (1.0d0 - f) * h0
    Bc = f * b1 + (1.0d0 - f) * b0

  end subroutine Lookup_3d


  !> B'g step inside p-block ip (PATO indexBprimeG): sub-blocks of decreasing
  !> B'g, the query clamped to the block's own range.
  pure subroutine Lookup_Block ( tab3, ip, bg, Tq, Bc, hw, clamped )
    implicit none
    type(bprime_3d_t), intent(in)    :: tab3
    integer,           intent(in)    :: ip
    real(R8),          intent(in)    :: bg, Tq
    real(R8),          intent(out)   :: Bc, hw
    logical,           intent(inout) :: clamped(4)
    ! Local
    integer  :: s0, s1, slo, shi
    real(R8) :: g, fg, c0, w0, c1, w1

    s0 = tab3%g_off(ip) + 1          ! largest B'g of the block
    s1 = tab3%g_off(ip+1)            ! smallest
    g  = bg
    if ( g < tab3%bg(s1) .or. g > tab3%bg(s0) ) clamped(2) = .true.
    g  = max( tab3%bg(s1), g )
    g  = min( tab3%bg(s0), g )
    slo = s0
    do while ( tab3%bg(slo) > g )
      slo = slo + 1
    enddo
    if ( tab3%bg(slo) == g ) then
      call Lookup_Sub( tab3, slo, Tq, Bc, hw, clamped )
      return
    endif
    shi = slo - 1
    fg  = ( g - tab3%bg(slo) ) / ( tab3%bg(shi) - tab3%bg(slo) )
    call Lookup_Sub( tab3, slo, Tq, c0, w0, clamped )
    call Lookup_Sub( tab3, shi, Tq, c1, w1, clamped )
    Bc = fg * c1 + (1.0d0 - fg) * c0
    hw = fg * w1 + (1.0d0 - fg) * w0

  end subroutine Lookup_Block


  !> T step inside sub-block s (PATO indexTi, increasing branch; rows are
  !> stored increasing): the first row at or above the clamped T.
  pure subroutine Lookup_Sub ( tab3, s, Tq, Bc, hw, clamped )
    implicit none
    type(bprime_3d_t), intent(in)    :: tab3
    integer,           intent(in)    :: s
    real(R8),          intent(in)    :: Tq
    real(R8),          intent(out)   :: Bc, hw
    logical,           intent(inout) :: clamped(4)
    ! Local
    integer  :: r0, r1, k
    real(R8) :: y, ft

    r0 = tab3%t_off(s) + 1
    r1 = tab3%t_off(s+1)
    y  = Tq
    if ( y > tab3%T(r1) .or. y < tab3%T(r0) ) clamped(3) = .true.
    y  = min( tab3%T(r1), max( tab3%T(r0), y ) )
    k  = r0
    do while ( tab3%T(k) < y )
      k = k + 1
    enddo
    if ( tab3%T(k) == y ) then
      Bc = tab3%bc(k)
      hw = tab3%hw(k)
      return
    endif
    ft = ( y - tab3%T(k-1) ) / ( tab3%T(k) - tab3%T(k-1) )
    Bc = ft * tab3%bc(k) + (1.0d0 - ft) * tab3%bc(k-1)
    hw = ft * tab3%hw(k) + (1.0d0 - ft) * tab3%hw(k-1)

  end subroutine Lookup_Sub


  ! ---------------------------------------------------------------- reader

  subroutine Read_Table ( path, tb )
    implicit none
    character(len=*),     intent(in)    :: path
    type(bprime_table_t), intent(inout) :: tb
    ! Local
    integer, parameter :: LL = 1024
    character(len=LL) :: line, t
    integer  :: u, ios, ln, nrow, ncol, n, i, a, b, lev, nlev
    real(R8), allocatable :: d(:,:)
    integer,  allocatable :: lnum(:)
    real(R8) :: v(7)
    logical  :: have_kind, have_axes

    tb%file = path
    open( newunit=u, file=trim(path), status='old', action='read', iostat=ios )
    if ( ios /= 0 ) then
      write(*,'(A)') ERR//trim(path)//': cannot open'
      error stop 1
    endif

    ! Pass 1: header lines and the row count
    ncol = 5 ; nrow = 0 ; ln = 0 ; have_kind = .false. ; have_axes = .false.
    do
      read( u, '(A)', iostat=ios ) line
      if ( ios /= 0 ) exit
      ln = ln + 1
      t = adjustl(line)
      if ( len_trim(t) == 0 ) cycle
      if ( t(1:2) == '//' ) then
        t = adjustl( t(3:) )
        if ( t(1:5) == 'kind:' ) then
          tb%kind = adjustl( t(6:) )
          select case ( trim(tb%kind) )
            case ( 'equilibrium', 'finite-rate', 'forced-equilibrium' )
            case default
              call Fail( ln, 'unknown kind "'//trim(tb%kind)//'" (equilibrium, finite-rate, forced-equilibrium)' )
          end select
          have_kind = .true.
        elseif ( t(1:5) == 'axes:' ) then
          select case ( Squeeze( t(6:) ) )
            case ( 'p Bprime-g T' )
              ncol = 5
            case ( 'rhoeueCM p Bprime-g T' )
              ncol = 6
            case default
              call Fail( ln, 'axes must be "p Bprime-g T" or "rhoeueCM p Bprime-g T"' )
          end select
          if ( nrow > 0 ) call Fail( ln, 'the axes line must precede the data' )
          have_axes = .true.
        endif
        cycle
      endif
      nrow = nrow + 1
    enddo
    if ( nrow < 2 ) call Fail( 0, 'fewer than two data rows' )
    if ( .not. have_kind ) write(*,'(A)') ' [WARNING] B'' table '//trim(path)//': no "// kind:" line, kind = unknown'
    tb%ndim = merge( 4, 3, ncol == 6 )

    ! Pass 2: the data rows, each exactly ncol reals
    allocate( d(ncol, nrow), lnum(nrow) )
    rewind( u )
    n = 0 ; ln = 0
    do
      read( u, '(A)', iostat=ios ) line
      if ( ios /= 0 ) exit
      ln = ln + 1
      t = adjustl(line)
      if ( len_trim(t) == 0 ) cycle
      if ( t(1:2) == '//' ) cycle
      if ( Count_Tokens(t) /= ncol ) &
        call Fail( ln, 'expected '//Itoa(ncol)//' numbers, found "'//trim(t)//'"' )
      read( t, *, iostat=ios ) v(1:ncol)
      if ( ios /= 0 ) call Fail( ln, 'not a row of numbers: "'//trim(t)//'"' )
      do i = 1, ncol
        if ( v(i) /= v(i) .or. abs(v(i)) > huge(1.0d0) ) call Fail( ln, 'value is not finite' )
      enddo
      n = n + 1
      d(:,n) = v(1:ncol)
      lnum(n) = ln
    enddo
    close( u )

    ! 4-D: split into levels of equal rhoe ue C_M, strictly increasing
    if ( ncol == 6 ) then
      nlev = 1
      do i = 2, nrow
        if ( d(1,i) /= d(1,i-1) ) then
          if ( d(1,i) < d(1,i-1) ) call Fail( lnum(i), 'rhoeueCM must increase from level to level' )
          nlev = nlev + 1
        endif
      enddo
      tb%nr = nlev
      allocate( tb%rc(nlev), tb%tab(nlev) )
      a = 1 ; lev = 0
      do while ( a <= nrow )
        b = a
        do while ( b < nrow )
          if ( d(1,b+1) /= d(1,a) ) exit
          b = b + 1
        enddo
        lev = lev + 1
        tb%rc(lev) = d(1,a)
        if ( tb%rc(lev) <= 0.0d0 ) call Fail( lnum(a), 'rhoeueCM must be > 0' )
        call Build_3d( d(2:6, a:b), lnum(a:b), tb%tab(lev) )
        a = b + 1
      enddo
      tb%lim(:,4) = [ tb%rc(1), tb%rc(nlev) ]
    else
      tb%nr = 1
      allocate( tb%tab(1) )
      call Build_3d( d, lnum, tb%tab(1) )
    endif
    tb%lim(1,1) = minval( d(ncol-4,:) ) ; tb%lim(2,1) = maxval( d(ncol-4,:) )
    tb%lim(1,2) = minval( d(ncol-3,:) ) ; tb%lim(2,2) = maxval( d(ncol-3,:) )
    tb%lim(1,3) = minval( d(ncol-1,:) ) ; tb%lim(2,3) = maxval( d(ncol-1,:) )

  contains

    subroutine Fail ( lnr, reason )
      implicit none
      integer,          intent(in) :: lnr
      character(len=*), intent(in) :: reason
      if ( lnr > 0 ) then
        write(*,'(A)') ERR//trim(path)//' line '//Itoa(lnr)//': '//reason
      else
        write(*,'(A)') ERR//trim(path)//': '//reason
      endif
      error stop 1
    end subroutine Fail

    !> rows (p, B'g, B'c, T, h_w) in file order -> the block structure
    subroutine Build_3d ( r, lr, t3 )
      implicit none
      real(R8),          intent(in)    :: r(:,:)
      integer,           intent(in)    :: lr(:)
      type(bprime_3d_t), intent(inout) :: t3
      ! Local
      integer :: m, i, j, k, np, ns, ip, is, i0, i1, j0, j1, nr
      logical :: up

      m = size(r,2)
      ! p-blocks
      np = 1
      do i = 2, m
        if ( r(1,i) /= r(1,i-1) ) then
          if ( r(1,i) > r(1,i-1) ) call Fail( lr(i), 'p must decrease from block to block' )
          np = np + 1
        endif
      enddo
      ! sub-blocks
      ns = 1
      do i = 2, m
        if ( r(1,i) /= r(1,i-1) ) then
          ns = ns + 1
        elseif ( r(2,i) /= r(2,i-1) ) then
          if ( r(2,i) > r(2,i-1) ) call Fail( lr(i), 'B''g must decrease from sub-block to sub-block' )
          ns = ns + 1
        endif
      enddo
      t3%np = np
      allocate( t3%p(np), t3%g_off(np+1), t3%bg(ns), t3%t_off(ns+1), t3%T(m), t3%bc(m), t3%hw(m) )
      ip = 0 ; is = 0 ; nr = 0
      t3%g_off(1) = 0 ; t3%t_off(1) = 0
      i0 = 1
      do while ( i0 <= m )
        i1 = i0
        do while ( i1 < m )
          if ( r(1,i1+1) /= r(1,i0) ) exit
          i1 = i1 + 1
        enddo
        ip = ip + 1
        t3%p(ip) = r(1,i0)
        if ( t3%p(ip) <= 0.0d0 ) call Fail( lr(i0), 'p must be > 0' )
        j0 = i0
        do while ( j0 <= i1 )
          j1 = j0
          do while ( j1 < i1 )
            if ( r(2,j1+1) /= r(2,j0) ) exit
            j1 = j1 + 1
          enddo
          is = is + 1
          t3%bg(is) = r(2,j0)
          if ( j1 - j0 + 1 < 2 ) call Fail( lr(j0), 'a (p, B''g) sub-table needs at least two T rows' )
          if ( r(4,j1) == r(4,j0) ) call Fail( lr(j0), 'T does not vary in this (p, B''g) sub-table' )
          up = r(4,j1) > r(4,j0)
          do k = j0 + 1, j1
            if ( up ) then
              if ( r(4,k) < r(4,k-1) ) call Fail( lr(k), 'T must not decrease inside an increasing (p, B''g) sub-table' )
            else
              if ( r(4,k) >= r(4,k-1) ) call Fail( lr(k), 'T must decrease strictly inside a decreasing (p, B''g) sub-table' )
            endif
          enddo
          do k = j0, j1
            nr = nr + 1
            if ( up ) then
              j = k
            else
              j = j1 - (k - j0)
            endif
            t3%T(nr)  = r(4,j)
            t3%bc(nr) = r(3,j)
            t3%hw(nr) = r(5,j)
          enddo
          t3%t_off(is+1) = nr
          j0 = j1 + 1
        enddo
        t3%g_off(ip+1) = is
        i0 = i1 + 1
      enddo
    end subroutine Build_3d

  end subroutine Read_Table


  integer function Count_Tokens ( s ) result ( n )
    implicit none
    character(len=*), intent(in) :: s
    integer :: i
    logical :: inword
    n = 0 ; inword = .false.
    do i = 1, len_trim(s)
      if ( s(i:i) == ' ' .or. s(i:i) == char(9) ) then
        inword = .false.
      elseif ( .not. inword ) then
        inword = .true.
        n = n + 1
      endif
    enddo
  end function Count_Tokens


  !> blanks collapsed to one, trimmed (for the axes line)
  function Squeeze ( s ) result ( o )
    implicit none
    character(len=*), intent(in) :: s
    character(len=len(s)) :: o
    integer :: i, k
    logical :: blank
    o = '' ; k = 0 ; blank = .true.
    do i = 1, len_trim(s)
      if ( s(i:i) == ' ' .or. s(i:i) == char(9) ) then
        if ( .not. blank ) then
          k = k + 1 ; o(k:k) = ' '
        endif
        blank = .true.
      else
        k = k + 1 ; o(k:k) = s(i:i)
        blank = .false.
      endif
    enddo
    o = trim(o)
  end function Squeeze


  function Itoa ( i ) result ( s )
    implicit none
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=16) :: b
    write(b,'(I0)') i
    s = trim(b)
  end function Itoa

end module FUSS_Lib_Bprime
