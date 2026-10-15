module FUSS_Series_Data_m
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  type, public :: time_series_type
    logical              :: exists=.false.
    logical              :: periodic=.false.
    character(len=5)     :: name
    integer              :: n
    real(8), allocatable :: time(:), var(:)
  contains
    procedure, pass(self) :: initialize
    procedure, pass(self) :: update
  end type time_series_type

  !> An n-column table in time (plan 08 section 2.1: the env table of an
  !> ablating wall, columns p_w rhoeUeCH0 h_r q_rad chem hconv T_edge).
  !> Linear in t between rows, evaluated as v_i + w (v_(i+1) - v_i): that form
  !> returns v_i exactly between equal rows, which (1 - w) v_i + w v_(i+1) does
  !> not (it misses by one ulp for 33 of 999 values of w at v = 2800), and a
  !> BC that delegates to another kernel must hand it bit-identical inputs.
  !> CLAMPED to the first row before t_1 and to the last after t_n, as PATO's
  !> constantBoundaryMappingModel does.
  type, public :: time_table_t
    integer :: n = 0                       ! columns after t
    integer :: nrow = 0
    real(8), allocatable :: time(:)        ! (nrow), strictly increasing
    real(8), allocatable :: var(:,:)       ! (n, nrow)
  contains
    procedure, pass(self) :: load => time_table_load
    procedure, pass(self) :: eval => time_table_eval
  end type time_table_t

contains

  !> Read `t v_1 .. v_n` rows; `//` lines are comments, blank lines skipped.
  !> At least two rows, t strictly increasing; any other line, a short row or
  !> a non-finite value stops the run naming the file and the line.
  subroutine time_table_load ( self, file, n )
    implicit none
    class(time_table_t), intent(inout) :: self
    character(len=*),    intent(in)    :: file
    integer,             intent(in)    :: n
    ! Local
    character(len=1024) :: line, t
    integer :: u, ios, ln, k, i
    real(8) :: v(n+1)
    real(8), allocatable :: rows(:,:)

    open( newunit=u, file=trim(file), status='old', action='read', iostat=ios )
    if ( ios /= 0 ) then
      write(*,'(A)') ' [ERROR] env table '//trim(file)//': cannot open'
      error stop 1
    endif
    allocate( rows(n+1, 64) )
    k = 0 ; ln = 0
    do
      read( u, '(A)', iostat=ios ) line
      if ( ios /= 0 ) exit
      ln = ln + 1
      t = adjustl(line)
      if ( len_trim(t) == 0 ) cycle
      if ( t(1:2) == '//' ) cycle
      read( t, *, iostat=ios ) v
      if ( ios /= 0 ) call Fail( 'expected '//Itoa(n+1)//' numbers (t and '//Itoa(n)//' columns)' )
      do i = 1, n + 1
        if ( v(i) /= v(i) .or. abs(v(i)) > huge(1.0d0) ) call Fail( 'value is not finite' )
      enddo
      if ( k > 0 ) then
        if ( v(1) <= rows(1,k) ) call Fail( 't must increase strictly from row to row' )
      endif
      if ( k == size(rows,2) ) call Grow()
      k = k + 1
      rows(:,k) = v
    enddo
    close( u )
    if ( k < 2 ) then
      write(*,'(A)') ' [ERROR] env table '//trim(file)//': fewer than two rows'
      error stop 1
    endif
    self%n = n
    self%nrow = k
    if ( allocated(self%time) ) deallocate( self%time )
    if ( allocated(self%var) )  deallocate( self%var )
    allocate( self%time(k), self%var(n,k) )
    self%time = rows(1,1:k)
    self%var  = rows(2:n+1,1:k)

  contains

    subroutine Fail ( reason )
      character(len=*), intent(in) :: reason
      write(*,'(A)') ' [ERROR] env table '//trim(file)//' line '//Itoa(ln)//': '//reason
      error stop 1
    end subroutine Fail

    subroutine Grow ()
      real(8), allocatable :: tmp(:,:)
      allocate( tmp(n+1, 2*size(rows,2)) )
      tmp(:,1:size(rows,2)) = rows
      call move_alloc( tmp, rows )
    end subroutine Grow

    function Itoa ( i ) result ( s )
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=12) :: b
      write(b,'(I0)') i
      s = trim(b)
    end function Itoa

  end subroutine time_table_load


  pure function time_table_eval ( self, t ) result ( v )
    implicit none
    class(time_table_t), intent(in) :: self
    real(8),             intent(in) :: t
    real(8) :: v(self%n)
    ! Local
    integer :: i
    real(8) :: w

    if ( t <= self%time(1) ) then
      v = self%var(:,1)
      return
    endif
    if ( t >= self%time(self%nrow) ) then
      v = self%var(:,self%nrow)
      return
    endif
    i = 1
    do while ( self%time(i+1) < t )
      i = i + 1
    enddo
    if ( self%time(i+1) == t ) then
      v = self%var(:,i+1)
      return
    endif
    w = ( t - self%time(i) ) / ( self%time(i+1) - self%time(i) )
    v = self%var(:,i) + w * ( self%var(:,i+1) - self%var(:,i) )

  end function time_table_eval


  subroutine initialize(self,file,bar)
    implicit none
    class(time_series_type), intent(inout) :: self
    character(len=*), intent(in)     :: file
    logical, intent(in)              :: bar
    character(len=256)               :: line
    integer :: ios, i
    integer :: UnitFree

    self%exists = .true.
    self%n = -1
    open(newunit=UnitFree, file = trim(file), status = 'OLD', action = 'READ', form = 'FORMATTED', iostat=ios)
    ios = 0
    do while (ios==0)
      read(UnitFree,'(A)',iostat=ios) line
      self%n = self%n+1
      if (index(line, 'periodic') > 0) then
        self%periodic = .true.
        self%n = self%n -1 ! Do not count this one
      endif
    enddo
    rewind(UnitFree)
    allocate(self%time(1:self%n))
    allocate(self%var(1:self%n))
    if (self%periodic) read(UnitFree,*) 
    do i = 1, self%n
      read(UnitFree,*) self%time(i), self%var(i)
    enddo
    close(UnitFree)
    
    if (bar) self%var = self%var

  end subroutine initialize


  pure function update(self,time) result(varout)
    implicit none
    class(time_series_type), intent(in) :: self
    real(8), intent(in)  :: time
    real(8)              :: varout
    real(8)              :: period, t_
    integer              :: i

    associate( time_ => self%time , var_ => self%var)
    varout = var_(size(var_))    ! a time no branch below takes (NaN): the last value
    t_ = time
    if (self%periodic) then
      period = time_(self%n) - time_(1)
      if (period > 0d0) t_ = time_(1) + modulo(t_ - time_(1), period)
    endif

    ! At or before the first time: the first value (before it the result was
    ! left undefined), as after the last time the last value.
    if (t_ <= time_(1)) then
        varout = var_(1)
        return
    endif

    do i = 2, size(time_)
      if (t_>time_(i-1) .and. t_<=time_(i)) then
        varout = (var_(i)-var_(i-1))/(time_(i)-time_(i-1))*(t_-time_(i-1))+var_(i-1)
        return
      endif
    enddo

    if (t_>time_(self%n)) varout = var_(size(var_))

    endassociate

  end function update

end module FUSS_Series_Data_m