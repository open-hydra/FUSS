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

contains

  !> Reads a time table. Every rank reads the same file, so an error stops every rank together (mpi_stop_all, exit code 1).
  subroutine initialize(self,file)
    use FUSS_Mod_MPI, only: mpi_stop_all
    implicit none
    class(time_series_type), intent(inout) :: self
    character(len=*), intent(in)     :: file
    character(len=256)               :: line
    integer :: ios, i
    integer :: UnitFree

    self%exists = .true.
    self%n = -1
    open(newunit=UnitFree, file = trim(file), status = 'OLD', action = 'READ', form = 'FORMATTED', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') ' [ERROR] time table "'//trim(file)//'" not found'
      call mpi_stop_all('time table not found')
    endif
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
    if (self%n < 1) then
      write(*,'(A)') ' [ERROR] time table "'//trim(file)//'" is empty'
      call mpi_stop_all('empty time table')
    endif
    do i = 1, self%n
      read(UnitFree,*,iostat=ios) self%time(i), self%var(i)
      if (ios /= 0) then
        write(*,'(A,I0)') ' [ERROR] time table "'//trim(file)//'": cannot read row ', i
        call mpi_stop_all('bad time table row')
      endif
      if (i > 1) then
        if (self%time(i) <= self%time(i-1)) then
          write(*,'(A,I0)') ' [ERROR] time table "'//trim(file)//'": times must increase, row ', i
          call mpi_stop_all('non-increasing time table')
        endif
      endif
    enddo
    close(UnitFree)

  end subroutine initialize


  pure function update(self,time) result(varout)
    implicit none
    class(time_series_type), intent(in) :: self
    real(8), intent(in)  :: time
    real(8)              :: varout
    real(8)              :: period, t_
    integer              :: i

    associate( time_ => self%time , var_ => self%var)
    t_ = time
    if (self%periodic) then
      period = time_(self%n) - time_(1)
      if (period > 0d0) t_ = time_(1) + modulo(t_ - time_(1), period)
    endif

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