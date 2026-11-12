!> Unit check of Read_Points_Table, the reader of the tables of points (INPUT/properties.dat; Load_Table), against a
!> list-directed read of the same file: every zone holds the rows of the file at the indices 1:N (Ni = N, Nj = Nk = 1), the
!> first column in mesh(1,i,1,1) and the others in vars(1:nvars,i,1,1), bit for bit, and a call that asks for one column
!> more is refused. Arguments: the table and its number of columns after the first.
program table_reader_check
  use iso_fortran_env, only: I8 => int64, R8 => real64
  use Lib_ORION_data,  only: ORION_data
  use FUSS_Load_Table, only: Read_Points_Table
  implicit none
  type(ORION_data)      :: orion
  character(len=4096)   :: file, line
  character(len=64)     :: arg
  real(R8), allocatable :: row(:)
  integer               :: nvars, unit, ios, z, i, rows, bad

  call get_command_argument ( 1, file )
  call get_command_argument ( 2, arg )
  read(arg,*) nvars
  if ( Read_Points_Table ( orion, nvars, trim(file) ) /= 0 ) then
    write(*,'(A)') 'TableReader: FAILED (Read_Points_Table refused '//trim(file)//')'
    error stop 1
  end if

  ! The file line by line: a line that starts with ZONE opens the next zone, a line that starts with a number is a row
  allocate ( row(0:nvars) )
  z = 0; i = 0; rows = 0; bad = 0
  open(newunit=unit, file=trim(file), status='old', action='read')
  do
    read(unit,'(A)',iostat=ios) line
    if ( ios /= 0 ) exit
    line = adjustl(line)
    if ( index(line,'ZONE') == 1 ) then
      if ( z > 0 ) call Check_Zone ( z, i )
      z = z + 1; i = 0
    else if ( scan(line(1:1),'0123456789+-.') == 1 ) then
      i = i + 1; rows = rows + 1
      read(line,*) row
      if ( z < 1 .or. z > size(orion%block) ) cycle
      if ( i > size(orion%block(z)%mesh,2) ) cycle
      if ( transfer(row(0),1_I8) /= transfer(orion%block(z)%mesh(1,i,1,1),1_I8) .or. &
           any( transfer(row(1:),1_I8,nvars) /= transfer(orion%block(z)%vars(:,i,1,1),1_I8,nvars) ) ) bad = bad + 1
    end if
  end do
  close(unit)
  call Check_Zone ( z, i )
  if ( size(orion%block) /= z ) bad = bad + 1
  if ( Read_Points_Table ( orion, nvars + 1, trim(file) ) == 0 ) bad = bad + 1
  write(*,'(A,I0,A,I0,2A)') ' ', z, ' zones, ', rows, ' rows: ', trim(file)
  if ( bad == 0 ) then
    write(*,'(A)') 'TableReader: all checks passed'
  else
    write(*,'(A,I0,A)') 'TableReader: FAILED (', bad, ' rows, zones or refusals differ from the file)'
    error stop 1
  end if

contains

  ! The layout of zone z, of n rows in the file
  subroutine Check_Zone ( z, n )
    integer, intent(in) :: z, n
    if ( z < 1 .or. z > size(orion%block) ) return
    associate ( blk => orion%block(z) )
      if ( blk%Ni /= n .or. blk%Nj /= 1 .or. blk%Nk /= 1 .or. any( lbound(blk%mesh) /= 1 ) .or. &
           any( ubound(blk%mesh) /= [1, n, 1, 1] ) .or. any( lbound(blk%vars) /= 1 ) .or. &
           any( ubound(blk%vars) /= [nvars, n, 1, 1] ) ) bad = bad + 1
    end associate
  end subroutine Check_Zone

end program table_reader_check
