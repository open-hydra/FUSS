module FUSS_IO_Probes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Parameters_m
  use FUSS_Base_Types_m, only: FUSS_vector_3D_type
  use FUSS_Config_Types_m, only: obj_io_probes
  use FUSS_Advanced_Types_m, only: FUSS_domain_type

  implicit none
  private
  public :: Setup_Probes, Write_Probes_Data, Close_Probes, Restart_Probes, Drop_Probe_Rows, Drop_Rows_After

  type :: real_ptr
    real(R8), pointer :: p
  end type real_ptr

  type :: obj_probe
    type(FUSS_vector_3D_type)        :: location
    integer, dimension(4)            :: ilocation
    integer                          :: nvar
    character(len=clen), allocatable :: names(:)
    type(real_ptr), allocatable      :: variables(:)
    real(R8), pointer                :: P
    real(R8)                         :: dtime
    integer                          :: ntime = -1           ! next multiple of dtime with a time-based row (-1: not set yet)
    integer                          :: diter
    integer                          :: unit
  contains
    private 
    procedure, pass(self) :: Place
  end type obj_probe

  integer, public :: nprobes
  type(obj_probe), allocatable, target :: probe(:)

contains

  subroutine Setup_Probes( domain, newrun )
    use IR_precision
    use strings, only: parse
    use FUSS_Mod_MPI, only: is_local_block, mpi_stop_all, mpi_stop_if_any
    implicit none
    type(FUSS_domain_type), intent(in)  :: domain
    logical, intent(in)                 :: newrun
    ! Local
    integer             :: i, ii, error
    character(len=clen) :: string(8)
    logical             :: open_failed
    character(len=hlen) :: open_detail

    nprobes = size(obj_io_probes)
    allocate(probe(1:nprobes))
    open_failed = .false.
    open_detail = ''

    do i = 1, nprobes

      obj_io_probes(i)%file = 'OUTPUT/'//trim(obj_io_probes(i)%file)//'.txt'

      ! Vars name
      string = ''
      probe(i)%nvar = 0
      call parse(obj_io_probes(i)%varnames,' ',string)
      do ii = 1, size(string)
        if (string(ii) /= '') then
          probe(i)%nvar = probe(i)%nvar + 1
        endif
      enddo
      allocate(probe(i)%names(1:probe(i)%nvar))
      do ii = 1, probe(i)%nvar
        probe(i)%names(ii) = trim(string(ii))
      enddo

      ! The variables are known names: every rank checks them, whether it owns the probe or not (a rank that does not own
      ! the probe would otherwise go on while the owner stops)
      do ii = 1, probe(i)%nvar
        if (trim(probe(i)%names(ii)) /= 'T') &
          call mpi_stop_all('unknown probe variable "'//trim(probe(i)%names(ii))//'" (accepted: T)')
      enddo

      ! Frequency
      probe(i)%dtime = obj_io_probes(i)%dtime
      probe(i)%diter = obj_io_probes(i)%diter

      ! Location
      probe(i)%location%c = obj_io_probes(i)%loc
      probe(i)%ilocation  = obj_io_probes(i)%iloc

      if (sum(obj_io_probes(i)%iloc)==0) call probe(i)%Place(domain)
      call Check_Location(probe(i), domain)

      ! Only the rank owning this probe's block sets up pointers and opens the file
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle

      call Assign_Variables(probe(i), domain)

      ! A restart appends to the probe file, creating it if missing (e.g. probes added
      ! or files moved away); a new run replaces it. Only the rank that owns the probe
      ! gets here: a failure is reduced over the ranks after the loop, and stops them all (an
      ! unopened unit would write to an arbitrary file).
      if (.not.newrun) then
        open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%file),status='UNKNOWN', &
             position='APPEND',action='WRITE',iostat=error)
      else
        open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%file),status='REPLACE',iostat=error)
      endif
      if (error/=0 .and. .not. open_failed) then
        open_failed = .true.
        open_detail = '[ERROR] cannot open probe file '//trim(obj_io_probes(i)%file)
      endif
  
    enddo

    call mpi_stop_if_any(open_failed, 'cannot open a probe file (see above)', open_detail)

  end subroutine Setup_Probes


  !> The cell whose centre is nearest to the probe location (needs the nodes of every block)
  subroutine Place (self, domain)
    implicit none
    class(obj_probe)                   :: self
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    integer :: i, j, k, b
    real(8) :: d0, d, c(3)

    d0 = huge(1d0)
    do b = 1, domain%nb
      associate ( nd => domain%blk(b)%node )
      do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        c = 0.125d0 * ( nd(i-1,j-1,k-1)%c + nd(i,j-1,k-1)%c + nd(i-1,j,k-1)%c + nd(i,j,k-1)%c &
                      + nd(i-1,j-1,k  )%c + nd(i,j-1,k  )%c + nd(i-1,j,k  )%c + nd(i,j,k  )%c )
        d = norm2(self%location%c - c)
        if (d<d0) then
          d0 = d
          self%ilocation = [b , i, j, k]
        endif
      enddo; enddo; enddo
      end associate
    enddo

  end subroutine Place


  !> index-position must name an existing cell
  subroutine Check_Location(probe, domain)
    use FUSS_Mod_MPI, only: mpi_stop_all
    use IR_precision, only: str
    implicit none
    type(obj_probe), intent(in)        :: probe
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    integer :: b, d

    b = probe%ilocation(1)
    if (b < 1 .or. b > domain%nb) &
      call mpi_stop_all('probe block '//trim(str(.true.,b))//' does not exist')
    do d = 1, 3
      if (probe%ilocation(1+d) < 1 .or. probe%ilocation(1+d) > domain%blk(b)%dim(d)) &
        call mpi_stop_all('probe index-position '//trim(str(.true.,probe%ilocation(1+d)))// &
                           ' outside block '//trim(str(.true.,b)))
    enddo

  end subroutine Check_Location


  subroutine Assign_Variables(probe, domain)
    use IR_precision
    use FUSS_Mod_MPI, only: mpi_stop_all
    implicit none
    type(obj_probe), intent(inout), target     :: probe
    type(FUSS_domain_type), intent(in), target :: domain
    ! Local
    integer :: v,s,b,i,j,k

    b = probe%ilocation(1)
    i = probe%ilocation(2)
    j = probe%ilocation(3)
    k = probe%ilocation(4)

    probe%P => domain%blk(b)%T(i,j,k)
    allocate(probe%variables(1:probe%nvar))

    do v = 1, probe%nvar
      nullify(probe%variables(v)%p)
      if (probe%names(v)=='T') probe%variables(v)%p => probe%P
      if (.not. associated(probe%variables(v)%p)) &
        call mpi_stop_all('unknown probe variable "'//trim(probe%names(v))//'" (accepted: T)')
    enddo

  end subroutine Assign_Variables


  !> Time counter of the time-based probe rows at a restart. Write_Probes_Data writes a time-based row when time >=
  !> dtime*ntime and then counts it (the next multiple that the time has not reached), so the rows fall on the multiples of
  !> dtime, at most one per step. A restart that goes on in time (resume) takes for ntime the first multiple its time has
  !> not reached, computed with the same product as Write_Probes_Data so that an exact multiple counts as written;
  !> otherwise, and in a steady run (time < 0), the counter is set by the first call, as in a new run.
  subroutine Restart_Probes ( time, resume )
    implicit none
    real(R8), intent(in) :: time
    logical,  intent(in) :: resume
    ! Local
    integer :: i

    do i = 1, nprobes
      probe(i)%ntime = -1
      if (.not. resume .or. time < 0.d0) cycle
      if (time / probe(i)%dtime >= 0.5_R8 * real(huge(1), R8)) cycle
      probe(i)%ntime = max(0, int(time / probe(i)%dtime) - 1)
      do while ( time >= probe(i)%dtime * probe(i)%ntime )
        probe(i)%ntime = probe(i)%ntime + 1
      enddo
    enddo

  end subroutine Restart_Probes


  !> Drop from the probe files the rows written after the solution a restart reads. For a restart that goes on from the
  !> counters of its solution (see Drop_Rows_After). Each rank handles the probes it owns: it closes the file, drops the
  !> rows, and opens it again at its end, as Setup_Probes does.
  subroutine Drop_Probe_Rows ( iteration, time )
    use FUSS_Mod_MPI, only: is_local_block, mpi_stop_if_any
    implicit none
    integer,  intent(in) :: iteration
    real(R8), intent(in) :: time
    ! Local
    integer :: i, error
    logical :: failed
    character(len=hlen) :: detail

    failed = .false.
    detail = ''
    do i = 1, nprobes
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
      close(probe(i)%unit)
      call Drop_Rows_After ( trim(obj_io_probes(i)%file), 1 + probe(i)%nvar, iteration, time )
      open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%file),status='UNKNOWN', &
           position='APPEND',action='WRITE',iostat=error)
      if (error/=0 .and. .not. failed) then
        failed = .true.
        detail = '[ERROR] cannot open probe file '//trim(obj_io_probes(i)%file)
      endif
    enddo
    call mpi_stop_if_any(failed, 'cannot open a probe file (see above)', detail)

  end subroutine Drop_Probe_Rows


  !> Drop the rows of a history file written after the solution a restart reads. residual-history.dat and the probe files
  !> get a row at the iterations, or at the times, that call for one. A run that went on past the solution a restart
  !> reads (it stopped without writing its final solution, and the restart reads an earlier numbered one) left in them rows
  !> that the restart writes again. The rows kept are those before the first one whose first value is an iteration
  !> greater than iteration, or a time later than time; that row and all the rows after it are dropped, so the file goes
  !> on as that of the run in one go. A row holds ncol values; written list-directed it may take more than one line, and a
  !> row starts on a new line: the values are counted to find where each row starts. A time is compared as printed, to
  !> 1e-12 relative: the row of the step of the solution may print its time rounded up in the last digit. A last row with
  !> fewer than ncol values (a file cut in the middle of a row) is dropped too. A line that starts with '#' (a header) is
  !> kept and holds no values. The file is truncated and closed, and ends with a newline (End_With_Newline); a missing
  !> file is left alone.
  subroutine Drop_Rows_After ( file, ncol, iteration, time )
    implicit none
    character(len=*), intent(in) :: file
    integer,          intent(in) :: ncol, iteration
    real(R8),         intent(in) :: time
    ! Local
    character(len=32768) :: line
    integer              :: unit, ios, nkeep, nrow, nval, k, k1, k2, it
    real(R8)             :: t
    logical              :: later

    open(newunit=unit,file=file,status='OLD',action='READWRITE',iostat=ios)
    if (ios /= 0) return
    nkeep = 0
    nrow  = 0
    nval  = 0
    do
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      ! A header line is kept and holds no values
      if (len_trim(line) > 0) then
        if (line(verify(line, ' '):verify(line, ' ')) == '#') then
          nkeep = nkeep + 1
          cycle
        endif
      endif
      if (mod(nval, ncol) == 0 .and. len_trim(line) > 0) then
        ! The line starts a row: its first value, line(k1:k2), is an iteration or a time
        k1 = verify(line, ' ')
        k2 = k1 + index(line(k1:), ' ') - 2
        later = .false.
        if (scan(line(k1:k2), '.eEdD') > 0) then
          read(line(k1:k2),*,iostat=ios) t
          if (ios == 0) later = t > time + abs(time) * 1.0e-12_R8
        else
          read(line(k1:k2),*,iostat=ios) it
          if (ios == 0) later = it > iteration
        endif
        if (later) exit
        ! Lines kept before the row that starts here
        nrow = nkeep
      endif
      ! Values on the line
      do k = 1, len_trim(line)
        if (line(k:k) /= ' ' .and. (k == 1 .or. line(max(k-1,1):max(k-1,1)) == ' ')) nval = nval + 1
      enddo
      nkeep = nkeep + 1
    enddo
    ! A last row with fewer than ncol values is cut: drop it from its first line
    if (mod(nval, ncol) /= 0) nkeep = nrow
    rewind(unit)
    do k = 1, nkeep
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
    enddo
    endfile(unit)
    close(unit)
    call End_With_Newline ( file )

  end subroutine Drop_Rows_After


  !> End a file with a newline if its last line has none. A history file cut by hand may end without a newline; the rows
  !> that a restart then writes in append mode would go on that line. Read and written as a stream of bytes.
  subroutine End_With_Newline ( file )
    use iso_fortran_env, only: int64
    implicit none
    character(len=*), intent(in) :: file
    ! Local
    integer        :: unit, ios
    integer(int64) :: nbytes
    character      :: last

    open(newunit=unit,file=file,access='STREAM',form='UNFORMATTED',status='OLD',action='READWRITE',iostat=ios)
    if (ios /= 0) return
    inquire(unit=unit,size=nbytes)
    if (nbytes > 0) then
      read(unit,pos=nbytes,iostat=ios) last
      if (ios == 0 .and. last /= new_line('a')) write(unit,pos=nbytes+1,iostat=ios) new_line('a')
    endif
    close(unit)

  end subroutine End_With_Newline


  !> Write the probe rows due at this iteration (each rank the probes it owns), and pass each row to the file system when it
  !> is written (flush): a run that stops before its end (killed, stopped at the time limit of a queue, aborted) leaves in
  !> the files every row it wrote, as a restart from its last solution expects (Drop_Rows_After); the rows held in the buffer
  !> of the unit would be lost with the process. A time-based row is written when the time
  !> has reached the next multiple of dtime, and at most one per step: the counter ntime then goes to the multiple after
  !> the one the time is in. The first call (also after a restart that gives no counter) sets it to the multiple after the
  !> time of the call.
  subroutine Write_Probes_Data( iter, time )
    use FUSS_Mod_MPI, only: is_local_block

    implicit none
    integer, intent(in) :: iter
    real(8), intent(in) :: time
    ! Local
    integer :: i, v

    do i = 1, nprobes
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
      if (probe(i)%ntime < 0) then   ! first call: the next multiple of dtime
        probe(i)%ntime = 1
        if (time / probe(i)%dtime < 0.5d0 * real(huge(1), R8)) probe(i)%ntime = max(1, int(time / probe(i)%dtime) + 1)
      endif
      if (mod(iter,probe(i)%diter)==0) then
        write(probe(i)%unit,*) iter, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
        flush(probe(i)%unit)
      elseif (time >= probe(i)%dtime*probe(i)%ntime) then
        write(probe(i)%unit,*) time, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
        flush(probe(i)%unit)
        probe(i)%ntime = probe(i)%ntime + 1
        if (time / probe(i)%dtime < 0.5d0 * real(huge(1), R8)) probe(i)%ntime = max(probe(i)%ntime, int(time / probe(i)%dtime) + 1)
      endif
    enddo

  end subroutine Write_Probes_Data


  !> Close the probe files at the end of the run: each rank closes the files of the probes it owns, which Setup_Probes or
  !> Drop_Probe_Rows opened, before the ranks finalize MPI (a unit left open is closed by the Fortran run-time library only
  !> when the process ends, after MPI_FINALIZE).
  subroutine Close_Probes ()
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    integer :: i

    do i = 1, nprobes
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
      close(probe(i)%unit)
    enddo

  end subroutine Close_Probes


end module FUSS_IO_Probes