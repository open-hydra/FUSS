module FUSS_IO_Probes
  use FUSS_Mod_MPI, only: fuss_abort
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Parameters_m
  use FUSS_Base_Types_m, only: FUSS_vector_3D_type
  use FUSS_Config_Types_m, only: obj_io_probes
  use FUSS_Advanced_Types_m, only: FUSS_domain_type

  implicit none
  private
  public :: Setup_Probes, Write_Probes_Data

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
    integer                          :: ntime
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
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(in)  :: domain
    logical, intent(in)                 :: newrun
    ! Local
    integer             :: i, ii, error, ntok
    character(len=clen), allocatable :: string(:)

    nprobes = size(obj_io_probes)
    allocate(probe(1:nprobes))

    do i = 1, nprobes

      ! The path has its own component: written into the 16-character name it
      ! was cut (`d01_lo` wrote OUTPUT/d01_lo.tx).
      obj_io_probes(i)%path = 'OUTPUT/'//trim(obj_io_probes(i)%file)//'.txt'

      ! Vars name. parse writes one element per token without a bound check, so
      ! the array is sized from the count (a fixed 8 let a ninth name write past it).
      ntok = Count_Tokens(obj_io_probes(i)%varnames)
      if ( allocated(string) ) deallocate(string)
      allocate(string(max(ntok, 1)))
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

      ! Frequency
      probe(i)%dtime = obj_io_probes(i)%dtime
      probe(i)%diter = obj_io_probes(i)%diter
      ! Output counter for the time-based trigger in Write_Probes_Data.
      ! Without this the component is undefined: dtime*ntime evaluates to
      ! garbage, the "time >= dtime*ntime" test never fires and the probe file
      ! is written empty (observed on test/transient/oscillating_Rod).
      probe(i)%ntime = 0

      ! Location
      probe(i)%location%c = obj_io_probes(i)%loc
      probe(i)%ilocation  = obj_io_probes(i)%iloc

      if (sum(obj_io_probes(i)%iloc)==0) call probe(i)%Place(domain)

      ! Only the rank owning this probe's block sets up pointers and opens the file
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle

      call Assign_Variables(probe(i), domain)

      ! A restart appends to the probe file, creating it if missing (e.g. probes added
      ! or files moved away); a new run replaces it. Only the rank that owns the probe
      ! gets here, so a failure aborts all ranks: an unopened unit would write to an
      ! arbitrary file.
      if (.not.newrun) then
        open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%path),status='UNKNOWN', &
             position='APPEND',action='WRITE',iostat=error)
      else
        open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%path),status='REPLACE',iostat=error)
      endif
      if (error/=0) then
        write(*,'(A)') ' [ERROR] cannot open probe file '//trim(obj_io_probes(i)%path)
        call fuss_abort ( 1 )
      endif
  
    enddo

  end subroutine Setup_Probes


  !> Number of blank-separated tokens of a string.
  pure integer function Count_Tokens ( text ) result ( n )
    implicit none
    character(len=*), intent(in) :: text
    integer :: k
    logical :: inside
    n = 0
    inside = .false.
    do k = 1, len_trim(text)
      if ( text(k:k) /= ' ' .and. text(k:k) /= char(9) ) then
        if ( .not. inside ) n = n + 1
        inside = .true.
      else
        inside = .false.
      endif
    enddo
  end function Count_Tokens


  subroutine Place (self, domain)
    implicit none
    class(obj_probe)                   :: self
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    integer :: i, j, k, b
    real(8) :: d0, d

    d0 = huge(1d0)
    do b = 1, domain%nb
      do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        d = norm2(self%location%c-domain%blk(b)%node(i,j,k)%c)
        if (d<d0) then
          d0 = d
          self%ilocation = [b , i, j, k]
        endif
      enddo; enddo; enddo
    enddo

  end subroutine Place


  subroutine Assign_Variables(probe, domain)
    use IR_precision
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
      select case ( trim(probe%names(v)) )
        case ('T')
          probe%variables(v)%p => probe%P
        case ('p_pg')
          ! the gas pressure of the type-2 closure (plan 07)
          if ( .not. allocated(domain%blk(b)%p_pg) ) then
            write(*,'(A)') ' [ERROR] probe variable p_pg needs [FUSS-Pyrolysis] gas-closure = type2'
            call fuss_abort ( 1 )
          endif
          probe%variables(v)%p => domain%blk(b)%p_pg(i,j,k)
        case default
          write(*,'(A)') ' [ERROR] unknown probe variable "'//trim(probe%names(v))//'" (known: T, p_pg)'
          call fuss_abort ( 1 )
      end select
    enddo

  end subroutine Assign_Variables


  subroutine Write_Probes_Data( iter, time )
    use FUSS_Mod_MPI, only: is_local_block

    implicit none
    integer, intent(in) :: iter
    real(8), intent(in) :: time
    ! Local
    integer :: i, v

    do i = 1, nprobes
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
      if (mod(iter,probe(i)%diter)==0) then
        write(probe(i)%unit,*) iter, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
      elseif (time >= probe(i)%dtime*probe(i)%ntime) then
        write(probe(i)%unit,*) time, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
        ! Advance to the next output slot. The max() makes this self-correcting
        ! on restart, where ntime starts at 0 but time does not -- otherwise the
        ! counter would have to catch up one call at a time, duplicating rows.
        probe(i)%ntime = max( probe(i)%ntime + 1, int(time/probe(i)%dtime) + 1 )
      endif
    enddo

  end subroutine Write_Probes_Data


end module FUSS_IO_Probes