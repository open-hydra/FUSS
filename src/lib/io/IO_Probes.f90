module FUSS_IO_Probes
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
    real(R8)                         :: tnext = -huge(1d0)   ! next time-based output
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
    use FUSS_Mod_MPI, only: is_local_block, mpi_abort_all
    implicit none
    type(FUSS_domain_type), intent(in)  :: domain
    logical, intent(in)                 :: newrun
    ! Local
    integer             :: i, ii, error
    character(len=clen) :: string(8)

    nprobes = size(obj_io_probes)
    allocate(probe(1:nprobes))

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
      ! gets here, so a failure aborts all ranks: an unopened unit would write to an
      ! arbitrary file.
      if (.not.newrun) then
        open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%file),status='UNKNOWN', &
             position='APPEND',action='WRITE',iostat=error)
      else
        open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%file),status='REPLACE',iostat=error)
      endif
      if (error/=0) call mpi_abort_all('cannot open probe file '//trim(obj_io_probes(i)%file))
  
    enddo

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
    use FUSS_Mod_MPI, only: mpi_abort_all
    use IR_precision, only: str
    implicit none
    type(obj_probe), intent(in)        :: probe
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    integer :: b, d

    b = probe%ilocation(1)
    if (b < 1 .or. b > domain%nb) &
      call mpi_abort_all('probe block '//trim(str(.true.,b))//' does not exist')
    do d = 1, 3
      if (probe%ilocation(1+d) < 1 .or. probe%ilocation(1+d) > domain%blk(b)%dim(d)) &
        call mpi_abort_all('probe index-position '//trim(str(.true.,probe%ilocation(1+d)))// &
                           ' outside block '//trim(str(.true.,b)))
    enddo

  end subroutine Check_Location


  subroutine Assign_Variables(probe, domain)
    use IR_precision
    use FUSS_Mod_MPI, only: mpi_abort_all
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
        call mpi_abort_all('unknown probe variable "'//trim(probe%names(v))//'" (accepted: T)')
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
      if (probe(i)%tnext < -0.5d0*huge(1d0)) &   ! first call (also after a restart): next grid point of dtime
        probe(i)%tnext = ( aint(time/probe(i)%dtime) + 1d0 ) * probe(i)%dtime
      if (mod(iter,probe(i)%diter)==0) then
        write(probe(i)%unit,*) iter, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
      elseif (time >= probe(i)%tnext) then
        write(probe(i)%unit,*) time, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
        probe(i)%tnext = max(probe(i)%tnext, aint(time/probe(i)%dtime)*probe(i)%dtime) + probe(i)%dtime
      endif
    enddo

  end subroutine Write_Probes_Data


end module FUSS_IO_Probes