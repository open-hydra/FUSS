module FUSS_Mod_Diagnostic
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Parameters_m

  implicit none
  private
  public :: Compute_Residual, Write_Diagnostic
  character(len=llen), private :: Dvarnames = '"E" "dt"'
  
contains


  !> The sum over the cells of the block is taken in a fixed order, the same for any number of threads: each thread sums
  !> whole i lines, and the sums of the lines are added in turn.
  subroutine Compute_Residual ( new, old, dt, n, average, total )
    use FUSS_Global_m
    implicit none
    integer, intent(in)     :: n(3)
    real(R8), intent(in)    :: new ( 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc )
    real(R8), intent(in)    :: old ( 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc )
    real(R8), intent(in)    :: dt ( n(1), n(2), n(3) )
    real(R8), intent(out)   :: average
    real(R8), intent(inout) :: total
    ! Local
    integer :: i, j, k
    real(R8) :: resn, residuo
    real(R8), allocatable :: line(:,:)   ! sum over each i line

    allocate ( line(n(2), n(3)) )

    !$omp parallel do private ( i, resn ) collapse (2)
    do k = 1, n(3)
    do j = 1, n(2)
      line(j,k) = 0d0
      do i = 1, n(1)

        resn = abs ( new(i,j,k) - old(i,j,k) ) ! delta(T)
        line(j,k) = line(j,k) + resn*resn  ! sum ( r(i)**2 )

      enddo
    enddo; enddo
    !$omp end parallel do

    residuo = 0d0
    do k = 1, n(3)
    do j = 1, n(2)
      residuo = residuo + line(j,k)
    enddo; enddo

    ! Average residual of this block
    average = sqrt ( residuo / float ( n(1)*n(2)*n(3) ) )
    ! Overall residual
    total = total + residuo
    deallocate ( line )

  end subroutine Compute_Residual


  !> Update the orion-field data and write them accordingly to the chosen format (vtk,tecplot)
  subroutine Write_Diagnostic( domain, IOfield, file )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_io
    use FUSS_Global_m
    use Lib_ORION_data
    use Lib_VTK
    use Lib_Tecplot
    use FUSS_Mod_MPI, only: mpi_is_root
    use FUSS_Mod_GhostExchange, only: gather_diagnostic_to_root, mpi_io_barrier
    use strings, only: parse
    implicit none
    type(FUSS_domain_type), intent(inout)  :: domain
    type(ORION_data), intent(inout)        :: IOfield
    character(llen), intent(in)            :: file
    ! Local
    character(len=llen) :: path
    character(len=llen) :: localpath_vtk
    character(len=llen) :: vtknames
    character(len=llen) :: name
    integer             :: E_IO, b, i, j, k
    character(len=clen) :: format(2), extension

    ! Gather R, dtlocal, beta from all ranks to root (collective)
    call gather_diagnostic_to_root(domain)

    if (mpi_is_root) then
      path = 'OUTPUT/'

      call parse(obj_io%sol_format,' ', format)
      
      ! Update IOfield variables with domain residuals
      do b = 1, size(IOfield%block)
        IOfield%block(b)%vars(1,:,:,:) = domain%blk(b)%r(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        ! Auxiliary variables: dt
        do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
          IOfield%block(b)%vars(1+1,i,j,k) = domain%blk(b)%dtlocal(i,j,k)
        enddo; enddo; enddo
      enddo

      ! Write the IOfield accordingly to the solution format
      select case(trim(format(1)))
      case('vtk')
        IOfield%vtk%format = trim(format(2))
        ! As the field files (IO_Solution): the .vtm lists the block files relative to its own folder (no leading '/' in
        ! the name), and the time is the solution time, minus the global iteration count in a steady run
        name = file
        if ( name(1:1) == '/' ) name = name(2:)
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        vtknames = unquoted(Dvarnames)
        E_IO = vtk_write_structured_multiblock(orion=IOfield,vtspath=trim(localpath_vtk)//trim(name), &
                                               vtmpath=trim(path)//trim(name),varnames=vtknames,time=IOfield%solutiontime)
      case('tecplot')
        IOfield%tec%format = trim(format(2))
        ! TecIO refuses a binary stream unless the name ends in .plt/.szplt, so
        ! the extension has to follow the format (as MOSE's writer already does).
        if (format(2)=='binary') then
          extension = '.szplt'
        else
          extension = '.tec'
        end if
        E_IO = tec_write_structured_multiblock(Nvars=1+1,orion=IOfield,varnames=Dvarnames,filename=trim(path)//trim(file)//trim(extension))
      end select
    end if

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine Write_Diagnostic

end module FUSS_Mod_Diagnostic
