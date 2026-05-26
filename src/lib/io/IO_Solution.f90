module FUSS_IO_Solution

  implicit none

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(r_solution_if), pointer, public :: Read_IC
  procedure(w_solution_if), pointer, public :: Write_solution

  !> Abstract interface relative to the finite-rate reactions source procedure
  abstract interface
    subroutine r_solution_if ( IOfield )
      use Lib_ORION_data
      implicit none
      type(ORION_data), intent(inout) :: IOfield
    end subroutine r_solution_if

    subroutine w_solution_if ( domain, IOfield, file )
      use FUSS_Advanced_Types_m
      use FUSS_Config_Types_m, only: obj_io
      use FUSS_Global_m
      use FUSS_Parameters_m
      use IR_Precision
      use Lib_ORION_data
      use Lib_VTK
      use Lib_Tecplot
      use FUSS_Lib_Solid
      implicit none
      type(FUSS_domain_type), intent(inout) :: domain
      type(ORION_data), intent(inout)    :: IOfield
      character(llen), intent(in)        :: file
    end subroutine w_solution_if
  end interface

contains

  subroutine Setup_Input_Solution()
    use FUSS_Config_Types_m, only: obj_sim_param, obj_io
    use FUSS_Global_m
    use FUSS_Parameters_m
    use IR_Precision
    implicit none
    logical         :: present
    integer         :: i
    character(llen) :: try
    character(6)    :: extension

    obj_io%ini_format=obj_io%sol_format

    if (obj_sim_param%newrun) then
      if (index(obj_io%ini_format,'vtk')>0) then
        extension = '.vtm'
        Read_IC  => Read_vtk_tec
      else
        if (index(obj_io%ini_format,'ascii')>0) then
          extension = '.tec'
        else
          extension = '.szplt'
        endif
        Read_IC  => Read_vtk_tec
      endif
      obj_io%nameinit   = 'INPUT/'//trim(FUSS_phase_prefix)//'ic'//extension
      obj_io%namesource = 'INPUT/'//trim(FUSS_phase_prefix)//'st'//extension
      inquire(file=obj_io%nameinit, exist=present)

    else

      if (index(obj_io%sol_format,'vtk')>0) then
        extension = '.vtm'
        Read_IC => Read_vtk_tec
      else
        if (index(obj_io%sol_format,'ascii')>0) then
          extension = '.tec'
        else
          extension = '.szplt'
        endif
        Read_IC => Read_vtk_tec
      end if
      obj_io%nameinit = 'OUTPUT/'//trim(FUSS_phase_prefix)//'field'//extension
      inquire(file=obj_io%nameinit, exist=present)

      if (.not.present) then
        i = 0
        do
          i = i+1
          try = 'OUTPUT/'//trim(FUSS_phase_prefix)//'field'//trim(str(.true.,i))//extension
          inquire(file=try,exist=present)
          if (present) then
            obj_io%nameinit = try
          else
            exit
          endif
        enddo
      endif

    endif

  end subroutine Setup_Input_Solution

  subroutine Read_vtk_tec ( IOfield )
    use FUSS_Config_Types_m, only: obj_sim_param, obj_io
    use FUSS_Parameters_m
    use Lib_ORION_data
    use Lib_VTK
    use Lib_Tecplot
    use strings, only: parse
    implicit none
    type(ORION_data), intent(inout) :: IOfield
    ! Local
    type(ORION_data) :: IOinit, IOsource
    integer          :: error_init, error_source, error_dim, b, nblocks
    character(clen)  :: format(2)

    error_init   = 0
    error_source = 0
    error_dim    = 0

    if ( obj_sim_param%newrun ) then
      call parse(obj_io%ini_format,' ', format)
    else
      call parse(obj_io%sol_format,' ', format)
    endif

    obj_io%IOtime = 0.d0
    
    select case(trim(format(1)))
    case('tecplot')
      IOfield%tec%format  = trim(format(2))
      IOinit%tec%format   = trim(format(2))
      IOsource%tec%format = trim(format(2))
      error_init   = tec_read_structured_multiblock(orion=IOinit,  filename=trim(obj_io%nameinit))
      error_source = tec_read_structured_multiblock(orion=IOsource,filename=trim(obj_io%namesource))
    case('vtk')
      IOfield%vtk%format  = trim(format(2))
      IOinit%vtk%format   = trim(format(2))
      IOsource%vtk%format = trim(format(2))
      error_init   = vtk_read_structured_multiblock(orion=IOinit,  vtmpath=obj_io%nameinit(1:len(trim(obj_io%nameinit))-4),vtspath='INPUT/vtk/field',time=obj_io%IOtime)
      error_source = vtk_read_structured_multiblock(orion=IOsource,vtmpath=obj_io%namesource(1:len(trim(obj_io%namesource))-4),vtspath='INPUT/vtk/field',time=obj_io%IOtime)
    end select

    if (error_init/=0) then
      obj_io%error_message = "[ERROR] reading input file "//trim(obj_io%nameinit)
      return
    endif

    ! Allocate IOfield
    nblocks = size(IOinit%block)
    allocate(IOfield%block(1:nblocks))
    do b = 1, nblocks
      allocate(IOfield%block(b)%mesh(1:3,0:IOinit%block(b)%Ni,0:IOinit%block(b)%Nj,0:IOinit%block(b)%Nk))
      allocate(IOfield%block(b)%vars(1:3,1:IOinit%block(b)%Ni,1:IOinit%block(b)%Nj,1:IOinit%block(b)%Nk))
    enddo

    ! Check dimensions IOinit and IOsource are compatible
    if (error_source == 0) then
      if (size(IOsource%block) /= size(IOinit%block)) error_dim = 1
      do b = 1, nblocks
        if (IOsource%block(b)%Ni /= IOinit%block(b)%Ni) error_dim = 1
        if (IOsource%block(b)%Nj /= IOinit%block(b)%Nj) error_dim = 1
        if (IOsource%block(b)%Nk /= IOinit%block(b)%Nk) error_dim = 1
      enddo
    endif
    if (error_dim == 1) then
      obj_io%error_message = "[ERROR] source file dimensions not compatible with init file"
      return
    endif

    ! Fill IOfield
    IOfield%solutiontime = IOinit%solutiontime
    do b = 1, nblocks
      IOfield%block(b)%Ni = IOinit%block(b)%Ni
      IOfield%block(b)%Nj = IOinit%block(b)%Nj
      IOfield%block(b)%Nk = IOinit%block(b)%Nk
      IOfield%block(b)%mesh = IOinit%block(b)%mesh
      IOfield%block(b)%vars(1:2,:,:,:) = IOinit%block(b)%vars(1:2,:,:,:)
      if (error_source == 0) then
        IOfield%block(b)%vars(3,:,:,:) = IOsource%block(b)%vars(3,:,:,:)
      else
        IOfield%block(b)%vars(3,:,:,:) = 0.0d0
      endif
    enddo

  end subroutine Read_vtk_tec


  !> Output setup
  subroutine Setup_Output_Solution ( IOfield )
    use IR_Precision
    use FUSS_Config_Types_m, only: obj_multigrid, obj_io
    use FUSS_Global_m
    use FUSS_Parameters_m
    use Lib_ORION_data
    use strings, only: parse
    implicit none
    type(ORION_data), intent(inout) :: IOfield(obj_multigrid%MGL)
    ! Local
    integer :: b, i, m

    ! IO Variables specification
    obj_io%Ovarnames=' "T" "matID" "qvol" "k" "rho" "cs" "h" '
    obj_io%Onvar = 1 + 6

    do m = 1, obj_multigrid%MGL
      IOfield(m)%vtk%node = .false.
      IOfield(m)%tec%node = .false.
      IOfield(m)%tec%bc = .false.
      if (index(obj_io%sol_format,'ascii')>0) then
        IOfield(m)%tec%extension = '.tec'
        IOfield(m)%tec%format = 'ascii'
        IOfield(m)%vtk%format = 'ascii'
      else
        IOfield(m)%tec%extension = '.szplt'
        IOfield(m)%tec%format = 'binary'
        IOfield(m)%vtk%format = 'binary'
      endif
    end do

    ! Concretize the sol subroutine
    Write_solution => Write_vtk_tec

    do m = 1, obj_multigrid%MGL
      do b = 1, Size(IOfield(m)%block)
        IOfield(m)%block(b)%name = 'Block'//trim(str(.true.,b))
      enddo
    enddo

    ! Reallocate IOfield vars if the number of variables read into the backup file is different from the solution one
    ! The reallocation is performed after reading the ICs and before the first solution is written!
    do m = 1, obj_multigrid%MGL
      if ( obj_io%Onvar /= Size(IOfield(m)%block(1)%vars,1) ) then
        do b = 1, Size ( IOfield(m)%block )
          IOfield(m)%block(b)%name = 'Block'//trim(str(.true.,b))
          deallocate ( IOfield(m)%block(b)%vars )
          allocate( IOfield(m)%block(b)%vars(1:obj_io%Onvar, &
                                             1:IOfield(m)%block(b)%Ni, &
                                             1:IOfield(m)%block(b)%Nj, &
                                             1:IOfield(m)%block(b)%Nk) )
        enddo
      endif
    enddo

  end subroutine Setup_Output_Solution

  subroutine Write_vtk_tec ( domain, IOfield, file )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_io
    use FUSS_Global_m
    use FUSS_Parameters_m
    use IR_Precision
    use Lib_ORION_data
    use Lib_VTK
    use Lib_Tecplot
    use FUSS_Lib_Solid
    use FUSS_Mod_MPI, only: mpi_is_root
    use FUSS_Mod_GhostExchange, only: gather_T_to_root, gather_qvol_to_root, gather_matID_to_root, mpi_io_barrier
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    type(ORION_data), intent(inout)       :: IOfield
    character(llen), intent(in)           :: file
    ! Local
    character(len=llen) :: path
    character(len=llen) :: localpath_vtk
    integer             :: E_IO, b, i, j, k
    real(8)             :: kappa, rho, cp, h

    ! Gather P from all ranks to root (collective — all ranks must call)
    call gather_T_to_root(domain)
    call gather_matID_to_root(domain)
    call gather_qvol_to_root(domain)

    if (mpi_is_root) then
      path = 'OUTPUT/'

      ! Update IOfield variables with domain primitives and other variables
      do b = 1, size(IOfield%block)
        IOfield%block(b)%vars(1,:,:,:) = domain%blk(b)%T(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        IOfield%block(b)%vars(2,:,:,:) = domain%blk(b)%matID(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        IOfield%block(b)%vars(3,:,:,:) = domain%blk(b)%qvol(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        ! Auxiliary variables
        do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
          call co_Kappa( domain%blk(b)%matID(i,j,k), domain%blk(b)%T(i,j,k), kappa )
          call co_Rho( domain%blk(b)%matID(i,j,k), domain%blk(b)%T(i,j,k), rho )
          call co_Cp( domain%blk(b)%matID(i,j,k), domain%blk(b)%T(i,j,k), cp )
          call co_H( domain%blk(b)%matID(i,j,k), domain%blk(b)%T(i,j,k), h )
          IOfield%block(b)%vars(4,i,j,k) = kappa
          IOfield%block(b)%vars(5,i,j,k) = rho
          IOfield%block(b)%vars(6,i,j,k) = cp
          IOfield%block(b)%vars(7,i,j,k) = h
        enddo; enddo; enddo
      enddo
    
      ! Write the IOfield accordingly to the solution format
      if (index(obj_io%sol_format,'vtk')>0) then
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        E_IO = vtk_write_structured_multiblock(orion=IOfield,vtspath=trim(localpath_vtk)//trim(file), &
                                                             vtmpath=trim(path)//trim(file),varnames=obj_io%Ovarnames,time=domain%time)
      else
        E_IO = tec_write_structured_multiblock(orion=IOfield,varnames=obj_io%Ovarnames,filename=trim(path)//trim(file)//trim(IOfield%tec%extension))
      end if
    end if

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine Write_vtk_tec

end module FUSS_IO_Solution