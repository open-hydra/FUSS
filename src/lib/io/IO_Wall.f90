module FUSS_IO_Wall
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Config_Types_m, only: obj_io, obj_io_bc
  use FUSS_Global_m
  use FUSS_Parameters_m
  use Lib_ORION_data
  implicit none
  ! Wall output, one container per multigrid level (index = domain%mg_level): allocated on
  ! every rank, filled and written on root only
  type(ORION_data), allocatable :: IOwall(:)

contains

  subroutine Initialize_Wall_File ( domain, extension )
    use FUSS_Advanced_Types_m
    use IR_Precision
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    character(len=*), intent(in)       :: extension
    ! Local
    integer :: b, f, bb, i, j, k, lv

    lv = domain%mg_level
    IOwall(lv)%tec%extension = extension
    obj_io%Onwall = 2
    obj_io%Owallnames = '"Tw" "qw"'

    allocate(IOwall(lv)%block(1:count(obj_io_bc%viscous_flag)))

    bb = 0
    do b = 1, domain % nb; do f = 1, 6
      if (obj_io_bc%viscous_flag( b , f )) then
        bb = bb + 1
        IOwall(lv)%block(bb)%name = 'B'//trim(str(.true.,b))//'F'//trim(str(.true.,f))
        select case(f)
        case(1)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:0,0:domain%blk(b)%dim(2),0:domain%blk(b)%dim(3)))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1,1:domain%blk(b)%dim(2),1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2)
              IOwall(lv) % block(bb) % mesh(1,0,j,k) = domain % blk(b) % node(0,j,k) % c(1)
              IOwall(lv) % block(bb) % mesh(2,0,j,k) = domain % blk(b) % node(0,j,k) % c(2)
              IOwall(lv) % block(bb) % mesh(3,0,j,k) = domain % blk(b) % node(0,j,k) % c(3)
              ! Section to fool compiler...
              if (IOwall(lv)%block(bb)%mesh(1,0,j,k)==0d0 .and. &
                  IOwall(lv)%block(bb)%mesh(2,0,j,k)==0d0 .and. &
                  IOwall(lv)%block(bb)%mesh(3,0,j,k)==0d0) &
              print*, IOwall(lv) % block(bb) % mesh(:,0,j,k)
          enddo; enddo
        case(2)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:0,0:domain%blk(b)%dim(2),0:domain%blk(b)%dim(3)))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1,1:domain%blk(b)%dim(2),1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2)
              IOwall(lv) % block(bb) % mesh(:,0,j,k) = domain % blk(b) % node(domain%blk(b)%dim(1),j,k) % c
          enddo; enddo
        case(3)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:0,0:domain%blk(b)%dim(3)))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1,1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do i = 0, domain%blk(b)%dim(1)
              IOwall(lv) % block(bb) % mesh(:,i,0,k) = domain % blk(b) % node(i,0,k) % c
          enddo; enddo
        case(4)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:0,0:domain%blk(b)%dim(3)))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1,1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do i = 0, domain%blk(b)%dim(1)
              IOwall(lv) % block(bb) % mesh(:,i,0,k) = domain % blk(b) % node(i,domain%blk(b)%dim(2),k) % c
          enddo; enddo
        case(5)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:domain%blk(b)%dim(2),0:0))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1:domain%blk(b)%dim(2),1))
          do j = 0, domain%blk(b)%dim(2); do i = 0, domain%blk(b)%dim(1)
              IOwall(lv) % block(bb) % mesh(:,i,j,0) = domain % blk(b) % node(i,j,0) % c
          enddo; enddo
        case(6)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:domain%blk(b)%dim(2),0:0))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1:domain%blk(b)%dim(2),1))
          do j = 0, domain%blk(b)%dim(2); do i = 0, domain%blk(b)%dim(1)
              IOwall(lv) % block(bb) % mesh(:,i,j,0) = domain % blk(b) % node(i,j,domain%blk(b)%dim(3)) % c
          enddo; enddo
        end select
        ! Cells of the face in each direction, 0 across it: the vtk writer of ORION (a978930) sizes the block
        ! file with them (the Tecplot writer takes the sizes from the mesh array); as in MOSE (IO_Wall.f90:100-104)
        IOwall(lv)%block(bb)%Ni = size(IOwall(lv)%block(bb)%mesh, 2) - 1
        IOwall(lv)%block(bb)%Nj = size(IOwall(lv)%block(bb)%mesh, 3) - 1
        IOwall(lv)%block(bb)%Nk = size(IOwall(lv)%block(bb)%mesh, 4) - 1
      endif
    enddo; enddo

  end subroutine Initialize_Wall_File


  !> Update the orion-field data and write them accordingly to the chosen format (vtk,tecplot)
  subroutine Write_Wall_Solution( domain, file )
    use Lib_VTK
    use Lib_Tecplot
    use FUSS_Advanced_Types_m
    use FUSS_Parameters_m
    use FUSS_Mod_MPI, only: mpi_is_root
    use FUSS_Mod_GhostExchange, only: mpi_io_barrier
    use strings, only: parse
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    character(len=llen), intent(in)       :: file
    ! Local
    character(len=llen) :: path
    character(len=llen) :: localpath_vtk
    character(len=llen) :: vtknames
    character(len=llen) :: name
    integer             :: E_IO, lv
    character(len=clen) :: format(2)

    lv = domain%mg_level

    ! Each rank computes wall properties for its local blocks and reduces to root.
    ! (collective — all ranks must call)
    call Extract_Properties_From_BC( domain )

    if (mpi_is_root) then
      path = 'OUTPUT/'
      call parse(obj_io%sol_format,' ', format)

      ! Write the IOwall accordingly to the solution format
      select case(trim(format(1)))
      case('vtk')
        IOwall(lv)%vtk%format = trim(format(2))
        ! The .vtm lists the block files relative to its own folder: the leading '/' of file is dropped, as for the field
        ! files (IO_Solution); the time is the solution time of the Tecplot outputs, minus the global iteration count in a
        ! steady run (it was the -1 of domain%time)
        name = file
        if ( name(1:1) == '/' ) name = name(2:)
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        vtknames = unquoted(obj_io%Owallnames)
        E_IO = vtk_write_structured_multiblock(orion=IOwall(lv),vtspath=trim(localpath_vtk)//trim(name), &
                                               vtmpath=trim(path)//trim(name),varnames=vtknames,time=IOwall(lv)%solutiontime)
      case('tecplot')
        IOwall(lv)%tec%format = trim(format(2))
        E_IO = tec_write_structured_multiblock(orion=IOwall(lv),varnames=obj_io%Owallnames,filename=trim(path)//trim(file)//trim(IOwall(lv)%tec%extension))
      end select
    end if

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine Write_Wall_Solution


  !> Compute wall BC properties on all ranks (local blocks only), reduce to root,
  !> and fill IOwall on root. Collective — all ranks must call.
  subroutine Extract_Properties_From_BC ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Global_m
    use FUSS_Mod_MPI, only: mpi_is_root, is_local_block, mpi_reduce_sum_r8_array
    use FUSS_Lib_BC_Fluxes_Wall, only: Wall_Temperature
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    real(R8), dimension(2) :: Ovar
    integer :: f, lower, upper, i, b, bb, Bm, Im, Jm, Km, Fm, Is, Js, Ks
    integer :: n_wall_total, wall_i, lv
    real(R8), allocatable :: wall_buf(:)

    lv = domain%mg_level

    ! Count total BC cells across all viscous faces
    n_wall_total = 0
    do b = 1, domain%nb
      do f = 1, 6
        if (obj_io_bc%viscous_flag(b,f)) n_wall_total = n_wall_total + domain%n_bf(b,f)
      end do
    end do

    allocate(wall_buf(2 * n_wall_total))
    wall_buf = 0d0

    ! Each rank computes Ovar only for its local blocks
    upper = 0; wall_i = 0
    do b = 1, domain%nb
      do f = 1, 6
        lower = upper + 1
        upper = upper + domain%n_bf(b,f)
        if (.not. obj_io_bc%viscous_flag(b,f)) cycle

        do i = lower, upper
          wall_i = wall_i + 1
          Bm = domain%bc(i)%b
          if (.not. is_local_block(Bm)) cycle

          Im = domain % bc(i) % i
          Jm = domain % bc(i) % j
          Km = domain % bc(i) % k
          Fm = domain % bc(i) % f
          Ovar = 0d0
          if (domain % bc(i) % type >= 301 .and. domain % bc(i) % type <= 305) then
            ! Time tables at the output time (the ghost fill leaves the last stage time)
            if (domain % bc(i) % BCtime % exists) then
              if (domain % bc(i) % type == 301) domain % bc(i) % qw = domain % bc(i) % BCtime % update(domain % time)
              if (domain % bc(i) % type == 302) domain % bc(i) % Tw = domain % bc(i) % BCtime % update(domain % time)
            endif
            ! ghosts=.false.: interior values only, current at any point of the step (the ghosts
            ! hold the last stage); one-sided tangential differences at the face edges
            call Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i), .false., Ovar(1), Ovar(2) )
          endif
          wall_buf(2*(wall_i-1)+1 : 2*wall_i) = Ovar
        enddo

      enddo ! face
    enddo ! block

    ! Reduce wall_buf to root. Each cell is owned by exactly one rank, so SUM gives
    ! correct values.
    call mpi_reduce_sum_r8_array(wall_buf, 2*n_wall_total)

    ! Root fills IOwall from the gathered buffer
    if (mpi_is_root) then
      upper = 0; wall_i = 0; bb = 0
      do b = 1, domain%nb
        do f = 1, 6
          lower = upper + 1
          upper = upper + domain%n_bf(b,f)
          if (.not. obj_io_bc%viscous_flag(b,f)) cycle
          bb = bb + 1

          do i = lower, upper
            wall_i = wall_i + 1
            Ovar = wall_buf(2*(wall_i-1)+1 : 2*wall_i)

            select case(f)
              case(1,2); Is = 1               ; Js = domain%bc(i)%j ; Ks = domain%bc(i)%k
              case(3,4); Is = domain%bc(i)% i ; Js = 1              ; Ks = domain%bc(i)%k
              case(5,6); Is = domain%bc(i)% i ; Js = domain%bc(i)%j ; Ks = 1
            end select
            IOwall(lv)%block(bb)%vars(1, Is, Js, Ks) = Ovar(1)
            IOwall(lv)%block(bb)%vars(2, Is, Js, Ks) = Ovar(2)
          enddo

        enddo ! face
      enddo ! block
    end if

    deallocate(wall_buf)

  end subroutine Extract_Properties_From_BC
  
end module FUSS_IO_Wall