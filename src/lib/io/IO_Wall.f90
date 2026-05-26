module FUSS_IO_Wall
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Config_Types_m, only: obj_io, obj_io_bc
  use FUSS_Global_m
  use FUSS_Parameters_m
  use Lib_ORION_data
  implicit none
  type(ORION_data) :: IOwall

contains

  subroutine Initialize_Wall_File ( domain, extension )
    use FUSS_Advanced_Types_m
    use IR_Precision
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    character(len=*), intent(in)       :: extension
    ! Local
    integer :: b, f, bb, i, j, k

    IOwall%tec%extension = extension
    obj_io%Onwall = 2
    obj_io%Owallnames = '"Tw" "qw"'

    allocate(IOwall%block(1:count(obj_io_bc%viscous_flag)))

    bb = 0
    do b = 1, domain % nb; do f = 1, 6
      if (obj_io_bc%viscous_flag( b , f )) then
        bb = bb + 1
        IOwall%block(bb)%name = 'B'//trim(str(.true.,b))//'F'//trim(str(.true.,f))
        select case(f)
        case(1)
          allocate(IOwall%block(bb)%mesh(1:3,0:0,0:domain%blk(b)%dim(2),0:domain%blk(b)%dim(3)))
          allocate(IOwall%block(bb)%vars(1:obj_io%Onwall,1,1:domain%blk(b)%dim(2),1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2)
              IOwall % block(bb) % mesh(1,0,j,k) = domain % blk(b) % node(0,j,k) % c(1)
              IOwall % block(bb) % mesh(2,0,j,k) = domain % blk(b) % node(0,j,k) % c(2)
              IOwall % block(bb) % mesh(3,0,j,k) = domain % blk(b) % node(0,j,k) % c(3)
              ! Section to fool compiler...
              if (IOwall%block(bb)%mesh(1,0,j,k)==0d0 .and. &
                  IOwall%block(bb)%mesh(2,0,j,k)==0d0 .and. &
                  IOwall%block(bb)%mesh(3,0,j,k)==0d0) &
              print*, IOwall % block(bb) % mesh(:,0,j,k)
          enddo; enddo
        case(2)
          allocate(IOwall%block(bb)%mesh(1:3,0:0,0:domain%blk(b)%dim(2),0:domain%blk(b)%dim(3)))
          allocate(IOwall%block(bb)%vars(1:obj_io%Onwall,1,1:domain%blk(b)%dim(2),1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2)
              IOwall % block(bb) % mesh(:,0,j,k) = domain % blk(b) % node(domain%blk(b)%dim(1),j,k) % c
          enddo; enddo
        case(3)
          allocate(IOwall%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:0,0:domain%blk(b)%dim(3)))
          allocate(IOwall%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1,1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do i = 0, domain%blk(b)%dim(1)
              IOwall % block(bb) % mesh(:,i,0,k) = domain % blk(b) % node(i,0,k) % c
          enddo; enddo
        case(4)
          allocate(IOwall%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:0,0:domain%blk(b)%dim(3)))
          allocate(IOwall%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1,1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do i = 0, domain%blk(b)%dim(1)
              IOwall % block(bb) % mesh(:,i,0,k) = domain % blk(b) % node(i,domain%blk(b)%dim(2),k) % c
          enddo; enddo
        case(5)
          allocate(IOwall%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:domain%blk(b)%dim(2),0:0))
          allocate(IOwall%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1:domain%blk(b)%dim(2),1))
          do j = 0, domain%blk(b)%dim(2); do i = 0, domain%blk(b)%dim(1)
              IOwall % block(bb) % mesh(:,i,j,0) = domain % blk(b) % node(i,j,0) % c
          enddo; enddo
        case(6)
          allocate(IOwall%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:domain%blk(b)%dim(2),0:0))
          allocate(IOwall%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1:domain%blk(b)%dim(2),1))
          do j = 0, domain%blk(b)%dim(2); do i = 0, domain%blk(b)%dim(1)
              IOwall % block(bb) % mesh(:,i,j,0) = domain % blk(b) % node(i,j,domain%blk(b)%dim(3)) % c
          enddo; enddo
        end select
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
    integer             :: E_IO
    character(len=clen) :: format(2)

    ! Each rank computes wall properties for its local blocks and reduces to root.
    ! (collective — all ranks must call)
    call Extract_Properties_From_BC( domain )

    if (mpi_is_root) then
      path = 'OUTPUT/'
      call parse(obj_io%sol_format,' ', format)

      ! Write the IOwall accordingly to the solution format
      select case(trim(format(1)))
      case('vtk')
        IOwall%vtk%format = trim(format(2))
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        E_IO = vtk_write_structured_multiblock(orion=IOwall,vtspath=trim(localpath_vtk)//trim(file), &
                                                             vtmpath=trim(path)//trim(file),varnames=obj_io%Owallnames,time=domain%time)
      case('tecplot')
        IOwall%tec%format = trim(format(2))
        E_IO = tec_write_structured_multiblock(orion=IOwall,varnames=obj_io%Owallnames,filename=trim(path)//trim(file)//trim(IOwall%tec%extension))
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
    use FUSS_Lib_BC_Fluxes_Wall_Heat
    use FUSS_Lib_BC_Fluxes_Wall_Temperature
    use FUSS_Lib_BC_Fluxes_Wall_HeatTransfer
    use FUSS_Lib_BC_Fluxes_Wall_Radiation
    use FUSS_Lib_BC_Fluxes_Wall_RadiationConvection
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    real(R8), dimension(2) :: Ovar
    integer :: f, lower, upper, i, b, bb, Bm, Im, Jm, Km, Fm, Is, Js, Ks
    integer :: n_wall_total, wall_i
    real(R8), allocatable :: wall_buf(:)

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
          select case ( domain % bc(i) % type )

            case (301) ! wall: prescribed heat flux
              call BC_Wall_Heat ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                  domain % bc(i) % qw, Ovar )

            case (302) ! wall: prescribed temperature
              call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), & 
                                         domain % bc(i) % Tw, Ovar )

            case (303) ! wall: prescribed hconv, Tref, qrad
              call BC_Wall_HeatTransfer ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                          domain % bc(i) % hconv, domain % bc(i) % Tref, domain % bc(i) % qw, Ovar )

            case (304) ! wall: prescribed epsilon, Tref
              call BC_Wall_Radiation ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                       domain % bc(i) % eps, domain % bc(i) % Tref, Ovar )
              
            case (305) ! wall: prescribed hconv, epsilon, Tref
              call BC_Wall_RadiationConvection ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                                 domain % bc(i) % hconv, domain % bc(i) % eps, domain % bc(i) % Tref, Ovar )

          end select
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
            IOwall%block(bb)%vars(1, Is, Js, Ks) = Ovar(1)
            IOwall%block(bb)%vars(2, Is, Js, Ks) = Ovar(2)
          enddo

        enddo ! face
      enddo ! block
    end if

    deallocate(wall_buf)

  end subroutine Extract_Properties_From_BC
  
end module FUSS_IO_Wall