module FUSS_Mod_Allocate_Data
  use, intrinsic :: iso_fortran_env, only : iostat_end
  
  implicit none
  private
  public :: Setup_Data_Structure, Allocate_Block, deallocate_remote_computation_data

contains

  subroutine Setup_Data_Structure ( domain, IOfield )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_io
    use FUSS_Global_m
    use Lib_ORION_data
    
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    type(ORION_data), intent(inout)       :: IOfield
    ! Local
    integer :: b, d, nblocks
   
    ! Domain is the FUSS-alias of IOfield
    nblocks = size ( IOfield%block )
    allocate( domain%blk( 1:nblocks ) )
    domain%nb = nblocks

    ! Define the ijk dimensions of each block
    do b = 1, nblocks
      domain%blk(b)%dim(1) = IOfield%block(b)%Ni
      domain%blk(b)%dim(2) = IOfield%block(b)%Nj
      domain%blk(b)%dim(3) = IOfield%block(b)%Nk
    enddo
    
    ! Check if number of simulation variables in orion-field matches FUSS expectation
    if ( size(IOfield%block(1)%vars, 1) < 3 ) then
      write(*,'(A)')         '[ERROR] Number of variables in IOfield does not match FUSS expectation.'
      write(*,'(A,I0,A,I0)') '        Expected: ', 3, ', Found: ', size(IOfield%block(1)%vars, 1)
      stop
    end if

    do b = 1, nblocks

      ! Allocate domain block
      call Allocate_Block ( domain%blk(b), domain%blk(b)%dim )
      
      ! Import domain-block nodes from orion-field
      do d = 1, 3
        domain%blk(b)%node(0:IOfield%block(b)%Ni,0:IOfield%block(b)%Nj,0:IOfield%block(b)%Nk)%c(d) &
        = IOfield%block(b)%mesh(d,0:,0:,0:)
      enddo
      ! Import domain-block primitives from orion-field variables
      domain%blk(b)%T(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk) &
      = IOfield%block(b)%vars(1, 1:,1:,1:)
      domain%blk(b)%matID(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk) &
      = IOfield%block(b)%vars(2, 1:,1:,1:)
      domain%blk(b)%qvol(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk) &
      = IOfield%block(b)%vars(3, 1:,1:,1:)
      domain%time = obj_io%IOtime

    enddo

  end subroutine Setup_Data_Structure


  subroutine Allocate_Block( blk, nijk )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_irs
    use FUSS_Global_m
    implicit none
    integer, intent(in)                  :: nijk(3)
    type(FUSS_block_type), intent(inout) :: blk
    ! Local
    integer :: ni, nj, nk

    ni = nijk(1) ; nj = nijk(2) ; nk = nijk(3)

    ! Metrics
    allocate( blk % node ( 0:ni, 0:nj, 0:nk ) )
    allocate( blk % M    ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % dl   ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % vol  ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % dir(1) % f (0:ni, nj, nk) )
    allocate( blk % dir(2) % f (ni, 0:nj, nk) )
    allocate( blk % dir(3) % f (ni, nj, 0:nk) )

    ! Prim and Residuals
    allocate( blk % T ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % TO, blk % R,  blk % matID, blk % qvol, mold = blk % T )

    ! Dt cell center with no ghost cells
    allocate( blk % dtlocal ( 1:ni, 1:nj, 1:nk) )

    ! Temp storage for residuals in IRS
    if ( obj_irs%enabled ) then
      allocate( blk % RS1, blk % RS2, mold = blk % R )
    end if

  end subroutine Allocate_Block
  

  subroutine deallocate_remote_computation_data(domain)
    use FUSS_Advanced_Types_m
    use FUSS_Mod_MPI, only: is_local_block, mpi_is_root

    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    integer :: b, d, i, c
    logical, allocatable :: needs_remote_T(:)

    ! Build mask of remote blocks whose T (and dir) must be kept:
    !  - chimera (102): donorID(:,1) can reference remote blocks
    allocate(needs_remote_T(domain%nb))
    needs_remote_T = .false.
    do i = 1, domain%nbound
      select case (domain%bc(i)%type)
        case (102) ! chimera
          if (allocated(domain%bc(i)%donorID)) then
            do c = 1, size(domain%bc(i)%donorID, 1)
              b = domain%bc(i)%donorID(c, 1)
              if (.not. is_local_block(b)) needs_remote_T(b) = .true.
            end do
          end if
      end select
    end do

    do b = 1, domain%nb
      if (is_local_block(b)) cycle

      ! Computation arrays — free on all ranks
      if (allocated(domain%blk(b)%TO))           deallocate(domain%blk(b)%TO)
      if (allocated(domain%blk(b)%R))            deallocate(domain%blk(b)%R)
      if (allocated(domain%blk(b)%RS1))          deallocate(domain%blk(b)%RS1)
      if (allocated(domain%blk(b)%RS2))          deallocate(domain%blk(b)%RS2)
      if (allocated(domain%blk(b)%dtlocal))      deallocate(domain%blk(b)%dtlocal)
      if (allocated(domain%blk(b)%qvol))         deallocate(domain%blk(b)%qvol)
      if (allocated(domain%blk(b)%matID))        deallocate(domain%blk(b)%matID)

      ! T — free unless this remote block is a chimera donor or manifold source
      if (.not. needs_remote_T(b)) then
        if (allocated(domain%blk(b)%T)) deallocate(domain%blk(b)%T)
      end if

      ! Metrics — free on non-root ranks only (root needs them for wall I/O)
      ! Keep dir on blocks needed for manifold (BC_Manifold reads blk(Bs)%dir%f%A)
      if (.not. mpi_is_root) then
        if (allocated(domain%blk(b)%node)) deallocate(domain%blk(b)%node)
        if (allocated(domain%blk(b)%M))    deallocate(domain%blk(b)%M)
        if (allocated(domain%blk(b)%dl))   deallocate(domain%blk(b)%dl)
        if (allocated(domain%blk(b)%vol))  deallocate(domain%blk(b)%vol)
        if (.not. needs_remote_T(b)) then
          do d = 1, 3
            if (allocated(domain%blk(b)%dir(d)%f)) deallocate(domain%blk(b)%dir(d)%f)
          end do
        end if
      end if
    end do

    deallocate(needs_remote_T)

  end subroutine deallocate_remote_computation_data


end module FUSS_Mod_Allocate_Data