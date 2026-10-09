module FUSS_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_BC
  public :: Print_BC_Summary

  ! BC type counters — populated by Read_BCfile, consumed by Print_BC_Summary
  integer :: nconnect, nwall, nsym, nchimera, ncoupled
  logical :: has_tdep_bc

contains

  subroutine Setup_BC ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_multigrid, obj_io_bc, obj_pyrolysis
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    ! Local
    integer :: m, error, i

    do m = 1, obj_multigrid%MGL

      !! Phase 1: Allocate and check
      call Allocate_BC ( domain(m) )
      
      if (allocated(obj_io_bc%viscous_flag)) deallocate(obj_io_bc%viscous_flag)
      allocate ( obj_io_bc%viscous_flag( 6 * domain(m) % nb, 6 ) )
      obj_io_bc%viscous_flag = .false.

      if (allocated(obj_io_bc%coupling_flag)) deallocate(obj_io_bc%coupling_flag)
      allocate ( obj_io_bc%coupling_flag( 6 * domain(m) % nb, 6 ) )    
      obj_io_bc%coupling_flag = .false.
      
      error = Check_BC ( domain(m) % nbound, m )
      if (error /= 0) cycle

      !! Phase 2: Read
      call Read_BCfile ( domain(m) % bc, domain(m) % n_bf, m )

      ! Pyrolysis (plan 06 section 3.1): the ghost chi stencil of every BC
      ! entry is allocated and zeroed here; plan 08 fills it at connections.
      if ( obj_pyrolysis%enabled ) then
        do i = 1, domain(m) % nbound
          if ( allocated(domain(m) % bc(i) % chig) ) deallocate( domain(m) % bc(i) % chig )
          allocate( domain(m) % bc(i) % chig( obj_pyrolysis%nsub ) )
          domain(m) % bc(i) % chig = 0.0d0
        enddo
      endif
      
    end do

  end subroutine Setup_BC
  

  subroutine Allocate_BC ( domain )
    use FUSS_Advanced_Types_m
    implicit none
    type(FUSS_domain_type) :: domain
    ! Local
    integer :: b, n, ni, nj, nk

    n = 0
    do b = 1, domain % nb
      ni = domain % blk(b) % dim(1)
      nj = domain % blk(b) % dim(2)
      nk = domain % blk(b) % dim(3)
      n = n + 2*nj*nk + 2*ni*nk + 2*nj*ni
    enddo

    domain % nbound = n
    allocate ( domain % bc ( domain % nbound ) )
    allocate ( domain % n_bf ( domain % nb, 6 ) )
    
  end subroutine Allocate_BC


  function Check_BC (n, level) result(ios)
    use FUSS_Config_Types_m, only: obj_io_bc
    use FUSS_Global_m,       only: FUSS_phase_prefix
    use IR_Precision,        only: str
    implicit none
    integer, intent(in)  :: n, level
    integer              :: ios
    integer              :: di(5), ti, unitfile, n_proof
    integer              :: c, ci, cii

    ios = 0

    ! Open file
    if (level == 1) then
      open(newunit=unitfile,file='INPUT/'//trim(FUSS_phase_prefix)//'bc.txt',status='old',iostat=ios,action='read')
    else
      open(newunit=unitfile,file='INPUT/'//trim(FUSS_phase_prefix)//'bc'//trim(str(.true.,level))//'.txt',status='old',iostat=ios,action='read')
    endif
    if (ios/=0) then
      obj_io_bc%error_message = '[ERROR] Boundary condition file not found for grid '//trim(str(.true.,level))
      return
    endif

    ! Cheak BC file consistency
    ios = 0; n_proof = -1
    do while (ios==0)
      read( unitfile,*,iostat=ios ) di(1), di(2), di(3), di(4), di(5), ti
      select case(ti)
      case(101,301,302,303,304,305,103)
        read( unitfile,*,iostat=ios )
      case(102)
        read( unitfile,*,iostat=ios ) ci, cii
        do c = 1, ci+cii
          read( unitfile,*,iostat=ios )
        enddo
      end select
      n_proof = n_proof + 1
    enddo

    if (n_proof /= n) then
      obj_io_bc%error_message = '[ERROR] Boundary conditions number ('//str(.true.,n_proof)//') is different than the one of the initial conditions ('//str(.true.,n)//')'
      close(unitfile)
      return
    endif

    ! Validation passed: reset ios (non-zero from EOF) to signal success
    ios = 0
    close(unitfile)

  end function Check_BC


  subroutine Read_BCfile ( bc, n_bf, level )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_io, obj_io_bc, obj_darcy
    use FUSS_Global_m
    use IR_Precision
    implicit none
    type(FUSS_bc_type), dimension(:), intent(inout) :: bc
    integer, intent(in)                             :: level
    integer, dimension(1:,1:), intent(inout)        :: n_bf
    ! Local
    integer :: cc, i, s
    integer :: unitfile, ios, cios
    character(len=32) :: BCfile
    character(len=256) :: rec

    cios = 0

    ! Open file
    if (level == 1) then
      open(newunit=unitfile,file='INPUT/'//trim(FUSS_phase_prefix)//'bc.txt',status='old',iostat=ios,action='read')
    else
      open(newunit=unitfile,file='INPUT/'//trim(FUSS_phase_prefix)//'bc'//trim(str(.true.,level))//'.txt',status='old',iostat=ios,action='read')
    endif
    if (ios/=0) then
      obj_io_bc%error_message = '[ERROR] Boundary condition file not found for grid '//trim(str(.true.,level))
      return
    endif

    ! Counters for specific BC types
    if (level == 1) then
      nconnect = 0
      nwall = 0
      nsym = 0
      nchimera = 0
      ncoupled = 0
      has_tdep_bc = .false.
    endif

    ! Counter for number of cells per face in each block
    n_bf = 0

    ! Read file
    do i = 1, size(bc)

      ! First line is equal for every BC type
      read( unitfile,*,iostat=ios ) bc(i)%b, bc(i)%i, bc(i)%j, bc(i)%k, bc(i)%f, bc(i)%type
      if (ios/=0) write(*,'(A)') '  Error in BC file'

      ! n_bf update
      n_bf( bc(i) % b, bc(i) % f ) = n_bf( bc(i) % b, bc(i) % f ) + 1

      ! Open-face gas pressure of plan 07 (read only under type 2): p-wall
      ! unless a 303/305 value line carries its own fourth real below.
      bc(i)%pw = obj_darcy%p_wall

      ! Second line is BC type-dependent
      select case( bc(i)%type )
        case(101)
          if (level == 1) nconnect = nconnect + 1
          read( unitfile,*,iostat=ios ) &
          bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22

        case(300)
          if (level == 1) nsym = nsym + 1
        
        case(301)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,*,iostat=ios) BCfile
          ! Time-varying properties: if the read string is not convertible to a real, the BC is time-dependent
          read(BCfile,*,iostat=cios) bc(i)%qw
          if (cios/=0) then
            call bc(i)%BCtime%initialize(file=BCfile,bar=.true.)
            has_tdep_bc = .true.
          endif
        
        case(302)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,*,iostat=ios) BCfile
          ! Time-varying properties: if the read string is not convertible to a real, the BC is time-dependent
          read(BCfile,*,iostat=cios) bc(i)%Tw
          if (cios/=0) then
            call bc(i)%BCtime%initialize(file=BCfile,bar=.true.)
            has_tdep_bc = .true.
          endif

        case(303)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          ! `hconv, qrad, Tref[, p_w]` (plan 07 section 4.3): the record is
          ! read into a buffer and parsed for four values first; a three-value
          ! line fails that read (end of the internal record) and is parsed
          ! again for three, exactly as before, so no existing case moves.
          ! A null fourth item (a trailing comma) leaves the -huge sentinel (a
          ! sentinel a user could type, such as -1, would be read as "not given").
          read(unitfile,'(A)',iostat=ios) rec
          bc(i)%pw = -huge(1.0d0)
          read(rec,*,iostat=cios) bc(i)%hconv, bc(i)%qrad, bc(i)%Tref, bc(i)%pw
          if ( cios /= 0 .or. bc(i)%pw == -huge(1.0d0) ) then
            read(rec,*,iostat=ios) bc(i)%hconv, bc(i)%qrad, bc(i)%Tref
            bc(i)%pw = obj_darcy%p_wall
          elseif ( bc(i)%pw <= 0.0d0 ) then
            write(*,'(A,I0,A)') ' [ERROR] bc.txt record ', i, ': the fourth value of a 303 line (p_w) must be > 0'
            error stop 1
          endif

        case(304)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,*,iostat=ios)  bc(i)%eps, bc(i)%Tref

        case(305)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,'(A)',iostat=ios) rec
          bc(i)%pw = -huge(1.0d0)
          read(rec,*,iostat=cios) bc(i)%hconv, bc(i)%eps, bc(i)%Tref, bc(i)%pw
          if ( cios /= 0 .or. bc(i)%pw == -huge(1.0d0) ) then
            read(rec,*,iostat=ios) bc(i)%hconv, bc(i)%eps, bc(i)%Tref
            bc(i)%pw = obj_darcy%p_wall
          elseif ( bc(i)%pw <= 0.0d0 ) then
            write(*,'(A,I0,A)') ' [ERROR] bc.txt record ', i, ': the fourth value of a 305 line (p_w) must be > 0'
            error stop 1
          endif

        case(102)
          if (level == 1) nchimera = nchimera + 1
          read( unitfile,*,iostat=ios ) (bc(i)%ni(cc),cc=1,2)
          allocate(bc(i)%donorID(1:sum(bc(i)%ni),1:4))
          allocate(bc(i)%volume_fraction(1:sum(bc(i)%ni)))
          do s = 1, bc(i)%ni(1)
            read( unitfile,*,iostat=ios ) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
          enddo
          do s = bc(i)%ni(1)+1, bc(i)%ni(1)+bc(i)%ni(2)
            read( unitfile,*,iostat=ios ) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
          enddo

        case(103)
          if (level == 1) ncoupled = ncoupled + 1
          obj_io_bc%coupling_flag( bc(i)%b , bc(i)%f ) = .true.
          read( unitfile,*,iostat=ios ) &
          bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          bc(i)%ext_flux = 0.0
      
      end select
      
    enddo

    if (nwall > 0) obj_io % write_wall = .true.

    close( unitfile )

  end subroutine Read_BCfile


  subroutine Print_BC_Summary ()
    implicit none

    write(*,*)
    write(*,'(A)') ' Boundary conditions'
    if (nconnect > 0) write(*,'(A,T35,I0)') '   Connection', nconnect
    if (nwall > 0)    write(*,'(A,T35,I0)') '   Viscous wall', nwall
    if (nsym > 0)     write(*,'(A,T35,I0)') '   Symmetry', nsym
    if (nchimera > 0) write(*,'(A,T35,I0)') '   Chimera', nchimera
    if (ncoupled > 0) write(*,'(A,T35,I0)') '   Coupled wall', ncoupled
    if (has_tdep_bc) write(*,'(A)') '   Time-dependent BC detected'

  end subroutine Print_BC_Summary

end module FUSS_IO_BC
