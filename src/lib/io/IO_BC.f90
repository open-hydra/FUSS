module FUSS_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_BC
  public :: Print_BC_Summary

  ! BC type counters — populated by Read_BCfile, consumed by Print_BC_Summary
  integer :: nconnect, nwall, nsym, nchimera, ncoupled, nablating, ncoupled_abl
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
      case(101,301,302,303,304,305,306,307,103)
        read( unitfile,*,iostat=ios )
      case(102,104)
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
    use FUSS_Config_Types_m, only: obj_io, obj_io_bc, obj_darcy, obj_ablation, obj_multigrid
    use FUSS_Global_m
    use FUSS_Lib_Bprime, only: Load_Bprime_Table
    use IR_Precision
    implicit none
    type(FUSS_bc_type), dimension(:), intent(inout) :: bc
    integer, intent(in)                             :: level
    integer, dimension(1:,1:), intent(inout)        :: n_bf
    ! Local
    integer :: cc, i, s, nzero_chim, ntok
    integer :: unitfile, ios, cios
    character(len=32) :: BCfile
    character(len=256) :: rec, bfile, efile
    character(len=1024) :: line306

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
      nablating = 0
      ncoupled_abl = 0
      has_tdep_bc = .false.
    endif

    ! Counter for number of cells per face in each block
    n_bf = 0
    nzero_chim = 0

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

        case(306)
          ! Ablating wall (plan 08 section 2.1): `lambda  T_inf  <bprime-file>
          ! <env-file>`, the files bare names under INPUT/, read through a
          ! buffer of its own (two 255-character names do not fit the shared
          ! rec buffer); `none` = no B' table, legal only when every env row
          ! has chem = 0 (checked at setup, Check_Ablation).
          if (level == 1) nwall = nwall + 1
          if (level == 1) nablating = nablating + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          ! refused here, at the first 306 record of the fine grid: a coarse
          ! grid would need its own bc file before any setup check could run
          if ( obj_multigrid%MGL > 1 ) then
            write(*,'(A)') ' [ERROR] the ablating wall (306) is not compatible with multigrid.'
            error stop 1
          endif
          ! One record, then an internal read: a list-directed read from the
          ! unit would run on into the next record when the line is short.
          bfile = '' ; efile = ''
          read(unitfile,'(A)',iostat=ios) line306
          if ( ios == 0 ) read(line306,*,iostat=ios) bc(i)%lambda, bc(i)%Tinf, bfile, efile
          if ( ios /= 0 .or. len_trim(efile) == 0 ) then
            write(*,'(A,I0,A)') ' [ERROR] bc.txt record ', i, &
              ': a 306 line is "lambda  T_inf  <bprime-file>  <env-file>" (bare names under INPUT/, or none for the table)'
            error stop 1
          endif
          if ( bc(i)%lambda < 0.0d0 .or. bc(i)%Tinf <= 0.0d0 ) then
            write(*,'(A,I0,A)') ' [ERROR] bc.txt record ', i, ': a 306 line needs lambda >= 0 and T_inf > 0'
            error stop 1
          endif
          call bc(i)%env%load( 'INPUT/'//trim(efile), 7 )
          bc(i)%ibprime = 0
          if ( trim(bfile) /= 'none' ) bc(i)%ibprime = Load_Bprime_Table( 'INPUT/'//trim(bfile) )
          if ( level == 1 ) then
            obj_ablation%enabled = .true.
            obj_ablation%n_faces = obj_ablation%n_faces + 1
          endif

        case(307)
          ! Coupled ablating wall (plan 08 section 2.5 ii): the 103
          ! connectivity record `bs is js ks fs d11 d12 d21 d22`, then `T_inf
          ! <bprime-file>` -- exactly 11 tokens; no lambda, a coupled face runs
          ! at lambda = 0 (the CFD's coefficient is already blown). The face is
          ! 306's balance with the env taken from the exchange (env_now).
          if (level == 1) nwall = nwall + 1
          if (level == 1) nablating = nablating + 1
          if (level == 1) ncoupled_abl = ncoupled_abl + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          obj_io_bc%coupling_flag( bc(i)%b , bc(i)%f ) = .true.
          if ( obj_multigrid%MGL > 1 ) then
            write(*,'(A)') ' [ERROR] the ablating wall (306, 307) is not compatible with multigrid.'
            error stop 1
          endif
          bfile = ''
          read(unitfile,'(A)',iostat=ios) line306
          ntok = 0
          if ( ios == 0 ) ntok = Count_Tokens( line306 )
          if ( ios == 0 .and. ntok == 11 ) read(line306,*,iostat=ios) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22, bc(i)%Tinf, bfile
          if ( ios /= 0 .or. ntok /= 11 ) then
            write(*,'(A,I0,A,I0,A)') ' [ERROR] bc.txt record ', i, ': a 307 line is "bs is js ks fs d11 d12 d21 d22  T_inf  '// &
              '<bprime-file>" (11 values, ', ntok, ' read): a coupled face takes no lambda, it runs at lambda = 0'
            error stop 1
          endif
          if ( bc(i)%Tinf <= 0.0d0 .or. trim(bfile) == 'none' ) then
            write(*,'(A,I0,A)') ' [ERROR] bc.txt record ', i, ': a 307 line needs T_inf > 0 and a B'' table '// &
              '(the chemistry of a coupled face is always on)'
            error stop 1
          endif
          bc(i)%lambda = 0.0d0
          bc(i)%ext_flux = 0.0
          bc(i)%ibprime = Load_Bprime_Table( 'INPUT/'//trim(bfile) )
          if ( level == 1 ) then
            obj_ablation%enabled = .true.
            obj_ablation%n_faces = obj_ablation%n_faces + 1
            obj_ablation%n_coupled = obj_ablation%n_coupled + 1
          endif

        case(102,104)
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
          if (bc(i)%type == 104) then
            ! 104: coupled wall as for 103 (coupling flux, wall output)
            obj_io_bc%coupling_flag( bc(i)%b , bc(i)%f ) = .true.
            bc(i)%ext_flux = 0.0
            ! Both ghost rows need MOSE donors, as checked on the MOSE side: without
            ! them the blended fluid state is 0 and cons2prim divides by rho = 0
            if (sum(bc(i)%volume_fraction(1:bc(i)%ni(1))) < 0.5d0 .or. &
                sum(bc(i)%volume_fraction(bc(i)%ni(1)+1:sum(bc(i)%ni))) < 0.5d0) then
              nzero_chim = nzero_chim + 1
            endif
          endif

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

    if (nzero_chim > 0) then
      block
        use FUSS_Mod_MPI, only: mpi_abort_all
        character(len=256) :: msg
        write(msg,'(I0,A,I0,A)') nzero_chim, ' coupled chimera (104) entries in the level-', level, &
          ' BC file have MOSE donor volume fractions summing to ~0 (no donors found by the BC '// &
          'builder). Their ghost states would be garbage; fix the BC file before running.'
        call mpi_abort_all(trim(msg))
      end block
    endif

  end subroutine Read_BCfile


  subroutine Print_BC_Summary ()
    use FUSS_Config_Types_m, only: obj_ablation
    use FUSS_Lib_Bprime,     only: Bprime_Count, Bprime_Describe
    implicit none
    integer :: it

    write(*,*)
    write(*,'(A)') ' Boundary conditions'
    if (nconnect > 0) write(*,'(A,T35,I0)') '   Connection', nconnect
    if (nwall > 0)    write(*,'(A,T35,I0)') '   Viscous wall', nwall
    if (nsym > 0)     write(*,'(A,T35,I0)') '   Symmetry', nsym
    if (nchimera > 0) write(*,'(A,T35,I0)') '   Chimera', nchimera
    if (ncoupled > 0) write(*,'(A,T35,I0)') '   Coupled wall', ncoupled
    if (nablating > 0) then
      if (nablating > ncoupled_abl) write(*,'(A,T35,I0)') '   Ablating wall (306)', nablating - ncoupled_abl
      if (ncoupled_abl > 0) write(*,'(A,T35,I0)') '   Coupled ablating wall (307)', ncoupled_abl
      do it = 1, Bprime_Count()
        write(*,'(A)') '     B'' table '//trim(Bprime_Describe(it))
      enddo
    endif
    if (has_tdep_bc) write(*,'(A)') '   Time-dependent BC detected'

  end subroutine Print_BC_Summary


  !> The number of list-directed values on a line: blank-, tab- or comma-
  !> separated tokens (bare names hold none of these, plan 08 section 2.1).
  pure function Count_Tokens ( line ) result ( n )
    implicit none
    character(len=*), intent(in) :: line
    integer :: n, k
    logical :: intok
    n = 0 ; intok = .false.
    do k = 1, len_trim(line)
      if ( line(k:k) == ' ' .or. line(k:k) == ',' .or. line(k:k) == char(9) ) then
        intok = .false.
      elseif ( .not. intok ) then
        n = n + 1 ; intok = .true.
      endif
    enddo
  end function Count_Tokens

end module FUSS_IO_BC
