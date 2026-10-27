module FUSS_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_BC
  public :: Print_BC_Summary

  ! BC type counters — populated by Read_BCfile, consumed by Print_BC_Summary
  integer :: nconnect, nwall, nsym, nchimera, ncoupled
  logical, public, protected :: has_tdep_bc = .false.   ! a 301/302 BC reads a time table
  character(len=64) :: bc_name = 'bc.txt'               ! the file being read, for the messages of its errors
  real(R8), parameter :: EMPTY  = -huge(1.0_R8)         ! a value of a property line that the read has not given
  integer,  parameter :: IEMPTY = -huge(1_I4)

contains

  subroutine Setup_BC ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_multigrid, obj_io_bc
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    ! Local
    integer :: m, error

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

    ! Check BC file consistency
    ios = 0; n_proof = -1
    do while (ios==0)
      read( unitfile,*,iostat=ios ) di(1), di(2), di(3), di(4), di(5), ti
      select case(ti)
      case(101,301,302,303,304,305,103)
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
    use FUSS_Config_Types_m, only: obj_io, obj_io_bc
    use FUSS_Global_m
    use IR_Precision
    implicit none
    type(FUSS_bc_type), dimension(:), intent(inout) :: bc
    integer, intent(in)                             :: level
    integer, dimension(1:,1:), intent(inout)        :: n_bf
    ! Local
    integer :: cc, i, s
    integer :: unitfile, ios, cios
    integer :: nline                         ! lines of the file read so far (Next_Line)
    integer :: nzero_chim
    character(len=256) :: BCfile
    character(len=:), allocatable :: pline   ! a line of the file, read whole
    integer :: lc

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
    if (level == 1) then
      bc_name = trim(FUSS_phase_prefix)//'bc.txt'
    else
      bc_name = trim(FUSS_phase_prefix)//'bc'//trim(str(.true.,level))//'.txt'
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
    nzero_chim = 0

    ! Read file: every line is read whole (Next_Line, which counts the lines of the file and skips the blank ones as a
    ! list-directed read does), then as its values; a line that is short or holds a value that is not a number stops
    ! every rank (BC_Line_Error: every rank reads the same file)
    nline = 0
    do i = 1, size(bc)

      ! First line is equal for every BC type
      call Next_Line ( unitfile, pline, ios, nline )
      if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%b, bc(i)%i, bc(i)%j, bc(i)%k, bc(i)%f, bc(i)%type
      if ( ios /= 0 ) call BC_Line_Error ( i, -1, nline, pline )

      ! n_bf update
      n_bf( bc(i) % b, bc(i) % f ) = n_bf( bc(i) % b, bc(i) % f ) + 1

      ! Second line is BC type-dependent
      select case( bc(i)%type )
        case(101)
          if (level == 1) nconnect = nconnect + 1
          bc(i)%bs = IEMPTY ; bc(i)%is = IEMPTY ; bc(i)%js = IEMPTY ; bc(i)%ks = IEMPTY ; bc(i)%fs = IEMPTY
          bc(i)%d11 = IEMPTY ; bc(i)%d12 = IEMPTY ; bc(i)%d21 = IEMPTY ; bc(i)%d22 = IEMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
          call Need_I ( i, bc(i)%type, nline, [ bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, &
                                                bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22 ], 'bs is js ks fs d11 d12 d21 d22' )

        case(300)
          if (level == 1) nsym = nsym + 1
        
        case(301)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
          BCfile = adjustl(pline)
          lc = len_trim(BCfile)
          if (lc > 0) then
            if (BCfile(lc:lc) == ',') BCfile(lc:lc) = ' '   ! the data lines end with a comma
          endif
          ! Time-varying properties: if the read string is not convertible to a real, the BC is time-dependent
          read(BCfile,*,iostat=cios) bc(i)%qw
          if (cios/=0) then
            call bc(i)%BCtime%initialize(file=BCfile)
            has_tdep_bc = .true.
          endif
        
        case(302)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
          BCfile = adjustl(pline)
          lc = len_trim(BCfile)
          if (lc > 0) then
            if (BCfile(lc:lc) == ',') BCfile(lc:lc) = ' '   ! the data lines end with a comma
          endif
          ! Time-varying properties: if the read string is not convertible to a real, the BC is time-dependent
          read(BCfile,*,iostat=cios) bc(i)%Tw
          if (cios/=0) then
            call bc(i)%BCtime%initialize(file=BCfile)
            has_tdep_bc = .true.
          endif

        case(303)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%hconv = EMPTY ; bc(i)%qrad = EMPTY ; bc(i)%Tref = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios )  bc(i)%hconv, bc(i)%qrad, bc(i)%Tref
          if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
          call Need_R ( i, bc(i)%type, nline, [ bc(i)%hconv, bc(i)%qrad, bc(i)%Tref ], 'hconv qrad Tref' )
        
        case(304)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%eps = EMPTY ; bc(i)%Tref = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios )  bc(i)%eps, bc(i)%Tref
          if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
          call Need_R ( i, bc(i)%type, nline, [ bc(i)%eps, bc(i)%Tref ], 'eps Tref' )

        case(305)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%hconv = EMPTY ; bc(i)%eps = EMPTY ; bc(i)%Tref = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios )  bc(i)%hconv, bc(i)%eps, bc(i)%Tref
          if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
          call Need_R ( i, bc(i)%type, nline, [ bc(i)%hconv, bc(i)%eps, bc(i)%Tref ], 'hconv eps Tref' )

        case(102,104)
          if (level == 1 .and. bc(i)%type == 102) nchimera = nchimera + 1
          bc(i)%ni = IEMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) (bc(i)%ni(cc),cc=1,2)
          if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
          call Need_I ( i, bc(i)%type, nline, bc(i)%ni, 'n1 n2' )
          if ( any(bc(i)%ni < 0) ) call BC_Value_Error ( i, bc(i)%type, nline, pline, 'the numbers of donors are negative' )
          allocate(bc(i)%donorID(1:sum(bc(i)%ni),1:4))
          allocate(bc(i)%volume_fraction(1:sum(bc(i)%ni)))
          bc(i)%donorID = IEMPTY
          bc(i)%volume_fraction = EMPTY
          do s = 1, bc(i)%ni(1)+bc(i)%ni(2)
            call Next_Line ( unitfile, pline, ios, nline )
            if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
            if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
            call Need_I ( i, bc(i)%type, nline, bc(i)%donorID(s,1:4), 'bD iD jD kD' )
            call Need_R ( i, bc(i)%type, nline, bc(i)%volume_fraction(s:s), 'vol_frac' )
          enddo
          if (bc(i)%type == 104) then
            ! Coupled wall as 103 (coupling flux, wall output); the donors are cells of the other solver
            if (level == 1) ncoupled = ncoupled + 1
            obj_io_bc%coupling_flag( bc(i)%b , bc(i)%f ) = .true.
            bc(i)%ext_flux = 0.0
          endif
          if (sum(bc(i)%volume_fraction(1:bc(i)%ni(1))) < 0.5d0 .or. &
              sum(bc(i)%volume_fraction(bc(i)%ni(1)+1:sum(bc(i)%ni))) < 0.5d0) then
            nzero_chim = nzero_chim + 1
          endif

        case(103)
          if (level == 1) ncoupled = ncoupled + 1
          obj_io_bc%coupling_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%bs = IEMPTY ; bc(i)%is = IEMPTY ; bc(i)%js = IEMPTY ; bc(i)%ks = IEMPTY ; bc(i)%fs = IEMPTY
          bc(i)%d11 = IEMPTY ; bc(i)%d12 = IEMPTY ; bc(i)%d21 = IEMPTY ; bc(i)%d22 = IEMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
          call Need_I ( i, bc(i)%type, nline, [ bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, &
                                                bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22 ], 'bs is js ks fs d11 d12 d21 d22' )
          bc(i)%ext_flux = 0.0
      
      end select
      
    enddo

    if (nwall > 0) obj_io % write_wall = .true.

    close( unitfile )

    if (nzero_chim > 0) then
      block
        use FUSS_Mod_MPI, only: mpi_stop_all
        character(len=256) :: msg
        write(msg,'(I0,A,I0,A)') nzero_chim, ' chimera (102/104) entries in the level-', level, &
          ' BC file have donor volume fractions summing below 0.5 (no donors found by the BC '// &
          'builder). Their ghost states would be garbage; fix the BC file before running.'
        call mpi_stop_all(trim(msg))
      end block
    endif

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


  !─────────────────────────────────────────────────────────────────────────────
  ! Next_Line: the next line of the boundary file that is not blank (Read_Line), as a list-directed read skips
  !   blank lines; nline counts every line read, blank ones included, so it is the line number in the file.
  subroutine Next_Line ( unit, line, ios, nline )
    implicit none
    integer, intent(in)                        :: unit
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out)                       :: ios
    integer, intent(inout)                     :: nline

    do
      call Read_Line ( unit, line, ios )
      if ( ios /= 0 ) exit
      nline = nline + 1
      if ( len_trim(line) > 0 ) exit
    end do

  end subroutine Next_Line


  !─────────────────────────────────────────────────────────────────────────────
  ! Read_Line: the next line of a formatted file, whole, whatever its length (read in pieces without advancing).
  !   ios < 0 at the end of the file; a last line without its end-of-line mark is still a line.
  subroutine Read_Line ( unit, line, ios )
    implicit none
    integer, intent(in)                        :: unit
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out)                       :: ios
    ! Local
    character(len=256) :: piece
    integer :: nread

    line = ''
    do
      read ( unit, '(A)', advance='no', iostat=ios, size=nread ) piece
      line = line // piece(1:nread)
      if ( ios /= 0 ) exit
    end do
    if ( is_iostat_eor(ios) ) ios = 0
    if ( is_iostat_end(ios) .and. len(line) > 0 ) ios = 0

  end subroutine Read_Line


  !─────────────────────────────────────────────────────────────────────────────
  ! BC_Line_Error: stop every rank on line nline of the boundary file, of record n, that cannot be read as the values it
  !   must hold (short, or a value that is not a number): its record line (btype < 0) or a property line of type btype;
  !   an empty line is the end of the file. Every rank reads the file, so every rank reaches this stop.
  subroutine BC_Line_Error ( n, btype, nline, line )
    use FUSS_Mod_MPI, only: mpi_stop_all
    implicit none
    integer, intent(in)          :: n, btype, nline
    character(len=*), intent(in) :: line
    ! Local
    character(len=160) :: head

    if ( len_trim(line) == 0 ) then
      write(head,'(3A,I0,A,I0)') '[ERROR] ', trim(bc_name), ' ends at line ', nline, ', before the end of record ', n
      call mpi_stop_all ( trim(head) )
    end if
    if ( btype < 0 ) then
      write(head,'(3A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, &
        '): the line is not block, i, j, k, face and type'
    else
      write(head,'(3A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, ', type ', btype, &
        '): the line is short or holds a value that is not a number'
    end if
    call mpi_stop_all ( trim(head)//': "'//trim(line(1:min(len(line),200)))//'"' )

  end subroutine BC_Line_Error


  !─────────────────────────────────────────────────────────────────────────────
  ! BC_Value_Error: stop every rank on a value of line nline of record n that reads but is not valid (why).
  subroutine BC_Value_Error ( n, btype, nline, line, why )
    use FUSS_Mod_MPI, only: mpi_stop_all
    implicit none
    integer, intent(in)          :: n, btype, nline
    character(len=*), intent(in) :: line, why
    ! Local
    character(len=160) :: head

    write(head,'(3A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, ', type ', btype, '): '
    call mpi_stop_all ( trim(head)//trim(why)//': "'//trim(line(1:min(len(line),200)))//'"' )

  end subroutine BC_Value_Error


  !─────────────────────────────────────────────────────────────────────────────
  ! Need_R, Need_I: stop every rank on the first of the values v of record n, type btype, line nline of the file, still
  !   EMPTY after the read of its property line (a null value, ",,"). names: their names, blank-separated, in the order
  !   of the line; the last name stands for the values after it.
  subroutine Need_R ( n, btype, nline, v, names )
    implicit none
    integer, intent(in)          :: n, btype, nline
    real(R8), intent(in)         :: v(:)
    character(len=*), intent(in) :: names
    ! Local
    integer :: k

    do k = 1, size(v)
      if ( v(k) == EMPTY ) call Empty_Field ( n, btype, nline, Name_Of ( names, k ) )
    end do

  end subroutine Need_R


  subroutine Need_I ( n, btype, nline, v, names )
    implicit none
    integer, intent(in)          :: n, btype, nline
    integer, intent(in)          :: v(:)
    character(len=*), intent(in) :: names
    ! Local
    integer :: k

    do k = 1, size(v)
      if ( v(k) == IEMPTY ) call Empty_Field ( n, btype, nline, Name_Of ( names, k ) )
    end do

  end subroutine Need_I


  !─────────────────────────────────────────────────────────────────────────────
  ! Empty_Field: stop every rank on the field what of the property line of record n, type btype, line nline of the
  !   file, left empty.
  subroutine Empty_Field ( n, btype, nline, what )
    use FUSS_Mod_MPI, only: mpi_stop_all
    implicit none
    integer, intent(in)          :: n, btype, nline
    character(len=*), intent(in) :: what
    ! Local
    character(len=240) :: msg

    write(msg,'(3A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, ', type ', btype, &
      '): the value of '//trim(what)//' is empty'
    call mpi_stop_all ( trim(msg) )

  end subroutine Empty_Field


  !─────────────────────────────────────────────────────────────────────────────
  ! Name_Of: the k-th blank-separated word of names, or the last one when there are fewer.
  pure function Name_Of ( names, k ) result ( w )
    implicit none
    character(len=*), intent(in) :: names
    integer, intent(in)          :: k
    character(len=len(names))    :: w
    ! Local
    integer :: p, q, m

    w = ' '
    p = 1
    m = 0
    do
      do while ( p <= len_trim(names) )
        if ( names(p:p) /= ' ' ) exit
        p = p + 1
      end do
      if ( p > len_trim(names) ) exit
      q = index ( names(p:), ' ' )
      if ( q == 0 ) then
        q = len(names) + 1
      else
        q = p + q - 1
      end if
      m = m + 1
      w = names(p:q-1)
      if ( m == k ) exit
      p = q
    end do

  end function Name_Of


end module FUSS_IO_BC
