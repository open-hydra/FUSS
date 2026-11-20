module FUSS_Load_Table
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Load_Table, Read_Points_Table

contains

  subroutine Load_Table()
    use FUSS_Config_Types_m, only: obj_table
    use FUSS_Global_m
    use FUSS_Parameters_m
    use Lib_Tecplot
    use Lib_ORION_data
    use strings, only: parse
    use FUSS_Finite_m, only: Not_Finite
    implicit none
    ! Local
    character(len=llen)  :: wmfile, propertiesfile
    integer              :: ios, i, unitFile, n, Ti1, Ti2, Tn, Tdev, c
    real(R8)             :: cap, dev, devmax, val
    character(len=12)    :: column
    character(llen)      :: wholestring, args(2)
    type(ORION_data)     :: orion

    !! ------------------------------------------------------
    !! Properties Table -------------------------------------
    !! ------------------------------------------------------
    obj_table%warning_message = 'none'
    obj_table%error_message   = 'none'
    obj_table%description     = 'none'

    wmfile = 'INPUT/'//trim(FUSS_phase_prefix)//'phase.txt'
    propertiesfile = 'INPUT/'//trim(FUSS_phase_prefix)//'properties.dat'

    open(newunit=unitFile,file=trim(wmfile),status='old',iostat=ios)
    if (ios/=0) then
      obj_table%error_message = '[ERROR] Phase file (phase.txt) not found'
    else
      ios = 0; n = -1
      read(unitFile,*) ! skip first line
      do while(ios==0)
        read(unitFile,'(A)',iostat=ios) wholestring
        n = n + 1
      enddo
      obj_table%n = n
      allocate(obj_table%name(1:n))
      rewind(unitFile)

      read(unitFile,*) !skip first line
      do i = 1, n
        read(unitFile,'(A)') wholestring
        call parse(wholestring,' ',args)
        obj_table%name(i) = trim(adjustl(args(1)))
      end do
      close(unitFile)
    endif

    ios = Read_Points_Table(orion,4,trim(propertiesfile))
    if (ios/=0) then
      obj_table%error_message = '[ERROR] Properties file (properties.dat) not found'
    else
      Ti1 = nint(orion%block(1)%mesh(1,1,1,1))
      Ti2 = Ti1 + orion%block(1)%Ni - 1
      obj_table%Tmin = Ti1
      obj_table%Tmax = Ti2
      allocate(obj_table%cp(1:obj_table%n,Ti1:Ti2))
      allocate(obj_table%rho(1:obj_table%n,Ti1:Ti2))
      allocate(obj_table%kappa(1:obj_table%n,Ti1:Ti2))
      allocate(obj_table%h(1:obj_table%n,Ti1:Ti2))
      do i = 1, obj_table%n
        obj_table%cp(i,Ti1:Ti2)    = orion%block(i)%vars(1,:,1,1)
        obj_table%rho(i,Ti1:Ti2)   = orion%block(i)%vars(2,:,1,1)
        obj_table%kappa(i,Ti1:Ti2) = orion%block(i)%vars(3,:,1,1)
        obj_table%h(i,Ti1:Ti2)     = orion%block(i)%vars(4,:,1,1)
      enddo

      ! Every Cp, Density, Conductivity and Energy of the table is finite, Cp and Density positive, Conductivity not
      ! negative (0 makes a material that conducts no heat): a value outside is an error of the file, found here and not
      ! at the first cell whose temperature reaches its row (or never, for a row no cell reaches). Not_Finite reads the
      ! bits: the sign is compared only on a finite value, as a build that traps the invalid operations requires
      table_check: do i = 1, obj_table%n
        do Tn = Ti1, Ti2
          do c = 1, 4
            select case (c)
            case (1); val = obj_table%cp(i,Tn);    column = 'Cp'
            case (2); val = obj_table%rho(i,Tn);   column = 'Density'
            case (3); val = obj_table%kappa(i,Tn); column = 'Conductivity'
            case (4); val = obj_table%h(i,Tn);     column = 'Energy'
            end select
            if ( Not_Finite(val) ) then
              write(obj_table%error_message,'(5A,I0,3A)') '[ERROR] ', trim(propertiesfile), ', material ', &
                trim(obj_table%name(i)), ', row T = ', Tn, ' K: ', trim(column), ' is not finite (NaN or infinity)'
              exit table_check
            endif
            if ( (c <= 2 .and. val <= 0d0) .or. (c == 3 .and. val < 0d0) ) then
              write(obj_table%error_message,'(5A,I0,3A,ES10.3,A)') '[ERROR] ', trim(propertiesfile), ', material ', &
                trim(obj_table%name(i)), ', row T = ', Tn, ' K: ', trim(column), ' = ', val, &
                merge(' must be > 0 ', ' must be >= 0', c <= 2)
              exit table_check
            endif
          enddo
        enddo
      enddo table_check
      ! A table with an error stops at set-up (Check_Table): neither the order of its Energy column nor the warning below
      if ( index(obj_table%error_message,'ERROR') > 0 ) return

      ! co_T bisects a non-decreasing Energy column (a NaN entry makes it a scan, as any other column)
      allocate(obj_table%h_monotone(1:obj_table%n))
      do i = 1, obj_table%n
        obj_table%h_monotone(i) = all( obj_table%h(i,Ti1+1:Ti2) >= obj_table%h(i,Ti1:Ti2-1) )
      enddo

      ! The Energy column must be h = int rho cp dT (cons advances h, prim divides by rho cp): warn if not
      do i = 1, obj_table%n
        devmax = 0d0; Tdev = Ti1
        do Tn = Ti1, Ti2 - 1
          cap = 0.5d0 * ( obj_table%rho(i,Tn)*obj_table%cp(i,Tn) + obj_table%rho(i,Tn+1)*obj_table%cp(i,Tn+1) )
          dev = abs( obj_table%h(i,Tn+1) - obj_table%h(i,Tn) - cap ) / cap
          if ( dev > devmax ) then
            devmax = dev; Tdev = Tn
          endif
        enddo
        if ( devmax > 1d-2 ) then
          write(obj_table%warning_message,'(A,A,A,F0.1,A,I0,A)') '[WARNING] properties.dat, material ', &
            trim(obj_table%name(i)), ': the Energy column differs from int rho cp dT (dh/dT off rho*cp by ', &
            100d0*devmax, ' % at ', Tdev, ' K); integration-variables = cons uses it as given'
        endif
      enddo
    endif

  end subroutine Load_Table
  
  ! A Tecplot table of points (a zone per material, a row per temperature) in the layout of ORION's former
  ! tec_read_points_multivars: block b holds Ni = N rows, Nj = Nk = 1, mesh(1,1:N,1,1) (the first column) and
  ! vars(1:nvars,1:N,1,1). Non-zero when ORION cannot read the file or a zone does not hold one coordinate and nvars
  ! variables at its points; orion is then left as ORION's reader left it
  function Read_Points_Table ( orion, nvars, filename ) result(ios)
    use IR_Precision,   only: R8P
    use Lib_Tecplot,    only: tec_read_structured_multiblock
    use Lib_ORION_data, only: ORION_data
    implicit none
    type(ORION_data), intent(inout) :: orion
    integer,          intent(in)    :: nvars
    character(len=*), intent(in)    :: filename
    integer                         :: ios
    ! Local
    real(R8P), allocatable :: mesh(:,:,:,:), vars(:,:,:,:)
    integer                :: b, n

    ! ORION's reader gives a line of N points as Ni = N-1, mesh(1:1,0:N-1,0:0,0:0), vars(1:nvar,0:N-1,0:0,0:0)
    ios = tec_read_structured_multiblock(orion=orion, filename=filename)
    if ( ios /= 0 ) return
    ios = 1
    do b = 1, size(orion%block)
      if ( size(orion%block(b)%mesh,1) /= 1 .or. size(orion%block(b)%mesh) /= size(orion%block(b)%mesh,2) ) return
      if ( .not. allocated(orion%block(b)%vars) ) return
      if ( size(orion%block(b)%vars,1) /= nvars .or. size(orion%block(b)%vars) /= nvars*size(orion%block(b)%mesh) ) return
    enddo
    do b = 1, size(orion%block)
      n = size(orion%block(b)%mesh)
      allocate( mesh(1,n,1,1), vars(nvars,n,1,1) )
      mesh = reshape( orion%block(b)%mesh, shape(mesh) )
      vars = reshape( orion%block(b)%vars, shape(vars) )
      call move_alloc( mesh, orion%block(b)%mesh )
      call move_alloc( vars, orion%block(b)%vars )
      orion%block(b)%Ni = n
      orion%block(b)%Nj = 1
      orion%block(b)%Nk = 1
    enddo
    ios = 0

  end function Read_Points_Table

end module FUSS_Load_Table