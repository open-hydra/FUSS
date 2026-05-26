module FUSS_Load_Table
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Load_Table

contains

  subroutine Load_Table()
    use FUSS_Config_Types_m, only: obj_table
    use FUSS_Global_m
    use FUSS_Parameters_m
    use Lib_Tecplot
    use Lib_ORION_data
    use strings, only: parse
    implicit none
    ! Local
    character(len=llen)  :: wmfile, propertiesfile
    integer              :: ios, i, unitFile, n, Ti1, Ti2
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

    ios = tec_read_points_multivars(orion,4,trim(propertiesfile))
    if (ios/=0) then
      obj_table%error_message = '[ERROR] Properties file (properties.dat) not found'
    else
      Ti1 = nint(orion%block(1)%mesh(1,1,1,1))
      Ti2 = Ti1 + orion%block(1)%Ni - 1
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
    endif

  end subroutine Load_Table
  
end module FUSS_Load_Table