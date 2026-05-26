module FUSS_Assign_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ir_precision
  
  implicit none
  private
  public :: Assign_Setup

contains

  subroutine Assign_Setup()
    use FUSS_Config_Types_m
    use FUSS_IO_Solution,         only: Setup_Input_Solution
    use FUSS_Mod_Newstate,        only: Assign_Integration_Variables
    implicit none

    ! Setting input solution
    call Setup_Input_Solution()

    ! Time
    if (obj_time_scheme%solver_type /= 'euler') then
      read(obj_time_scheme%solver_type(3:3), *) obj_time_scheme%n_rk
    else
      obj_time_scheme%n_rk = 1
    end if
    if (obj_irs%beta>0d0) obj_irs%enabled = .true. 
    call Assign_Integration_Variables()

    !! Descriptions, warnings and errors

    ! Time scheme
    if (obj_time_scheme%solver_type == 'euler') then
      obj_time_scheme%description = 'Explicit Euler'
    else if (obj_time_scheme%solver_type == 'RK2') then
      obj_time_scheme%description = 'Second-order Runge-Kutta'
    else if (obj_time_scheme%solver_type == 'RK3') then
      obj_time_scheme%description = 'Third-order Runge-Kutta'
    end if
    if (obj_time_scheme%time_accurate) then
      obj_time_scheme%description = trim(obj_time_scheme%description)//' with time-accurate switch enabled'
    end if
    if (obj_irs%enabled) then
      obj_irs%description = 'Beta set to '//trim(str(.true.,real(obj_irs%beta)))
    end if

  end subroutine Assign_Setup

end module FUSS_Assign_Setup