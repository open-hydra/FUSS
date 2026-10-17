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
    use FUSS_Field_List_m,        only: Field_List_Init
    use FUSS_Mod_Pyrolysis,       only: Pyrolysis_Assign_Setup
    use FUSS_Mod_Darcy,           only: Darcy_Assign_Setup
    implicit none

    ! Setting input solution
    call Setup_Input_Solution()

    ! Time
    if (obj_time_scheme%solver_type /= 'euler') then
      read(obj_time_scheme%solver_type(3:3), *) obj_time_scheme%n_rk
    else
      obj_time_scheme%n_rk = 1
    end if
    ! obj_irs%enabled is promoted from irs-beta in FUSS_Wrap_Setup, BEFORE
    ! Check_Input, so the guards see it. Do not move it back here.
    call Assign_Integration_Variables()

    ! Mesh motion. Derived here rather than registered directly, so that the
    ! INI exposes one readable keyword (`law`) instead of a bare on/off flag
    ! that could disagree with it. Must run before Setup_Data_Structure, which
    ! allocates blk%node_old only when this is true.
    obj_mesh_motion%enabled = ( trim(obj_mesh_motion%law) /= 'static' )

    ! Field list (plan 06 section 3.1): the pyrolysis rows need nsub, known
    ! once the material files are loaded (before this routine); the converters
    ! and the derived outflow direction come with them. Must run before
    ! Setup_Data_Structure, whose Allocate_Block walks the list.
    ! obj_darcy%enabled (gas-closure = type2) is derived in FUSS_Wrap_Setup
    ! right after Read_Inifile, before Check_Input, for the same reason as
    ! obj_irs%enabled: the guards must see it.
    call Field_List_Init( obj_pyrolysis%nsub, obj_pyrolysis%enabled, obj_darcy%enabled )
    if ( obj_pyrolysis%enabled ) call Pyrolysis_Assign_Setup()
    if ( obj_darcy%enabled )     call Darcy_Assign_Setup()

    ! Space: face metric of the interior diffusive fluxes
    obj_space_scheme%inverse_metric = ( trim(obj_space_scheme%diffusive_metric) /= 'mean' )

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

    ! Mesh motion
    if (obj_mesh_motion%enabled) then
      obj_mesh_motion%description = 'ALE, law = '//trim(obj_mesh_motion%law)
    else
      obj_mesh_motion%description = 'static mesh'
    end if

  end subroutine Assign_Setup

end module FUSS_Assign_Setup