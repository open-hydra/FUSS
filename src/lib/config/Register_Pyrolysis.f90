!> Registry entries of the [FUSS-Pyrolysis] section (plan 06 §3.5).
!>
!> The section carries only run-level scalars. Every per-phase and per-reaction
!> constant lives in INPUT/<prefix>pyrolysis.dat (loaded by FUSS_Load_Pyrolysis),
!> so no array key has to be pre-sized here. The feature switch is the PRESENCE
!> of the section, detected in Scan_Ini before the registry exists; the keys
!> below are registered unconditionally, so an absent section leaves every one
!> of them at its default and the runtime re-read of input.ini stays silent.
!> All four are immutable mid-run (none is in RUNTIME_MUTABLE).
module FUSS_Read_Pyrolysis
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Config_Types_m
  use FUSS_Input_Registry

  implicit none
  private
  public :: Register_Pyrolysis

contains

  subroutine Register_Pyrolysis ()
    use FUSS_Parameters_m
    implicit none
    character(len=llen) :: section

    section = trim(codename)//'-Pyrolysis'

    call reg%add( trim(section), 'gas-closure', obj_pyrolysis%gas_closure, 'type1', &
                  'Pyrolysis-gas closure: type1 = quasi-steady 1-D outflow (type2 = Darcy, plan 07)', &
                  'type1 ,  type2', .false. )

    call reg%add( trim(section), 'outflow-direction', obj_pyrolysis%outflow_direction, 'none', &
                  'Index direction along which the pyrolysis gas leaves the (single) block', &
                  'none ,  +i ,  -i ,  +j ,  -j ,  +k ,  -k', .false. )

    call reg%add( trim(section), 'p-gas', obj_pyrolysis%p_gas, '101325.0', &
                  'Prescribed pyrolysis-gas pressure [Pa] at which the gas table is evaluated (type1)', &
                  '> 0', .false. )

    call reg%add( trim(section), 'hbar-guard-report', obj_pyrolysis%hbar_guard_report, '.false.', &
                  'Print at setup the materials whose hbar guard fired (rho_v = rho_c, no net decomposition)', &
                  'logical', .false. )

  end subroutine Register_Pyrolysis

end module FUSS_Read_Pyrolysis
