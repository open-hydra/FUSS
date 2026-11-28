!> The checked interface of the per-stage user hook.
!>
!> Explicit_Step calls the hook once per Runge-Kutta stage as
!> External_Function(domain), from inside its OpenMP parallel region. Until
!> this module existed the dummy was declared `external`, so a procedure with
!> the wrong signature (hydra's argument-less Dummy_Function, for one) was
!> accepted silently and happened to be harmless only because every body
!> ignored its arguments. Declaring the dummy as procedure(external_function_if)
!> in FUSS_solve and Explicit_Step makes the mismatch a compile-time error;
!> test/run_hook_interface_check.sh demonstrates both outcomes.
!>
!> The interface lives in its own module because FUSS_Wrap_Solve already uses
!> FUSS_Mod_Explicit: importing it from either of those would be circular.
module FUSS_Hook_Interface_m
  implicit none
  private
  public :: external_function_if

  abstract interface
    subroutine external_function_if ( domain )
      use FUSS_Advanced_Types_m, only: FUSS_domain_type
      implicit none
      type(FUSS_domain_type), intent(inout) :: domain
    end subroutine external_function_if
  end interface

end module FUSS_Hook_Interface_m
