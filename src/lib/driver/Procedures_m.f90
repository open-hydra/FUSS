module FUSS_Procedures_m
  use FUSS_Wrap_Setup
  use FUSS_Wrap_Solve
  use FUSS_Wrap_Postprocess

  implicit none

  type :: FUSS_type

  contains
    procedure, nopass  :: setup => FUSS_setup
    procedure, nopass  :: solve => FUSS_solve
    procedure, nopass  :: postprocess => FUSS_postprocess
  end type FUSS_type

end module FUSS_Procedures_m