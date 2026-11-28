program FUSS_program
#if defined (_OPENMP)
  use omp_lib
#endif
  use FUSS_Advanced_Types_m, only: FUSS_simulation_type
  use FUSS_Config_Types_m,   only: obj_sim_param
  use FUSS_Procedures_m,     only: FUSS_type
  use FUSS_Mod_MPI
#ifdef USE_MPI
  use FUSS_Mod_GhostExchange, only: cleanup_ghost_schedule
#endif
  implicit none
  type(FUSS_type)            :: FUSS
  type(FUSS_simulation_type) :: simulation

  ! Initialize MPI environment (no-op if USE_MPI is not defined)
  call mpi_init_env()

#if defined (_OPENMP)
  !$omp parallel
  obj_sim_param%nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  if (mpi_is_root) then
    write(*,'(A)')    ' Parallel execution'
    write(*,'(A)')    ' OpenMP:'
    write(*,'(A,I4)') ' -  Number of threads --> ', obj_sim_param%nthreads
  end if
#else
  if (mpi_is_root) write(*,'(A)')    ' Serial execution'
  obj_sim_param%nthreads = 1
#endif

#ifdef USE_MPI
  if (mpi_is_root) then
    write(*,'(A)')      ' MPI:'
    write(*,'(A,I4)')   ' -  Number of ranks   --> ', mpi_size_
  end if
#endif

  ! Solving with FUSS
  call FUSS%setup( simulation )

  obj_sim_param%TODO = 1
  do while ( obj_sim_param%TODO <= 2 )
    call FUSS%solve( simulation, Dummy_Function )
    if ( obj_sim_param%TODO <= 2 ) call FUSS%postprocess( simulation )
  enddo
  
  call FUSS%postprocess( simulation )

    ! Free persistent MPI requests before finalizing
#ifdef USE_MPI
  call cleanup_ghost_schedule()
#endif

  ! Finalize MPI environment
  call mpi_finalize_env()

contains

    !> User hook invoked by the solver once per Runge-Kutta stage.
    !!
    !! Explicit_Step calls this as External_Function(domain) -- see
    !! Mod_Explicit.f90 -- so it MUST accept the domain argument. It previously
    !! took none, which made the call non-conforming (it happened to be
    !! harmless only because the dummy body ignores everything).
    !!
    !! Two constraints on any real implementation:
    !!   1. It is called from INSIDE an !$omp parallel region (opened in
    !!      Explicit_Step), so use !$omp do -- never !$omp parallel do.
    !!   2. It is called once per RK stage, not once per time step.
    !!
    !! The solver declares this hook as procedure(external_function_if)
    !! (src/lib/base/Hook_Interface_m.f90), so a wrong signature is a compile
    !! error; test/run_hook_interface_check.sh demonstrates it. hydra's single
    !! argument-less Dummy_Function, passed to both FUSS%solve and MOSE%solve,
    !! must be split in two when hydra moves its FUSS pointer past this commit.
    !!
    !! When pyrolysis is enabled the stage work (Pyrolysis_Stage) has already
    !! added the gas-advection flux to blk%R and taken the energy-budget sums
    !! from it before this hook runs: a hook must not modify R in that case.
    subroutine Dummy_Function ( domain )
      use FUSS_Advanced_Types_m, only: FUSS_domain_type
      implicit none
      type(FUSS_domain_type), intent(inout) :: domain
      ! Deliberately empty: no user-defined operations in the standalone solver.
      ! Reference the argument so compilers do not warn about it being unused.
      if (.false.) domain%iter = domain%iter
    end subroutine Dummy_Function

end program FUSS_program