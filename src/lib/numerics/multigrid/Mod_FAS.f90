module FUSS_Mod_FAS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: FAS_Solve

contains

  !> One FAS cycle rooted at the running level (grid sequencing keeps handing finer roots: FMG).
  !> The root level's smoothing iterations are normal reported iterations; every coarser visit is
  !> silent (inner). The cycle stops early when the root declares the end of the run or of the level.
  subroutine FAS_Solve ( simulation, External_Function )
    use FUSS_Advanced_Types_m, only: FUSS_simulation_type
    use FUSS_Config_Types_m,   only: obj_multigrid
    implicit none
    type(FUSS_simulation_type), intent(inout) :: simulation
    external :: External_Function
    ! Local
    integer :: top

    top = obj_multigrid%MG_level
    call FAS_Level ( simulation, External_Function, top, top )
    obj_multigrid%MG_level = top   ! the postprocess works on the root level

  end subroutine FAS_Solve


  recursive subroutine FAS_Level ( simulation, External_Function, level, top )
    use FUSS_Advanced_Types_m, only: FUSS_simulation_type
    use FUSS_Config_Types_m,   only: obj_multigrid, obj_sim_param
    use FUSS_Mod_Explicit,     only: Explicit_Step
    implicit none
    type(FUSS_simulation_type), intent(inout) :: simulation
    external :: External_Function
    integer, intent(in) :: level, top
    ! Local
    integer :: it, g, nsmooth, ngamma

    ! Pre-smoothing; the coarsest level of a cycle takes its own budget. A root that is also
    ! the coarsest level (the first stint of the FMG descent) has no cycle below and smooths
    ! the whole pre+post budget here, never zero: pre-iter = 0 (legal, post-iter >= 1) would
    ! otherwise freeze the level and the never-allocated residual norm crashes the first
    ! residual-history write.
    nsmooth = obj_multigrid%pre_iter
    if (level == obj_multigrid%MGL) then
      if (level /= top) then
        nsmooth = obj_multigrid%coarse_iter
      else
        nsmooth = max(obj_multigrid%pre_iter + obj_multigrid%post_iter, 1)
      endif
    endif

    obj_multigrid%MG_level = level
    do it = 1, nsmooth
      call Explicit_Step ( simulation%domain, External_Function, inner = (level /= top) )
      if (level == top) then
        if (obj_sim_param%TODO == 3 .or. obj_multigrid%change_MG) return
      endif
    enddo

    if (level == obj_multigrid%MGL) return

    call FAS_To_Coarse ( simulation, level )

    ngamma = 1
    if (trim(obj_multigrid%cycle) == 'w') ngamma = 2
    do g = 1, ngamma
      call FAS_Level ( simulation, External_Function, level+1, top )
    enddo

    call FAS_Correct ( simulation, level )

    ! Post-smoothing
    obj_multigrid%MG_level = level
    do it = 1, obj_multigrid%post_iter
      call Explicit_Step ( simulation%domain, External_Function, inner = (level /= top) )
      if (level == top) then
        if (obj_sim_param%TODO == 3 .or. obj_multigrid%change_MG) return
      endif
    enddo

  end subroutine FAS_Level


  !> Restrict the state, remember it in h, and build the forcing of level+1:
  !> F = sum of the children's volume-integrated steady residuals (which carry the finer level's
  !> own forcing) minus the coarse one at the restricted state.
  subroutine FAS_To_Coarse ( simulation, level )
    use FUSS_Advanced_Types_m, only: FUSS_simulation_type
    use FUSS_Mod_Explicit,     only: Evaluate_Residual
    use FUSS_Mod_Multigrid,    only: Restriction
    use FUSS_Lib_Multigrid,    only: fas_store_h0, fas_build_forcing
    use FUSS_Mod_MPI,          only: is_local_block
    implicit none
    type(FUSS_simulation_type), intent(inout) :: simulation
    integer, intent(in) :: level
    ! Local
    integer :: b

    call Evaluate_Residual ( simulation%domain(level) )

    simulation%domain(level+1)%fas_forcing = .false.
    call Restriction ( Fine=simulation%domain(level), Coarse=simulation%domain(level+1) )
    do b = 1, simulation%domain(level+1)%nb
      if (.not. is_local_block(b)) cycle
      call fas_store_h0 ( simulation%domain(level+1)%blk(b)%T, simulation%domain(level+1)%blk(b)%matID, &
                          simulation%domain(level+1)%blk(b)%fas_h0, simulation%domain(level+1)%blk(b)%dim )
    enddo

    call Evaluate_Residual ( simulation%domain(level+1) )

    do b = 1, simulation%domain(level+1)%nb
      if (.not. is_local_block(b)) cycle
      call fas_build_forcing ( simulation%domain(level  )%blk(b)%R, simulation%domain(level  )%blk(b)%qvol, &
                               simulation%domain(level  )%blk(b)%vol, &
                               simulation%domain(level+1)%blk(b)%R, simulation%domain(level+1)%blk(b)%qvol, &
                               simulation%domain(level+1)%blk(b)%vol, &
                               simulation%domain(level+1)%blk(b)%fas_f, &
                               simulation%domain(level  )%blk(b)%dim, simulation%domain(level+1)%blk(b)%dim )
    enddo
    simulation%domain(level+1)%fas_forcing = .true.

  end subroutine FAS_To_Coarse


  !> Prolong the coarse correction (in h) onto level, then drop the coarse forcing
  subroutine FAS_Correct ( simulation, level )
    use FUSS_Advanced_Types_m, only: FUSS_simulation_type
    use FUSS_Lib_Multigrid,    only: coarse2fine_corr
    use FUSS_Mod_MPI,          only: is_local_block
    implicit none
    type(FUSS_simulation_type), intent(inout) :: simulation
    integer, intent(in) :: level
    ! Local
    integer :: b

    do b = 1, simulation%domain(level)%nb
      if (.not. is_local_block(b)) cycle
      call coarse2fine_corr ( simulation%domain(level  )%blk(b)%T,  simulation%domain(level  )%blk(b)%matID, &
                              simulation%domain(level+1)%blk(b)%T,  simulation%domain(level+1)%blk(b)%matID, &
                              simulation%domain(level+1)%blk(b)%fas_h0, &
                              simulation%domain(level  )%blk(b)%dim, simulation%domain(level+1)%blk(b)%dim )
    enddo
    simulation%domain(level+1)%fas_forcing = .false.

  end subroutine FAS_Correct

end module FUSS_Mod_FAS
