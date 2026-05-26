module FUSS_Mod_Explicit
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Explicit_Step

contains

  subroutine Explicit_Step ( domain, External_Function )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m
    use FUSS_Global_m
    use FUSS_Mod_dt,         only: Set_Global_dt, Compute_dt
    use FUSS_Lib_Ghost,      only: Fill_Ghost_Cell
    use FUSS_Mod_Fluxes,     only: Fluxes
    use FUSS_Mod_BC_Fluxes,  only: BC_Fluxes
    use FUSS_Mod_Newstate,   only: RK_Newstate
    use FUSS_Mod_Diagnostic, only: Compute_Residual
    use FUSS_Mod_MPI, only: is_local_block, mpi_reduce_sum_r8, &
                           mpi_is_root, mpi_bcast_logical, mpi_bcast_integer
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    external :: External_Function
    ! Local
    logical  :: endsim, iosim, endmg
    integer  :: i_rk, b, level
    real(R8) :: average

    level = obj_multigrid%MG_level
    obj_multigrid%change_MG = .false.
    domain(level) % iter = domain(level) % iter + 1
    obj_sim_param%iter_from_call = obj_sim_param%iter_from_call + 1
    obj_sim_param%iter_general   = obj_sim_param%iter_general + 1

    if (obj_sim_param%HYDRA_time_accurate) then
      call Set_Global_dt ( domain(level) )
      domain(level) % time = domain(level) % time + domain(level) % dtglobal
    else
      domain(level) % dtglobal = 1d5
      call Compute_dt ( domain(level), obj_time_scheme%vnn, obj_time_scheme%rampa_vnn_iter )  ! Compute local and minimum time step
      if ( obj_time_scheme%time_accurate ) then
        call Set_Global_dt ( domain(level) )  ! Time-accurate: apply global minimum time step
        domain(level) % time = domain(level) % time + domain(level) % dtglobal
      endif
    endif

    !$omp parallel
    call Copy_State ( domain(level) )

    rk: do i_rk = 1, obj_time_scheme%n_RK
        
      call Fill_Ghost_Cell ( domain(level) )               ! Fill ghost cells
      call Fluxes ( domain(level) )                        ! Diffusive fluxes 
      call BC_Fluxes ( domain(level) )                     ! Boundary fluxes

      call External_Function ( domain(level) )             ! External function (e.g. source terms)

      call RK_Newstate ( domain(level), i_rk )             ! State update

    enddo rk
    !$omp end parallel

    ! Save residuals
    if ( (obj_sim_param%iter_from_call == 1) .or. (mod (domain(level) % iter, obj_io%res_diter) == 0) ) then
      obj_sim_param%residuotot = 0.d0
      do b = 1, domain(level) % nb
        if (.not. is_local_block(b)) cycle
        call Compute_Residual ( new=domain(level)%blk(b)%T, &
                                old=domain(level)%blk(b)%TO, &
                                dt=domain(level)%blk(b)%dtlocal, &
                                n=domain(level)%blk(b)%dim, &
                                average=average, &
                                total=obj_sim_param%residuotot )
      enddo
      call mpi_reduce_sum_r8(obj_sim_param%residuotot)
      if (mpi_is_root) obj_sim_param%residuotot = sqrt ( obj_sim_param%residuotot ) ! L2 norm time derivative
    endif
    
    ! Determine simulation control flags on root, then broadcast to all ranks
    if (mpi_is_root) then
      if (level == 1) then
        endsim = ( obj_sim_param%iter_from_call >= domain(1) % itermax ) &
            .or. ( obj_sim_param%residuotot <= obj_sim_param%res_threshold ) &
            .or. ( domain(1) % time >= obj_sim_param%time_threshold )
      else
        endmg = ( obj_sim_param%iter_from_call >= domain(level) % itermax ) &
            .or. ( obj_sim_param%residuotot <= obj_sim_param%res_threshold ) &
            .or. ( domain(level) % time >= obj_sim_param%time_threshold )
      endif

      iosim  = ( mod (domain(level) % iter, obj_io%sol_diter) == 0) &
            .or. ( domain(level) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime )

      if ( endsim ) then
        obj_sim_param%TODO = 3
      elseif ( iosim .or. endmg ) then
        obj_sim_param%TODO = 2
        if (endmg .and. obj_multigrid%MGL > 1) then
          obj_multigrid%change_MG = .true.
        endif
      else
        obj_sim_param%TODO = 1
      endif
    end if

    ! Broadcast simulation control from root to all ranks
    call mpi_bcast_integer(obj_sim_param%TODO)
    call mpi_bcast_logical(obj_multigrid%change_MG)

  end subroutine Explicit_Step


  subroutine Copy_State ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, j, k, b

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        
        domain % blk(b) % TO (i,j,k) = domain % blk(b) % T (i,j,k)
      
      enddo; enddo; enddo
      !$omp end do
    enddo

  end subroutine Copy_State

end module FUSS_Mod_Explicit