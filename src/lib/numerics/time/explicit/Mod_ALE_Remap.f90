!>@brief Conservative remap of the state onto the moved mesh (ALE).
!>
!> WHY A REMAP RATHER THAN AN EXTRA TERM IN THE RK STAGES
!> -----------------------------------------------------
!> The mesh moves ONCE per time step, but a Runge-Kutta step evaluates its
!> residual once per STAGE. Folding the swept-volume flux and the cell-volume
!> change into the stage residual therefore applies a per-step quantity three
!> times, each with a different RK weight and a different intermediate state.
!>
!> That was measured, not guessed: with the terms inside the RK stages, global
!> energy was off by 5.7e-4 relative and -- the diagnostic that settles it --
!> the error did NOT shrink when dt was halved or quartered at fixed end time
!> (5.678e-4, 5.719e-4, 5.740e-4). A temporal-consistency artifact would have
!> converged; a constant error means the accounting itself is wrong.
!>
!> Splitting the two operations fixes it exactly:
!>
!>   1. REMAP (here). Transfer the state from the old mesh to the new one:
!>
!>        h_new = [ h_old * V_old  +  sum_f h_upwind * dV_swept(f) ] / V_new
!>
!>      applied ONCE per step, outside the RK stages.
!>
!>   2. CONDUCT. The ordinary static-mesh Runge-Kutta step, on the new mesh,
!>      completely unchanged. This is why Lib_Newstate needs no ALE branch.
!>
!> PROPERTIES
!>
!>   Exactly conservative. Summing over cells, sum(V_new h_new) = sum(V_old
!>   h_old) + sum_cells sum_f h_upwind dV. Two cells sharing a face see equal
!>   and opposite dV (verified in test_morph_block to 2.4e-15) and, because the
!>   upwind value is selected on the SIGN of dV, they select the SAME h. The
!>   double sum therefore telescopes to zero over interior faces.
!>
!>   Free-stream preserving. For uniform h the bracket is h*(V_old + (V_new -
!>   V_old)) = h*V_new by the discrete GCL, so h_new = h exactly.
!>
!>   Inert when the mesh does not move: dV == 0 and V_new == V_old give
!>   h_new = h_old, bit-for-bit.
module FUSS_Mod_ALE_Remap
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Remap_State_ALE

contains

  !> Conservative transfer of the state onto the moved mesh.
  !>
  !> Call ONCE per time step, after the mesh has moved and before Copy_State,
  !> from inside an OpenMP parallel region (it uses `!$omp do`).
  !>
  !> Ghost temperatures must already be valid: the upwind value at a block
  !> interface is read from the ghost cell.
  subroutine Remap_State_ALE ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_mesh_motion
    use FUSS_Lib_Solid,      only: co_H, co_T
    use FUSS_Mod_MPI,        only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: i, j, k, b, f, ii, jj, kk
    real(R8) :: h, h_up, hsum, dV, vnew, vold

    if ( .not. obj_mesh_motion%enabled ) return

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle

      ! PASS 1 -- compute the remapped volumetric enthalpy for every cell,
      ! parking it in TO. Reading T and writing TO keeps every cell's upwind
      ! lookup on the SAME pre-remap state; writing T in place here would let a
      ! cell see its neighbour's already-remapped value, the two sides of a face
      ! would disagree on the upwind enthalpy, and conservation would break.
      ! (TO is free: Copy_State overwrites it immediately after this routine.)
      !$omp do collapse(3) schedule(static) private(i,j,k,f,ii,jj,kk,h,h_up,hsum,dV,vnew,vold)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)

        vnew = domain%blk(b)%vol(i,j,k)
        vold = domain%blk(b)%vol_old(i,j,k)

        call co_H( domain%blk(b)%matID(i,j,k), domain%blk(b)%T(i,j,k), h )
        hsum = h * vold

        do f = 1, 6
          dV = domain%blk(b)%dV_swept(f,i,j,k)
          if ( dV == 0.0d0 ) cycle

          if ( dV > 0.0d0 ) then
            ! Face swept outward: the volume gained used to belong to the
            ! neighbour, so it arrives carrying the neighbour's enthalpy.
            ii = i; jj = j; kk = k
            select case (f)
              case (1); ii = i - 1
              case (2); ii = i + 1
              case (3); jj = j - 1
              case (4); jj = j + 1
              case (5); kk = k - 1
              case (6); kk = k + 1
            end select
            call co_H( domain%blk(b)%matID(ii,jj,kk), domain%blk(b)%T(ii,jj,kk), h_up )
          else
            ! Face swept inward: this cell loses volume, carrying its own.
            h_up = h
          endif

          hsum = hsum + h_up * dV
        enddo

        domain%blk(b)%TO(i,j,k) = hsum / vnew

      enddo; enddo; enddo
      !$omp end do

      ! PASS 2 -- convert back to temperature.
      !$omp do collapse(3) schedule(static) private(i,j,k)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        call co_T( domain%blk(b)%matID(i,j,k), domain%blk(b)%TO(i,j,k), &
                   domain%blk(b)%T(i,j,k) )
      enddo; enddo; enddo
      !$omp end do

    enddo

  end subroutine Remap_State_ALE

end module FUSS_Mod_ALE_Remap
