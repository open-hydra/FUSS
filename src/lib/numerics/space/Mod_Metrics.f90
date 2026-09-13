module FUSS_Mod_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

contains

  !> Build the geometry for every block and every ghost-cell interface.
  !>
  !> The cell-level geometry (volume, metric tensor, cell lengths, face normals
  !> and areas) is computed by the MORPH component through Adapter_Morph. The
  !> arithmetic is unchanged from the routines that used to live in
  !> Lib_Metrics.f90, so results are bit-identical; see plan 09 for why the
  !> geometry was split out.
  !>
  !> The BC ghost-cell metrics below still use Lib_Metrics because they need the
  !> block-connectivity topology. Plan 09 section 7 recommends moving them into
  !> MORPH too, passing connectivity as plain integer descriptors; that is a
  !> separate step and is deliberately not done here, so that this change can be
  !> verified as inert on its own.
  subroutine Setup_Metrics ( domain )
    use FUSS_Base_Types_m
    use FUSS_Advanced_Types_m
    use FUSS_Global_m
    use FUSS_Lib_Metrics, only: delthe
    use FUSS_Adapter_Morph, only: Adapter_Mesh_Type, Adapter_Block_Metrics
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: b
    integer :: ci, cj, ck
    logical :: ok
    character(len=256) :: message

    ! Mesh classification. These were globals written directly by
    ! Check_Mesh_Type; MORPH returns them instead of reaching into solver state.
    call Adapter_Mesh_Type ( domain%blk(1), ndir, delthe )

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      call Adapter_Block_Metrics ( domain%blk(b), domain%nb, b, ok, message, ci, cj, ck )

      if ( .not. ok ) then
        write(*,'(A)')            ' [ERROR] mesh geometry is invalid'
        write(*,'(A,A)')          '         ', trim(message)
        write(*,'(A,4(I0,A))')    '         at block ', b, ', cell (', ci, ',', cj, ',', ck, ')'
        error stop 2
      endif
    enddo

    call Setup_BC_Metrics ( domain )

  end subroutine Setup_Metrics


  !> Build the ghost-cell metrics for every boundary entry.
  !>
  !> Split out of Setup_Metrics so that Update_Mesh can redo it every step: these
  !> are functions of the block geometry, and on a moving mesh leaving them at
  !> their t = 0 values means the geometry the solver describes at the boundary
  !> is not the geometry it is solving on.
  !>
  !> bc%Mg IS LOAD-BEARING. It is read in exactly one place, BC_Connection (bc
  !> types 101/102), where it supplies the NEIGHBOUR's metric in the
  !> conductivity-and-metric weighted interface temperature
  !>
  !>     T_int = (k1 T1 D1 + k2 T2 D2) / (k1 D1 + k2 D2)
  !>
  !> which is exactly the flux-continuity interface temperature between two cell
  !> centres at dx1/2 and dx2/2. Get D2 wrong and the interface temperature is
  !> wrong.
  !>
  !> Measured, on test/numerics-features/conjugate_two_layer -- two materials,
  !> dx2/dx1 = 4, exact solution piecewise linear:
  !>
  !>     correct                 max |T - T_exact| = 4.5e-09 K
  !>     every bc%Mg halved      max |T - T_exact| = 2.9e+01 K
  !>
  !> Be careful how you test this. The same sabotage on multimat_Plate changes
  !> NOTHING, over 214 steps, bit for bit -- because every interface there has
  !> the same cell spacing on both sides, and with D1 = D2 the weights reduce to
  !> the conductivity ratio and the neighbour's metric drops out. That is the one
  !> configuration in which this field looks dead, and until 2026-09-10 it was
  !> the only multi-block configuration the suite contained.
  !>
  !> Hence the per-step refresh: as soon as mesh motion changes the spacing near
  !> a block interface, a stale bc%Mg is a wrong interface temperature. Plan 08's
  !> receding surface compresses exactly those cells.
  !>
  !> bc%dlg and bc%volg are a different matter -- they are written here and read
  !> NOWHERE in the solver. They are left in place deliberately rather than
  !> commented out: they are `intent(out)` dummies of BC_Connect_Metrics and
  !> friends, so suppressing the assignments would leave them undefined, which is
  !> the same uninitialised-read defect class already fixed three times in this
  !> codebase. If they are ever removed, remove the arguments too.
  !>
  !> This does not affect a static mesh: Update_Mesh returns before calling it
  !> when motion is disabled.
  subroutine Setup_BC_Metrics ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Lib_Metrics, only: BC_Connect_Metrics, BC_Symmetry_Metrics, &
                                BC_Extrapolate_Metrics
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: i
    integer :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s

    !$omp parallel
    ! Create the nodes for gc layers of ghost cell
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, domain % nbound
      Bm = domain % bc(i) % b
      Im = domain % bc(i) % i
      Jm = domain % bc(i) % j
      Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
      select case ( domain % bc(i) % type )
        case(101) ! block connection
          Bs = domain % bc(i) % bs
          Is = domain % bc(i) % is
          Js = domain % bc(i) % js
          Ks = domain % bc(i) % ks
          Fs = domain % bc(i) % fs
          d11s = domain % bc(i) % d11
          d12s = domain % bc(i) % d12
          d21s = domain % bc(i) % d21
          d22s = domain % bc(i) % d22
          call BC_Connect_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                    Is, Js, Ks, Fs, domain % blk(Bs), d11s, d12s, d21s, d22s, &
                                    domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg)
        case(300) ! symmetry
          call BC_Symmetry_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                     domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg )
        case default  ! Rientra qui anche chimera
          call BC_Extrapolate_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                        domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg )
      end select
    enddo
    !$omp end parallel

  end subroutine Setup_BC_Metrics


  !> Advance the mesh one step and refresh all geometry. Called once per time
  !> step, NOT once per Runge-Kutta stage: the swept volumes describe the motion
  !> over the whole step, and re-applying the law per stage would move the mesh
  !> n_RK times as far.
  !>
  !> Must be called OUTSIDE any OpenMP parallel region -- the MORPH kernels open
  !> their own, exactly as Setup_Metrics does.
  !>
  !> `dt` here only decides WHERE the mesh is placed. The discrete GCL holds for
  !> whatever motion actually occurred, because the swept volumes are computed
  !> from the old and new node positions themselves, not from dt. So an imperfect
  !> dt estimate costs accuracy, never conservation.
  subroutine Update_Mesh ( domain, t, dt )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_mesh_motion
    use FUSS_Adapter_Morph,  only: Adapter_Update_Mesh
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    real(R8),               intent(in)    :: t, dt
    ! Local
    integer  :: b, ci, cj, ck, si, sj, sk
    logical  :: ok
    real(R8) :: gcl_rel, sweep_ratio
    character(len=256) :: message

    if ( .not. obj_mesh_motion%enabled ) return

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle

      call Adapter_Update_Mesh ( domain%blk(b), domain%nb, b, t, dt, &
                                 gcl_rel, sweep_ratio, si, sj, sk, &
                                 ok, message, ci, cj, ck )

      if ( .not. ok ) then
        write(*,'(A)')         ' [ERROR] mesh update failed'
        write(*,'(A,A)')       '         ', trim(message)
        write(*,'(A,4(I0,A))') '         at block ', b, ', cell (', ci, ',', cj, ',', ck, ')'
        error stop 2
      endif

      ! Assert the discrete GCL every step. A violation means the geometry and
      ! the ALE state update disagree about how much volume was swept, which
      ! corrupts the solution silently -- so it is fatal, not a warning.
      if ( gcl_rel > obj_mesh_motion%gcl_tol ) then
        write(*,'(A)')          ' [ERROR] discrete Geometric Conservation Law violated'
        write(*,'(A,ES12.4,A,ES12.4)') '         relative residual ', gcl_rel, &
                                       ' exceeds tolerance ', obj_mesh_motion%gcl_tol
        write(*,'(A,I0)')       '         block ', b
        error stop 2
      endif

      ! Moving-mesh stability. The GCL says the geometry is self-consistent; it
      ! says nothing about the mesh having moved too far in one step. Once a
      ! cell sweeps out more than its own volume the remap's coefficient on that
      ! cell's own state goes negative and the update extrapolates. dt is chosen
      ! from conduction alone and carries no mesh-velocity limit, so nothing
      ! upstream prevents this.
      if ( sweep_ratio >= 1.0d0 ) then
        write(*,'(A)')          ' [ERROR] the mesh moved too far in one time step'
        write(*,'(A,ES12.4,A)') '         a cell swept out ', sweep_ratio, ' times its own volume;'
        write(*,'(A)')          '         the conservative remap is unstable at or above 1.'
        write(*,'(A,4(I0,A))')  '         worst at block ', b, ', cell (', si, ',', sj, ',', sk, ')'
        write(*,'(A)')          '         Reduce the mesh velocity or the time step (lower vnn).'
        error stop 2
      endif
    enddo

    ! The ghost-cell metrics are functions of the geometry that has just moved,
    ! so they have to be rebuilt too. See Setup_BC_Metrics.
    call Setup_BC_Metrics ( domain )

  end subroutine Update_Mesh

end module FUSS_Mod_Metrics
