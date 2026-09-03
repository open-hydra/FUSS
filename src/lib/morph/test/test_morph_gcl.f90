!>@brief Unit test for the discrete Geometric Conservation Law.
!>
!> Proves, without involving the solver, that
!>     V^{n+1} - V^n  ==  sum_f dV_swept(f)
!> to floating-point round-off for arbitrary node motion of a distorted cell.
!>
!> This is deliberately a standalone program: the GCL identity is the single
!> most correctness-critical thing in the moving-mesh phase, and it is far
!> cheaper and sharper to test here than through a full solver run.
!>
!> Build (no CMake needed):
!>   gfortran -O2 -fcheck=bounds -o test_morph_gcl \
!>       src/lib/morph/base/Morph_Types_m.f90 \
!>       src/lib/morph/metrics/Morph_Metrics.f90 \
!>       src/lib/morph/metrics/Morph_GCL.f90 \
!>       src/lib/morph/test/test_morph_gcl.f90
program test_morph_gcl
  use iso_fortran_env, only: R8 => real64
  use Morph_Types_m
  use Morph_Metrics, only: Morph_Metric_Tensor
  use Morph_GCL

  implicit none

  integer, parameter :: NTRIAL = 200
  real(R8), parameter :: TOL_REL = 1.0e-13_R8   ! round-off budget for ~1e2 flops

  real(R8) :: x_old(3,8), x_new(3,8)
  real(R8) :: v_old, v_new, dV(6), resid, rel, worst_rel
  real(R8) :: seed
  integer  :: t, n, c, nfail
  logical  :: ok

  ! Legacy-formula comparison
  type(morph_vec3_t)   :: Nd(8)
  type(morph_tens3_t)  :: M
  type(morph_vec3_t)   :: dl
  real(R8)             :: vol_legacy, gap, worst_gap
  type(morph_status_t) :: st

  nfail = 0
  worst_rel = 0.0_R8
  worst_gap = 0.0_R8
  seed = 0.123456789_R8

  write(*,'(A)') '== Morph GCL unit test =='
  write(*,'(A,I0,A)') '   trials: ', NTRIAL, ' distorted cells with arbitrary motion'
  write(*,'(A,ES10.2)') '   relative tolerance: ', TOL_REL
  write(*,*)

  do t = 1, NTRIAL

    ! ---- A unit cell, then displaced corners so the cell is genuinely warped
    !      (non-planar faces), not a parallelepiped.
    call unit_cell( x_old )
    do n = 1, 8
      do c = 1, 3
        x_old(c,n) = x_old(c,n) + 0.30_R8 * ( rnd(seed) - 0.5_R8 )
      enddo
    enddo

    ! ---- Arbitrary motion: every node moves independently, in all three
    !      directions, by a large fraction of the cell size. This is much more
    !      demanding than the smooth stretching a real recession law produces.
    do n = 1, 8
      do c = 1, 3
        x_new(c,n) = x_old(c,n) + 0.25_R8 * ( rnd(seed) - 0.5_R8 )
      enddo
    enddo

    v_old = Morph_Cell_Volume_Signed( x_old )
    v_new = Morph_Cell_Volume_Signed( x_new )
    call Morph_Swept_Volumes_Cell( x_old, x_new, dV )

    resid = abs( ( v_new - v_old ) - sum(dV) )
    rel   = resid / max( abs(v_new), abs(v_old) )
    worst_rel = max( worst_rel, rel )
    if ( rel > TOL_REL ) then
      nfail = nfail + 1
      if ( nfail <= 3 ) then
        write(*,'(A,I4)')       '   FAIL trial ', t
        write(*,'(A,ES23.15)')  '     V_old      = ', v_old
        write(*,'(A,ES23.15)')  '     V_new      = ', v_new
        write(*,'(A,ES23.15)')  '     dV_new-old = ', v_new - v_old
        write(*,'(A,ES23.15)')  '     sum(dV)    = ', sum(dV)
        write(*,'(A,ES23.15)')  '     rel resid  = ', rel
      endif
    endif

    ! ---- Signed vs legacy absolute-value volume, on the untangled old cell.
    do n = 1, 8
      Nd(n)%c = x_old(:,n)
    enddo
    call st%clear()
    call Morph_Metric_Tensor( Nd(1),Nd(2),Nd(3),Nd(4),Nd(5),Nd(6),Nd(7),Nd(8), &
                              M, dl, vol_legacy, st )
    gap = Morph_GCL_Volume_Gap( x_old, vol_legacy )
    worst_gap = max( worst_gap, gap )

  enddo

  write(*,'(A,ES10.2)') '   worst relative GCL residual      : ', worst_rel
  write(*,'(A,ES10.2)') '   worst signed-vs-legacy volume gap: ', worst_gap
  write(*,*)

  ok = ( nfail == 0 )
  if ( ok ) then
    write(*,'(A)') '   RESULT: PASS -- discrete GCL holds to round-off'
  else
    write(*,'(A,I0,A,I0,A)') '   RESULT: FAIL -- ', nfail, ' of ', NTRIAL, ' trials exceeded tolerance'
  endif

  ! ---- Degenerate check: zero motion must give exactly zero swept volume.
  call unit_cell( x_old )
  x_new = x_old
  call Morph_Swept_Volumes_Cell( x_old, x_new, dV )
  if ( any( dV /= 0.0_R8 ) ) then
    write(*,'(A,ES10.2)') '   RESULT: FAIL -- zero motion gave nonzero swept volume, max ', maxval(abs(dV))
    ok = .false.
  else
    write(*,'(A)') '   zero-motion swept volumes are exactly zero: PASS'
  endif

  ! ---- Uniform translation must give zero net volume change.
  call unit_cell( x_old )
  x_new = x_old
  x_new(1,:) = x_new(1,:) + 0.37_R8
  x_new(2,:) = x_new(2,:) - 0.11_R8
  call Morph_Swept_Volumes_Cell( x_old, x_new, dV )
  if ( abs(sum(dV)) > TOL_REL ) then
    write(*,'(A,ES10.2)') '   RESULT: FAIL -- rigid translation changed volume by ', sum(dV)
    ok = .false.
  else
    write(*,'(A)') '   rigid translation conserves volume: PASS'
  endif

  ! ---- Sign convention: a uniformly expanding cell must give POSITIVE volume
  !      and positive swept volumes on every face. This is what catches an
  !      inward/outward orientation mistake, which the GCL residual cannot see
  !      (it is satisfied by a globally wrong sign just as well).
  call unit_cell( x_old )
  v_old = Morph_Cell_Volume_Signed( x_old )
  if ( v_old <= 0.0_R8 ) then
    write(*,'(A,F12.8)') '   RESULT: FAIL -- unit cube has non-positive volume ', v_old
    ok = .false.
  else
    write(*,'(A,F12.8,A)') '   unit cube volume = ', v_old, ' (expected +1): PASS'
  endif

  x_new = 1.5_R8 * x_old            ! uniform dilation about the origin
  call Morph_Swept_Volumes_Cell( x_old, x_new, dV )
  if ( any( dV < 0.0_R8 ) .or. sum(dV) <= 0.0_R8 ) then
    write(*,'(A)')        '   RESULT: FAIL -- expansion produced a negative swept volume'
    write(*,'(A,6F10.5)') '     dV = ', dV
    ok = .false.
  else
    write(*,'(A,F12.8,A)') '   expansion: all swept volumes positive, sum = ', sum(dV), ': PASS'
  endif

  if ( .not. ok ) error stop 1

contains

  subroutine unit_cell ( x )
    real(R8), intent(out) :: x(3,8)
    ! Node ordering: 1=(0,0,0) 2=(0,0,1) 3=(0,1,0) 4=(0,1,1)
    !                5=(1,0,0) 6=(1,0,1) 7=(1,1,0) 8=(1,1,1)
    x(:,1) = [0.0_R8, 0.0_R8, 0.0_R8]
    x(:,2) = [0.0_R8, 0.0_R8, 1.0_R8]
    x(:,3) = [0.0_R8, 1.0_R8, 0.0_R8]
    x(:,4) = [0.0_R8, 1.0_R8, 1.0_R8]
    x(:,5) = [1.0_R8, 0.0_R8, 0.0_R8]
    x(:,6) = [1.0_R8, 0.0_R8, 1.0_R8]
    x(:,7) = [1.0_R8, 1.0_R8, 0.0_R8]
    x(:,8) = [1.0_R8, 1.0_R8, 1.0_R8]
  end subroutine unit_cell

  !> Deterministic pseudo-random sequence in (0,1). Reproducible across runs and
  !> platforms, so a failure is always reproducible.
  real(R8) function rnd ( s )
    real(R8), intent(inout) :: s
    s = mod( 8121.0_R8 * s + 28411.0_R8, 134456.0_R8 )
    rnd = s / 134456.0_R8
  end function rnd

end program test_morph_gcl
