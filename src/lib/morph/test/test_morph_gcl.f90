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
!>       src/lib/morph/metrics/Morph_GCL.f90 \
!>       src/lib/morph/metrics/Morph_Metrics.f90 \
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

    ! ---- How much did fixing the volume formula actually change?
    ! Morph_Metric_Tensor now returns the signed volume, so comparing it against
    ! itself would prove nothing. legacy_volume() below is a local copy of the
    ! five-tetrahedron abs() sum that FUSS used before, kept HERE (not in the
    ! component) purely so this test can quantify the correction.
    vol_legacy = legacy_volume( x_old )
    gap = Morph_GCL_Volume_Gap( x_old, vol_legacy )
    worst_gap = max( worst_gap, gap )

  enddo

  write(*,'(A,ES10.2)') '   worst relative GCL residual      : ', worst_rel
  write(*,'(A,ES10.2)') '   worst signed-vs-legacy volume gap: ', worst_gap
  write(*,'(A)')        '     (that gap is the size of the correction from fixing the'
  write(*,'(A)')        '      old five-tetrahedron abs() volume; it is 0 for planar faces)'
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

  !> FUSS's ORIGINAL cell volume: five tetrahedra, absolute value of each.
  !> Kept here only so the test can report how large the correction was. Do not
  !> use this anywhere else -- it is exact for planar-faced cells, carries an
  !> O(face warp) error otherwise, and cannot detect an inverted cell.
  real(R8) function legacy_volume ( x ) result ( vol )
    real(R8), intent(in) :: x(3,8)
    real(R8) :: vx(8), vy(8), vz(8)

    ! The original remapped the corner order before summing.
    vx(1)=x(1,1); vy(1)=x(2,1); vz(1)=x(3,1)
    vx(2)=x(1,5); vy(2)=x(2,5); vz(2)=x(3,5)
    vx(3)=x(1,3); vy(3)=x(2,3); vz(3)=x(3,3)
    vx(4)=x(1,7); vy(4)=x(2,7); vz(4)=x(3,7)
    vx(5)=x(1,2); vy(5)=x(2,2); vz(5)=x(3,2)
    vx(6)=x(1,6); vy(6)=x(2,6); vz(6)=x(3,6)
    vx(7)=x(1,4); vy(7)=x(2,4); vz(7)=x(3,4)
    vx(8)=x(1,8); vy(8)=x(2,8); vz(8)=x(3,8)

    ! tvol is a sibling internal procedure, not nested inside this one:
    ! Fortran does not allow an internal procedure to contain further internal
    ! procedures.
    vol = tvol(vx,vy,vz,1,2,3,5) + tvol(vx,vy,vz,2,4,3,8) &
        + tvol(vx,vy,vz,5,8,6,2) + tvol(vx,vy,vz,5,7,8,3) &
        + tvol(vx,vy,vz,5,8,2,3)

  end function legacy_volume


  real(R8) function tvol ( vx, vy, vz, ia, ib, ic, id ) result ( v )
    real(R8), intent(in) :: vx(8), vy(8), vz(8)
    integer,  intent(in) :: ia, ib, ic, id

    v = abs(((vx(ib)-vx(ia))* &
      ((vy(ic)-vy(ia))*(vz(id)-vz(ia))-(vy(id)-vy(ia))*(vz(ic)-vz(ia)))+ &
                              (vy(ib)-vy(ia))* &
      ((vx(id)-vx(ia))*(vz(ic)-vz(ia))-(vx(ic)-vx(ia))*(vz(id)-vz(ia)))+ &
                              (vz(ib)-vz(ia))* &
      ((vx(ic)-vx(ia))*(vy(id)-vy(ia))-(vx(id)-vx(ia))*(vy(ic)-vy(ia)))) &
      /6.d0)

  end function tvol


  !> Deterministic pseudo-random sequence in (0,1). Reproducible across runs and
  !> platforms, so a failure is always reproducible.
  real(R8) function rnd ( s )
    real(R8), intent(inout) :: s
    s = mod( 8121.0_R8 * s + 28411.0_R8, 134456.0_R8 )
    rnd = s / 134456.0_R8
  end function rnd

end program test_morph_gcl
