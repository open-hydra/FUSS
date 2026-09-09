!>@brief Block-level test of the MORPH update loop.
!>
!> test_morph_gcl.f90 proves the GCL identity for a single cell. This program
!> proves it survives the full per-step path -- Morph_Init, then repeated
!> Morph_Update with a real motion law over a real multi-cell block -- which is
!> where bookkeeping errors (stale time levels, wrong roll order, allocation
!> reuse) actually appear.
!>
!> It also pins the two properties the solver will depend on:
!>   * the STATIC law must leave the geometry exactly unchanged and produce
!>     exactly zero swept volumes, which is what makes the zero-motion
!>     regression (plan 05 section 5.2) meaningful;
!>   * the PRESCRIBED law must keep the GCL residual at round-off for hundreds
!>     of steps of non-uniform, non-axis-aligned motion.
program test_morph_block
  use iso_fortran_env, only: R8 => real64
  use Morph_API

  implicit none

  integer,  parameter :: NI = 12, NJ = 10, NK = 8
  integer,  parameter :: NSTEP = 200
  real(R8), parameter :: DT = 1.0e-3_R8
  real(R8), parameter :: TOL_REL = 1.0e-12_R8

  type(morph_vec3_t), allocatable :: node(:,:,:), node_old(:,:,:), node_ref(:,:,:)
  type(morph_geom_t)   :: geom
  type(morph_status_t) :: st
  type(morph_motion_static_t)     :: law_static
  type(morph_motion_prescribed_t) :: law_moving

  integer  :: dim(3), s
  real(R8) :: t, r_abs, r_rel, worst_rel
  real(R8) :: min_vol, min_jac, vol_sum0, vol_sum
  logical  :: ok

  dim = [NI, NJ, NK]
  ok = .true.

  allocate( node    (0:NI, 0:NJ, 0:NK) )
  allocate( node_old(0:NI, 0:NJ, 0:NK) )
  allocate( node_ref(0:NI, 0:NJ, 0:NK) )

  call build_mesh( node )
  node_ref = node

  write(*,'(A)')          '== Morph block-level update test =='
  write(*,'(A,3(I0,A))')  '   block: ', NI, ' x ', NJ, ' x ', NK, ' cells'
  write(*,'(A,I0)')       '   steps: ', NSTEP
  write(*,*)

  ! =========================================================================
  ! 1. STATIC law: geometry must not change at all.
  ! =========================================================================
  call Morph_Init( node, dim, geom, st )
  if ( .not. st%ok() ) then
    write(*,'(A,A)') '   FAIL Morph_Init: ', trim(st%message); ok = .false.
  endif
  vol_sum0 = sum( geom%vol(1:NI,1:NJ,1:NK) )

  t = 0.0_R8
  do s = 1, 20
    call Morph_Update( node, node_old, dim, t, DT, law_static, geom, st )
    if ( .not. st%ok() ) then
      write(*,'(A,A)') '   FAIL static update: ', trim(st%message); ok = .false.; exit
    endif
    t = t + DT
  enddo

  if ( max_node_diff( node, node_ref ) > 0.0_R8 ) then
    write(*,'(A,ES10.2)') '   FAIL: static law moved the mesh, max |dx| = ', &
                          max_node_diff( node, node_ref )
    ok = .false.
  else
    write(*,'(A)') '   static law leaves nodes exactly unchanged: PASS'
  endif

  if ( any( geom%dV_swept /= 0.0_R8 ) ) then
    write(*,'(A,ES10.2)') '   FAIL: static law gave nonzero swept volume, max ', &
                          maxval(abs(geom%dV_swept))
    ok = .false.
  else
    write(*,'(A)') '   static law gives exactly zero swept volumes: PASS'
  endif

  vol_sum = sum( geom%vol(1:NI,1:NJ,1:NK) )
  if ( vol_sum /= vol_sum0 ) then
    write(*,'(A,ES23.15,A,ES23.15)') '   FAIL: static volume drifted ', vol_sum0, ' -> ', vol_sum
    ok = .false.
  else
    write(*,'(A,ES23.15,A)') '   static total volume exactly constant (', vol_sum, '): PASS'
  endif

  ! =========================================================================
  ! 2. PRESCRIBED law: GCL must hold every step.
  ! =========================================================================
  node = node_ref
  call Morph_Init( node, dim, geom, st )

  ! Non-uniform, three-dimensional, not aligned with any computational
  ! direction, and with the three components out of phase.
  law_moving%name  = 'prescribed-test'
  law_moving%amp   = [ 0.012_R8, 0.009_R8, 0.007_R8 ]
  law_moving%kx    = [ 3.1_R8,  1.7_R8,  2.3_R8 ]
  law_moving%ky    = [ 2.2_R8,  3.5_R8,  1.1_R8 ]
  law_moving%kz    = [ 1.3_R8,  2.9_R8,  3.7_R8 ]
  law_moving%omega = 7.0_R8

  worst_rel = 0.0_R8
  t = 0.0_R8
  do s = 1, NSTEP
    call Morph_Update( node, node_old, dim, t, DT, law_moving, geom, st )
    if ( .not. st%ok() ) then
      write(*,'(A,I0,A,A,A,3I4)') '   FAIL moving update at step ', s, ': ', &
            trim(st%message), ' at cell', st%i, st%j, st%k
      ok = .false.; exit
    endif

    call Morph_GCL_Residual( geom, dim, r_abs, r_rel )
    worst_rel = max( worst_rel, r_rel )
    t = t + DT
  enddo

  write(*,'(A,ES10.2)') '   worst relative GCL residual over all steps: ', worst_rel
  if ( worst_rel > TOL_REL ) then
    write(*,'(A)') '   FAIL: GCL residual above tolerance'; ok = .false.
  else
    write(*,'(A)') '   GCL holds every step under prescribed motion: PASS'
  endif

  ! =========================================================================
  ! 3. FACE PAIRING -- the property the ALE flux's conservation rests on.
  !
  ! Cell (i,j,k)'s face 2 (i-high) and cell (i+1,j,k)'s face 1 (i-low) are the
  ! SAME physical face. For the solver's swept-volume fluxes to telescope, the
  ! two cells must see equal and opposite swept volumes. If they do not, energy
  ! is created or destroyed at every interior face and no amount of care in the
  ! state update can recover it.
  ! =========================================================================
  block
    real(R8) :: pair_max, denom
    integer  :: ii, jj, kk
    pair_max = 0.0_R8
    do kk = 1, NK
    do jj = 1, NJ
    do ii = 1, NI
      denom = max( abs(geom%vol(ii,jj,kk)), 1.0e-300_R8 )
      if ( ii < NI ) pair_max = max( pair_max, &
           abs( geom%dV_swept(2,ii,jj,kk) + geom%dV_swept(1,ii+1,jj,kk) ) / denom )
      if ( jj < NJ ) pair_max = max( pair_max, &
           abs( geom%dV_swept(4,ii,jj,kk) + geom%dV_swept(3,ii,jj+1,kk) ) / denom )
      if ( kk < NK ) pair_max = max( pair_max, &
           abs( geom%dV_swept(6,ii,jj,kk) + geom%dV_swept(5,ii,jj,kk+1) ) / denom )
    enddo; enddo; enddo
    write(*,'(A,ES10.2)') '   worst face-pairing asymmetry (relative): ', pair_max
    if ( pair_max > TOL_REL ) then
      write(*,'(A)') '   FAIL: shared faces disagree on their swept volume;'
      write(*,'(A)') '         the ALE flux cannot be conservative.'
      ok = .false.
    else
      write(*,'(A)') '   shared faces see equal and opposite swept volumes: PASS'
    endif
  end block

  call Morph_Quality_Check( node, dim, min_vol, min_jac, st )
  write(*,'(A,ES12.4,A,ES12.4)') '   final mesh: min signed vol = ', min_vol, &
                                 '   min corner Jacobian = ', min_jac
  if ( .not. st%ok() ) then
    write(*,'(A,A)') '   FAIL quality: ', trim(st%message); ok = .false.
  else
    write(*,'(A)') '   mesh remained untangled: PASS'
  endif

  write(*,*)
  if ( ok ) then
    write(*,'(A)') '   RESULT: PASS'
  else
    write(*,'(A)') '   RESULT: FAIL'
    error stop 1
  endif

contains

  !> Largest absolute node displacement between two meshes. Exact zero is the
  !> pass criterion for the static law, so no tolerance is applied here.
  real(R8) function max_node_diff ( a, b ) result ( d )
    type(morph_vec3_t), intent(in) :: a(0:,0:,0:), b(0:,0:,0:)
    integer :: i, j, k

    d = 0.0_R8
    do k = lbound(a,3), ubound(a,3)
    do j = lbound(a,2), ubound(a,2)
    do i = lbound(a,1), ubound(a,1)
      d = max( d, maxval( abs( a(i,j,k)%c - b(i,j,k)%c ) ) )
    enddo; enddo; enddo

  end function max_node_diff


  !> A deliberately non-uniform, non-orthogonal block: geometric stretching in
  !> every direction plus a shear, so faces are genuinely warped. A uniform
  !> Cartesian box would let several classes of metric bug pass unnoticed.
  subroutine build_mesh ( nd )
    type(morph_vec3_t), intent(out) :: nd(0:,0:,0:)
    integer  :: i, j, k
    real(R8) :: xi, et, ze

    do k = 0, NK
    do j = 0, NJ
    do i = 0, NI
      xi = real(i,R8) / real(NI,R8)
      et = real(j,R8) / real(NJ,R8)
      ze = real(k,R8) / real(NK,R8)

      ! stretched, then sheared
      nd(i,j,k)%c(1) = 0.40_R8 * xi**1.20_R8 + 0.05_R8 * et * ze
      nd(i,j,k)%c(2) = 0.30_R8 * et**0.85_R8 + 0.04_R8 * xi * ze
      nd(i,j,k)%c(3) = 0.25_R8 * ze**1.10_R8 + 0.03_R8 * xi * et
    enddo; enddo; enddo

  end subroutine build_mesh

end program test_morph_block
