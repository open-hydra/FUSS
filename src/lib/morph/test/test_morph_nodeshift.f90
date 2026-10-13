!>@brief Unit test of the nodeshift law (plan 08 section 2.3, stage B).
!>
!> One graded 10 x 4 x 1 block. Per receding face (i-high, i-low, j-high), with a
!> displacement that varies across the face, over many steps:
!>   * the discrete GCL residual stays at round-off (<= 1e-14 relative);
!>   * the receding face sweeps exactly -A*delta in total for a uniform normal
!>     displacement delta (checked on the i-high face of a straight block);
!>   * the nodes of the opposite (pinned) face are bitwise unchanged;
!>   * every cell of a line shrinks by the same fraction (the grading is kept);
!>   * a zero displacement leaves every node bitwise unchanged and every swept
!>     volume exactly zero;
!>   * a displacement equal to the line length returns MORPH_ERR_NEGVOL, and the
!>     law never stops the program.
program test_morph_nodeshift
  use iso_fortran_env, only: R8 => real64
  use Morph_API

  implicit none

  integer,  parameter :: NI = 10, NJ = 4, NK = 1
  real(R8), parameter :: LX = 0.05_R8, LY = 0.02_R8, LZ = 0.01_R8
  real(R8), parameter :: TOL_GCL = 1.0e-14_R8

  type(morph_vec3_t), allocatable :: node(:,:,:), node_old(:,:,:), node0(:,:,:)
  type(morph_geom_t)   :: geom
  type(morph_status_t) :: st
  type(morph_motion_nodeshift_t) :: law
  integer  :: dim(3), s, i, j, k, nfail
  real(R8) :: r_abs, r_rel, worst, sweep, area, delta, frac0, frac, t
  logical  :: pinned_ok, upd_ok

  dim = [NI, NJ, NK]
  nfail = 0
  allocate( node(0:NI,0:NJ,0:NK), node_old(0:NI,0:NJ,0:NK), node0(0:NI,0:NJ,0:NK) )

  write(*,'(A)') '== Morph nodeshift law test =='

  ! ---- 1. i-high face recedes uniformly along -x: sweep, pinned face, grading
  call build_mesh( node )
  node0 = node
  call Morph_Init( node, dim, geom, st )
  call check( st%ok(), 'Morph_Init' )
  law%dir = 1 ; law%side = 2
  allocate( law%disp(0:NJ, 0:NK) )
  delta = 1.0e-4_R8
  do k = 0, NK ; do j = 0, NJ
    law%disp(j,k)%c = [ -delta, 0.0_R8, 0.0_R8 ]
  enddo ; enddo
  worst = 0.0_R8 ; t = 0.0_R8 ; upd_ok = .true.
  do s = 1, 50
    call Morph_Update( node, node_old, dim, t, 1.0_R8, law, geom, st )
    upd_ok = upd_ok .and. st%ok()
    call Morph_GCL_Residual( geom, dim, r_abs, r_rel )
    worst = max( worst, r_rel )
    if ( s == 1 ) then
      sweep = sum( geom%dV_swept(2, NI, 1:NJ, 1:NK) )
      area  = LY * LZ
      call check( abs( sweep + area * delta ) <= 1.0e-12_R8 * area * delta, 'surface sweep = -A delta' )
      write(*,'(A,ES12.4,A,ES12.4)') '   surface sweep ', sweep, '  expected ', -area * delta
    endif
    t = t + 1.0_R8
  enddo
  call check( upd_ok, 'every update ok, i-high' )
  call check( worst <= TOL_GCL, 'GCL residual, i-high' )
  write(*,'(A,ES10.2)') '   max relative GCL residual, i-high: ', worst
  pinned_ok = .true.
  do k = 0, NK ; do j = 0, NJ
    if ( any( node(0,j,k)%c /= node0(0,j,k)%c ) ) pinned_ok = .false.
  enddo ; enddo
  call check( pinned_ok, 'pinned face i = 0 bitwise unchanged' )
  ! uniform compression: every cell on line (j, k) = (1, 0) keeps its share of the line
  frac0 = ( node0(1,1,0)%c(1) - node0(0,1,0)%c(1) ) / ( node0(NI,1,0)%c(1) - node0(0,1,0)%c(1) )
  frac  = ( node(1,1,0)%c(1)  - node(0,1,0)%c(1)  ) / ( node(NI,1,0)%c(1)  - node(0,1,0)%c(1)  )
  call check( abs( frac - frac0 ) <= 1.0e-13_R8, 'grading kept (first cell share of the line)' )
  call check( abs( node(NI,1,0)%c(1) - ( node0(NI,1,0)%c(1) - 50 * delta ) ) <= 1.0e-15_R8, &
              'surface moved by 50 delta' )
  deallocate( law%disp )

  ! ---- 2. i-low face recedes with a displacement varying across the face
  call build_mesh( node )
  node0 = node
  call Morph_Init( node, dim, geom, st )
  law%dir = 1 ; law%side = 1
  allocate( law%disp(0:NJ, 0:NK) )
  do k = 0, NK ; do j = 0, NJ
    law%disp(j,k)%c = [ 1.0e-4_R8 * ( 1.0_R8 + 0.3_R8 * j ), 0.0_R8, 0.0_R8 ]
  enddo ; enddo
  worst = 0.0_R8 ; upd_ok = .true.
  do s = 1, 50
    call Morph_Update( node, node_old, dim, t, 1.0_R8, law, geom, st )
    upd_ok = upd_ok .and. st%ok()
    call Morph_GCL_Residual( geom, dim, r_abs, r_rel )
    worst = max( worst, r_rel )
  enddo
  call check( upd_ok, 'every update ok, i-low' )
  call check( worst <= TOL_GCL, 'GCL residual, i-low, varying displacement' )
  write(*,'(A,ES10.2)') '   max relative GCL residual, i-low:  ', worst
  pinned_ok = .true.
  do k = 0, NK ; do j = 0, NJ
    if ( any( node(NI,j,k)%c /= node0(NI,j,k)%c ) ) pinned_ok = .false.
  enddo ; enddo
  call check( pinned_ok, 'pinned face i = NI bitwise unchanged' )
  deallocate( law%disp )

  ! ---- 3. j-high face recedes (lines along j)
  call build_mesh( node )
  node0 = node
  call Morph_Init( node, dim, geom, st )
  law%dir = 2 ; law%side = 2
  allocate( law%disp(0:NI, 0:NK) )
  do k = 0, NK ; do i = 0, NI
    law%disp(i,k)%c = [ 0.0_R8, -2.0e-4_R8, 0.0_R8 ]
  enddo ; enddo
  worst = 0.0_R8 ; upd_ok = .true.
  do s = 1, 20
    call Morph_Update( node, node_old, dim, t, 1.0_R8, law, geom, st )
    upd_ok = upd_ok .and. st%ok()
    call Morph_GCL_Residual( geom, dim, r_abs, r_rel )
    worst = max( worst, r_rel )
  enddo
  call check( upd_ok, 'every update ok, j-high' )
  call check( worst <= TOL_GCL, 'GCL residual, j-high' )
  write(*,'(A,ES10.2)') '   max relative GCL residual, j-high: ', worst
  pinned_ok = .true.
  do k = 0, NK ; do i = 0, NI
    if ( any( node(i,0,k)%c /= node0(i,0,k)%c ) ) pinned_ok = .false.
  enddo ; enddo
  call check( pinned_ok, 'pinned face j = 0 bitwise unchanged' )

  ! ---- 4. zero displacement: bitwise inert, exact zero swept volumes
  node0 = node
  do k = 0, NK ; do i = 0, NI
    law%disp(i,k)%c = 0.0_R8
  enddo ; enddo
  call Morph_Update( node, node_old, dim, t, 1.0_R8, law, geom, st )
  call check( st%ok(), 'update, zero displacement' )
  pinned_ok = .true.
  do k = 0, NK ; do j = 0, NJ ; do i = 0, NI
    if ( any( node(i,j,k)%c /= node0(i,j,k)%c ) ) pinned_ok = .false.
  enddo ; enddo ; enddo
  call check( pinned_ok .and. all( geom%dV_swept == 0.0_R8 ), &
              'zero displacement leaves the mesh and the swept volumes exactly unchanged' )

  ! ---- 5. a displacement equal to the line length: a status, not a stop
  do k = 0, NK ; do i = 0, NI
    law%disp(i,k)%c = [ 0.0_R8, -( node(i,NJ,k)%c(2) - node(i,0,k)%c(2) ), 0.0_R8 ]
  enddo ; enddo
  call Morph_Update( node, node_old, dim, t, 1.0_R8, law, geom, st )
  call check( st%code == MORPH_ERR_NEGVOL, 'disp = L_line returns MORPH_ERR_NEGVOL' )
  write(*,'(A,I0,A,A)') '   consumed line: status ', st%code, ', ', trim(st%message)
  deallocate( law%disp )

  write(*,*)
  if ( nfail == 0 ) then
    write(*,'(A)') 'ALL NODESHIFT CHECKS PASSED'
  else
    write(*,'(I0,A)') nfail, ' NODESHIFT CHECK(S) FAILED'
    error stop 1
  endif

contains

  subroutine check ( cond, what )
    logical,          intent(in) :: cond
    character(len=*), intent(in) :: what
    if ( cond ) then
      write(*,'(A,A)') '   PASS ', what
    else
      write(*,'(A,A)') '   FAIL ', what
      nfail = nfail + 1
    endif
  end subroutine check

  !> Graded in x (ratio 0.1 last/first, as the ATC meshes), uniform in y and z.
  subroutine build_mesh ( nd )
    type(morph_vec3_t), intent(inout) :: nd(0:,0:,0:)
    real(R8) :: x(0:NI), r, d1
    integer  :: ii, jj, kk
    r  = 0.1_R8 ** ( 1.0_R8 / ( NI - 1 ) )
    d1 = LX * ( 1.0_R8 - r ) / ( 1.0_R8 - r ** NI )
    x(0) = 0.0_R8
    do ii = 1, NI
      x(ii) = x(ii-1) + d1 * r ** ( ii - 1 )
    enddo
    x(NI) = LX
    do kk = 0, NK ; do jj = 0, NJ ; do ii = 0, NI
      nd(ii,jj,kk)%c = [ x(ii), LY * jj / NJ, LZ * kk / NK ]
    enddo ; enddo ; enddo
  end subroutine build_mesh

end program test_morph_nodeshift
