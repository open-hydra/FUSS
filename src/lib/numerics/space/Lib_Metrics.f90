module FUSS_Lib_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Base_Types_m
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Parameters_m

  implicit none
  real(R8), public :: delthe = 0d0   ! grid axisymmetric angle
  real(R8), public :: skew_factor_max = 1d0   ! largest time-step skew factor of the cells of the finest level (Skew_Factor)
  ! A degenerate or inverted cell found by Compute_Metric_Tensor, inside the OpenMP threads of Setup_Metrics: the run stops
  ! after the parallel region, on every rank (First_Bad_Metric, Stop_On_Bad_Metric)
  logical, public, protected :: bad_metric = .false.
  private :: First_Bad_Metric, Cross3

contains

  !> Mesh type, from the nodes of every block.
  !>  - ndir is the number of index directions with more than one cell: 3 if a block has more than one cell in k,
  !>    otherwise 2 if a block has more than one cell in j, otherwise 1. Until October 2026 only block 1 decided: a 3-D
  !>    mesh whose block 1 is one cell thick in k ran as 2-D, with block 1 decoupled across its k faces.
  !>  - delthe is the angle of the axisymmetric wedge of a 2-D mesh, 0 for a plane 2-D mesh and for 1-D and 3-D meshes
  !>    (FUSS has no rotational periodicity, type 200, so a 3-D mesh has no sector). A 2-D mesh is a wedge when, in every
  !>    block, the k = 1 node plane is the k = 0 plane turned about the x axis by the same angle; it is plane when no
  !>    block is so turned. Every node is checked (Node_Planes); blocks that disagree stop the run.
  !>  - The angle is the azimuth turned by the node (0,1) of the block, or by the node farthest from the axis when that
  !>    one lies on it. The angles of two blocks agree when they differ by at most ETA of the angle.
  !>  - The checks are relative to the move of the nodes per cell in k (K_Move_Per_Cell; the thickness of a 2-D block),
  !>    not to the size of the coordinates: the coordinates of a mesh are written with a finite number of digits, and
  !>    their round-off, relative to that move, is what a turn or a translation leaves over. A check relative to the
  !>    size of the coordinates, at 1e-6 of it, fails a wedge of 5 degrees written with 5 significant digits, which is
  !>    then taken for a plane mesh, and passes a plane slab far from the x axis, which is then taken for a wedge. A 2-D
  !>    block nearer a turn than a translation, or nearer a translation than a turn, but not to within ETA of its
  !>    thickness, stops the run: its nodes are written with too few digits to tell its type, or its k planes are neither
  !>    parallel nor turned about the x axis.
  subroutine Check_Mesh_Type( domain )
    use FUSS_Mod_MPI, only: mpi_stop_all
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    real(R8), parameter :: TOL = 1d-6, ETA = 1d-2, PI = 3.141592653589793d0
    integer  :: b, nb, nturn
    real(R8) :: sr, ang, ang_b, fit, fit_move
    logical  :: turned, moved, nearer_turn, nearer_move
    character(len=512) :: msg   ! room for the longest message and the block number

    nb = domain % nb

    ndir = 1
    do b = 1, nb
      if ( domain % Blk(b) % dim(2) > 1 ) ndir = max ( ndir, 2 )
      if ( domain % Blk(b) % dim(3) > 1 ) ndir = 3
    enddo

    delthe = 0d0
    if ( ndir /= 2 ) return

    ! Size of the distance from the x axis, to tell a probe node on the axis
    sr = 0d0
    do b = 1, nb
      associate ( node => domain % Blk(b) % node )
        sr = max ( sr, maxval ( sqrt ( node(:,:,:) % c(2)**2 + node(:,:,:) % c(3)**2 ) ) )
      end associate
    enddo

    ! Wedge if every block is turned by the same angle; plane if none is
    nturn = 0
    do b = 1, nb
      call Node_Planes ( domain % Blk(b), 1, sr, TOL, ETA, turned, ang_b, moved, fit, nearer_turn, fit_move, &
                         nearer_move )
      if ( nearer_turn .and. .not. turned ) then
        write(msg,'(A,I0,A,ES10.3,A,ES9.2,A)') '[ERROR] Mesh type: block ', b, ' is nearer a wedge of angle ', ang_b, &
          ' rad than a plane block, but its nodes match the turn only to ', fit, ' of its thickness (1e-2 needed):'// &
          ' write the node coordinates with more significant digits'
        call mpi_stop_all ( trim(msg) )
      end if
      if ( nearer_move .and. .not. moved ) then
        write(msg,'(A,I0,A,ES9.2,A)') '[ERROR] Mesh type: block ', b, ' is nearer a plane block than a wedge, but its'// &
          ' k = 1 node plane is the k = 0 plane translated only to ', fit_move, ' of its thickness (1e-2 needed):'// &
          ' its k planes must be parallel, or the coordinates written with more significant digits'
        call mpi_stop_all ( trim(msg) )
      end if
      if ( .not. turned ) cycle
      nturn = nturn + 1
      if ( nturn == 1 ) then
        ang = ang_b
      elseif ( abs ( ang_b - ang ) > ETA * abs ( ang ) ) then
        write(msg,'(A,I0,A,ES12.5,A,ES12.5,A)') '[ERROR] Mesh type: block ', b, ' is a wedge of angle ', ang_b, &
          ' rad, another block one of ', ang, ' rad'
        call mpi_stop_all ( trim(msg) )
      end if
    enddo
    if ( nturn == nb ) then
      delthe = ang
    elseif ( nturn > 0 ) then
      write(msg,'(A,I0,A,I0,A)') '[ERROR] Mesh type: ', nturn, ' of the ', nb, &
        ' blocks of the 2-D mesh are axisymmetric wedges (k = 1 plane turned about the x axis), the others are not'
      call mpi_stop_all ( trim(msg) )
    end if

  contains

    !> Compare the node planes k = 0 and k = ktop of a block. cell is the move of the nodes per cell in k measured along
    !> the k lines (K_Move_Per_Cell; in 2-D, ktop = 1, the thickness of the block). turned: every node of the plane ktop
    !> is the node of the plane 0 turned about the x axis by ang, with |ang| > tol, to within eta of cell, and closer
    !> than as a translation. moved: every node is the node of the plane 0 translated by the move of the node (0,0), to
    !> within eta of cell, and closer than as a turn, or with |ang| <= tol, which is no turn: the turn by a null angle is
    !> the translation by zero. Comparing the two keeps a wedge from passing for a translation and a slab far from the
    !> axis, which a small turn nearly fits, from passing for a wedge. A block whose probe node keeps its azimuth while
    !> its plane k = 1 is not the plane k = 0 translated to eta (planes that touch or cross, or a plane k = 1 scaled
    !> radially about the axis) stops here: without the null angle it would pass for plane. fit, fit_move: the misfits of
    !> the turn and of the translation over cell; nearer_turn: |ang| > tol and the turn fits better than the translation;
    !> nearer_move: the translation fits better than the turn, or |ang| <= tol.
    subroutine Node_Planes ( blk, ktop, sr, tol, eta, turned, ang, moved, fit, nearer_turn, fit_move, nearer_move )
      type(FUSS_block_type), intent(in) :: blk
      integer, intent(in)               :: ktop
      real(R8), intent(in)              :: sr, tol, eta
      logical, intent(out)              :: turned, moved, nearer_turn, nearer_move
      real(R8), intent(out)             :: ang, fit, fit_move
      ! Local
      integer  :: i, j, i0, j0
      real(R8) :: theta1, theta2, ca, sa, rmax, r, a(3), c(3), shift(3), turn_off, move_off, cell

      associate ( node => blk % node )
        ! Probe node: (0,1), or the node of the plane 0 farthest from the axis when (0,1) lies on it
        i0 = 0
        j0 = 1
        if ( sqrt ( node(0,1,0)%c(2)**2 + node(0,1,0)%c(3)**2 ) <= tol * sr ) then
          rmax = -1d0
          do j = 0, blk % dim(2)
            do i = 0, blk % dim(1)
              r = sqrt ( node(i,j,0)%c(2)**2 + node(i,j,0)%c(3)**2 )
              if ( r > rmax ) then
                rmax = r
                i0 = i
                j0 = j
              end if
            enddo
          enddo
        end if
        theta1 = atan2 ( node(i0,j0,0)%c(3), node(i0,j0,0)%c(2) )
        theta2 = atan2 ( node(i0,j0,ktop)%c(3), node(i0,j0,ktop)%c(2) )
        ang = theta2 - theta1
        if ( ang >  pi ) ang = ang - 2d0 * pi
        if ( ang < -pi ) ang = ang + 2d0 * pi

        ca = cos ( ang )
        sa = sin ( ang )
        shift = node(0,0,ktop)%c - node(0,0,0)%c
        turn_off = 0d0
        move_off = 0d0
        do j = 0, blk % dim(2)
          do i = 0, blk % dim(1)
            a = node(i,j,0)%c
            c = node(i,j,ktop)%c
            turn_off = max ( turn_off, abs ( c(1) - a(1) ), abs ( c(2) - ( a(2) * ca - a(3) * sa ) ), &
                                       abs ( c(3) - ( a(2) * sa + a(3) * ca ) ) )
            move_off = max ( move_off, maxval ( abs ( c - a - shift ) ) )
          enddo
        enddo
        cell = K_Move_Per_Cell ( blk, ktop )
        nearer_turn = abs ( ang ) > tol .and. turn_off < move_off
        nearer_move = move_off < turn_off .or. abs ( ang ) <= tol
        turned = nearer_turn .and. turn_off <= eta * cell
        moved  = nearer_move .and. move_off <= eta * cell
        fit      = huge ( 1d0 )
        fit_move = huge ( 1d0 )
        if ( cell > 0d0 ) fit = turn_off / cell
        if ( cell > 0d0 ) fit_move = move_off / cell
      end associate

    end subroutine Node_Planes

  end subroutine Check_Mesh_Type


  !> The move of the nodes of a block per cell in k, from the node plane 0 to the plane ktop, measured along the k
  !> lines: for every node of the plane 0 the sum, over the cells of its k line, of its move across each cell (the
  !> largest component), and the largest of these sums over ktop. With one cell in k (ktop = 1) it is the move from the
  !> plane 0 to the plane 1, bit for bit. The scale of the checks of Check_Mesh_Type.
  pure function K_Move_Per_Cell ( blk, ktop ) result ( cell )
    type(FUSS_block_type), intent(in) :: blk
    integer, intent(in)               :: ktop
    real(R8)                          :: cell
    ! Local
    integer  :: i, j, k
    real(R8) :: s

    cell = 0d0
    do j = 0, blk % dim(2)
      do i = 0, blk % dim(1)
        s = 0d0
        do k = 1, ktop
          s = s + maxval ( abs ( blk % node(i,j,k) % c - blk % node(i,j,k-1) % c ) )
        enddo
        cell = max ( cell, s )
      enddo
    enddo
    cell = cell / ktop

  end function K_Move_Per_Cell


  subroutine Compute_Norm_Area( b )
    implicit none
    class(FUSS_block_type) :: b  
    ! Local
    real(R8) :: Ai, snix, sniy, sniz, Aj, snjx, snjy, snjz, Ak, snkx, snky, snkz
    integer      :: i, j, k, im, jm, km
    real(R8) :: d1(3), d2(3), d3(3)
    real(R8) :: snixx, sniyy, snizz, snjxx, snjyy, snjzz, snkxx, snkyy, snkzz
    real(R8) :: scal, signi, signj, signk

    ! Peliminary operations
    im = b%dim(1)
    jm = b%dim(2)
    km = b%dim(3)
  
    ! compute sign of normal vectors to the intefaces
  
    ! i direction
    d1 = b%node(1,1,1)%c - b%node(1,0,0)%c
    d2 = b%node(1,0,1)%c - b%node(1,1,0)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snixx=d3(1)
    sniyy=d3(2)
    snizz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(1,0,1)%c + b%node(1,0,0)%c + b%node(1,1,0)%c )
    d2 = 0.25d0*( b%node(0,1,1)%c + b%node(0,0,1)%c + b%node(0,0,0)%c + b%node(0,1,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snixx+d3(2)*sniyy+d3(3)*snizz
  
    signi=sign(1.d0,scal)

    ! j direction
    d1 = b%node(0,1,1)%c - b%node(1,1,0)%c
    d2 = b%node(1,1,1)%c - b%node(0,1,0)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snjxx=d3(1)
    snjyy=d3(2)
    snjzz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(0,1,1)%c + b%node(0,1,0)%c + b%node(1,1,0)%c )
    d2 = 0.25d0*( b%node(1,0,1)%c + b%node(0,0,1)%c + b%node(0,0,0)%c + b%node(1,0,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snjxx+d3(2)*snjyy+d3(3)*snjzz
  
    signj=sign(1.d0,scal)

    ! k direction
    d1 = b%node(0,1,1)%c - b%node(1,0,1)%c
    d2 = b%node(0,0,1)%c - b%node(1,1,1)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snkxx=d3(1)
    snkyy=d3(2)
    snkzz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(0,1,1)%c + b%node(0,0,1)%c + b%node(1,0,1)%c )
    d2 = 0.25d0*( b%node(1,1,0)%c + b%node(0,1,0)%c + b%node(0,0,0)%c + b%node(1,0,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snkxx+d3(2)*snkyy+d3(3)*snkzz
  
    signk=sign(1.d0,scal)
    
    ! compute metrics: n, A
    
    !$omp parallel private (d1,d2,d3,i,j,k,snix,sniy,sniz,Ai,Aj,snjx,snjy,snjz,Ak,snkx,snky,snkz)
    
    ! i direction
    !$omp do collapse(3)
    do k = 1, km ; do j = 1, jm ; do i = 0, im
  
      d1 = b%node(i,j,k)%c - b%node(i,j-1,k-1)%c
      d2 = b%node(i,j-1,k)%c - b%node(i,j,k-1)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Ai = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      snix = d3(1)/Ai*signi
      sniy = d3(2)/Ai*signi
      sniz = d3(3)/Ai*signi
  
      if (Ai == 0d0) then
        snix = 0d0
        sniy = 0d0
        sniz = 0d0
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(1)%f(i,j,k)%a = Ai
      b%dir(1)%f(i,j,k)%n = [ snix, sniy, sniz ]
      
    enddo; enddo; enddo
  
    ! j direction
    !$omp do collapse(3)
    do k = 1, km ; do j = 0, jm ; do  i = 1, im
  
      d1 = b%node(i-1,j,k)%c - b%node(i,j,k-1)%c
      d2 = b%node(i,j,k)%c - b%node(i-1,j,k-1)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Aj = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      snjx = d3(1)/Aj*signj
      snjy = d3(2)/Aj*signj
      snjz = d3(3)/Aj*signj
  
      if (Aj == 0d0) then
        snjx = 0d0
        snjy = 0d0
        snjz = 0d0
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(2)%f(i,j,k)%A = Aj
      b%dir(2)%f(i,j,k)%n = [ snjx, snjy, snjz ]
  
    enddo; enddo; enddo
  
    ! k direction
    !$omp do collapse(3)
    do k = 0, km ; do j = 1, jm ; do i = 1, im
  
      d1 = b%node(i-1,j,k)%c - b%node(i,j-1,k)%c
      d2 = b%node(i-1,j-1,k)%c - b%node(i,j,k)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Ak = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      snkx = d3(1)/Ak*signk
      snky = d3(2)/Ak*signk
      snkz = d3(3)/Ak*signk
  
      if (Ak == 0d0) then
        snkx = 0d0
        snky = 0d0
        snkz = 0d0
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(3)%f(i,j,k)%A = Ak
      b%dir(3)%f(i,j,k)%n = [ snkx, snky, snkz ]
  
    enddo; enddo; enddo
    !$omp end parallel
  
  end subroutine Compute_Norm_Area


  !> True for the first degenerate or inverted cell of this rank, in any OpenMP thread: the flag is shared and set in a
  !> critical section, so that one thread writes the cell. A STOP there would end the process while the other threads go
  !> on (with two threads the exit code can be 0, the message of the cell lost or written twice).
  logical function First_Bad_Metric ( ) result ( first )
    !$omp critical (fuss_bad_metric)
    first = .not. bad_metric
    bad_metric = .true.
    !$omp end critical (fuss_bad_metric)
  end function First_Bad_Metric


  !> Every rank stops, exit code 1, when one of them found a degenerate or inverted cell (bad_metric). Each rank computes
  !> the cells of its own blocks, so the ranks take the decision together (the largest of a flag) before the stop.
  !> Setup_Metrics calls it on every rank, outside its parallel regions.
  subroutine Stop_On_Bad_Metric ( )
    use FUSS_Mod_MPI, only: mpi_stop_if_any
    implicit none
    call mpi_stop_if_any ( bad_metric, 'degenerate, inverted or folded cell (see above)', '' )
  end subroutine Stop_On_Bad_Metric


  subroutine Compute_Metric_Tensor( N1, N2, N3, N4, N5, N6, N7, N8, M, dl, vol, det1, loc, ratio )
    implicit none
    type(FUSS_vector_3D_type), intent(in)  :: N1, N2, N3, N4, N5, N6, N7, N8
    type(FUSS_tensor_3D_type), intent(out) :: M
    type(FUSS_vector_3D_type), intent(out) :: dl
    real(R8), intent(out)             :: vol
    real(R8), optional, intent(in)    :: det1     ! ghost cell: det of the first cell, whose sign and size it must keep
                                                  ! (0: no such check: the second layer, the first behind a face that
                                                  ! is not a chimera face, and behind a face of zero area)
    integer,  optional, intent(in)    :: loc(5)   ! block, cell i j k, face of a ghost cell (0 otherwise): stop message
    real(R8), optional, intent(in)    :: ratio    ! ghost cell: r = dl_2/dl_1 of the first two cells (stop message)
    ! Local
    integer  :: h
    real(R8) :: det, A(3,3), cofactor(3,3), scale
    logical  :: bad, first, check_vol
    real(R8) :: r(3,8), ea(3,3), v
    real(R8), volatile :: rv(3,8)   ! the nodes relative to N1, rounded once as written (the volume, below)

    ! Compute average cell dimension. A is M^-1, inverse of metric tensor.
    ! M1 = [ xcs, ycs, zcs ; xet, yet, zet ; xzi, yzi, zzi ].
    
    A(1,:) = N5 % c - N1 % c + &
             N6 % c - N2 % c + &
             N7 % c - N3 % c + &
             N8 % c - N4 % c
  
    A(2,:) = N3 % c - N1 % c + &
             N4 % c - N2 % c + &
             N7 % c - N5 % c + &
             N8 % c - N6 % c
             
    A(3,:) = N2 % c - N1 % c + &
             N4 % c - N3 % c + &
             N6 % c - N5 % c + &
             N8 % c - N7 % c

    A = 0.25d0 * A

    ! Determinant of M^-1.
    det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

    ! Volume of the cell: the integral of det J of the trilinear map of its nodes, whose faces are the bilinear faces of
    ! Compute_Norm_Area, shared with the cells next door (docs/theory/numerics.md, Grid Metrics):
    !   V = 1/6 sum_d e_d . (S_d- + S_d+),
    ! e_d the mean edge of direction d and S_d-, S_d+ the area vectors, (diagonal x diagonal)/2, of the two faces of
    ! direction d, both along +d. It is exact on every cell, twisted faces included, and the volumes of the cells sum to the
    ! volume enclosed by the boundary faces (the five tetrahedra it replaces were not, and were wrong on the cells with
    ! twisted faces). Signed: negative on a left-handed cell, stored in absolute value below.
    ! Only differences of nodes enter, so the rounding is that of the size of the cell, not of its coordinates, as long as
    ! the compiler rounds the differences as written. With -ffast-math (gfortran -fassociative-math, ifx -fp-model fast) it
    ! may rewrite a difference of differences across statements, (q - c1) - (c2 - c1) as q - c2 (test/unit/
    ! volatile_order_check.f90). A VOLATILE object is defined and referenced each time as the program says (Fortran 2018,
    ! 8.5.19), so the nodes relative to N1 are stored in rv and read back into r, and every later expression has terms of
    ! the size of the cell only. A pure procedure cannot hold VOLATILE locals (C1589).
    rv(:,1) = 0d0
    rv(:,2) = N2 % c - N1 % c
    rv(:,3) = N3 % c - N1 % c
    rv(:,4) = N4 % c - N1 % c
    rv(:,5) = N5 % c - N1 % c
    rv(:,6) = N6 % c - N1 % c
    rv(:,7) = N7 % c - N1 % c
    rv(:,8) = N8 % c - N1 % c
    r = rv
    ea(:,1) = 0.25d0 * ( r(:,5) + ( r(:,6) - r(:,2) ) + ( r(:,7) - r(:,3) ) + ( r(:,8) - r(:,4) ) )
    ea(:,2) = 0.25d0 * ( r(:,3) + ( r(:,4) - r(:,2) ) + ( r(:,7) - r(:,5) ) + ( r(:,8) - r(:,6) ) )
    ea(:,3) = 0.25d0 * ( r(:,2) + ( r(:,4) - r(:,3) ) + ( r(:,6) - r(:,5) ) + ( r(:,8) - r(:,7) ) )
    ! The faces, corners in the cyclic order that turns about +d: i N1 N3 N4 N2 and N5 N7 N8 N6, j N1 N2 N6 N5 and
    ! N3 N4 N8 N7, k N1 N5 N7 N3 and N2 N6 N8 N4; twice their area vectors are (3rd - 1st) x (4th - 2nd)
    v =     dot_product ( ea(:,1), Cross3 ( r(:,4), r(:,2) - r(:,3) ) + Cross3 ( r(:,8) - r(:,5), r(:,6) - r(:,7) ) )
    v = v + dot_product ( ea(:,2), Cross3 ( r(:,6), r(:,5) - r(:,2) ) + Cross3 ( r(:,8) - r(:,3), r(:,7) - r(:,4) ) )
    v = v + dot_product ( ea(:,3), Cross3 ( r(:,7), r(:,3) - r(:,5) ) + Cross3 ( r(:,8) - r(:,2), r(:,4) - r(:,6) ) )
    v = v / 12d0

    ! Degenerate cell: |det| against the product of the lengths of the rows of A, its bound (Hadamard's inequality:
    ! equal for an orthogonal cell, whatever its stretching). A ghost cell of the first layer behind a chimera face must
    ! also keep the sign and the size of det of the first cell, to 1e-12 (a block may be left-handed: Compute_Norm_Area
    ! takes the orientation of the normals from its first cell): a row of A that rounding leaves tiny but still in its
    ! direction passes the first test, not this one. A mesh error, or ghost cells extrapolated from a wall-thin first cell
    ! whose coordinates were read back in single precision.
    ! This runs inside the OpenMP threads of Setup_Metrics, where a STOP ends the process while the other threads go on:
    ! a bad cell raises bad_metric, the first thread of the rank that finds one writes it, M is the identity (finite for
    ! what the region still computes), and Setup_Metrics stops every rank after the region (Stop_On_Bad_Metric).
    scale = norm2 ( A(1,:) ) * norm2 ( A(2,:) ) * norm2 ( A(3,:) )
    bad = .not. ( abs(det) > 1d-12 * scale )
    if ( present(det1) ) then
      if ( det1 /= 0d0 ) bad = bad .or. .not. ( det * det1 > 1d-12 * det1**2 )
    end if
    ! The volume has the sign of det: det J has the sign of the block in every point of a cell that is not turned over,
    ! and so has its integral. A cell turned over away from its centre can keep det at the centre and give a volume of
    ! the other sign or zero. Not checked where det1 = 0 (the second ghost layer, the first behind a face that is not a
    ! chimera face, and behind a face of zero area), whose size alone is checked.
    check_vol = .true.
    if ( present(det1) ) check_vol = det1 /= 0d0
    if ( check_vol ) bad = bad .or. .not. ( v * det > 0d0 )
    first = .false.
    if ( bad ) then
      first = First_Bad_Metric ( )
      M % c = 0d0
      do h = 1, 3
        M % c(h,h) = 1d0
      enddo
    end if
    if ( first ) then
      if ( present(loc) ) then
        if ( loc(5) > 0 ) then
          write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] Metric tensor: block ', loc(1), ', cell', loc(2:4), &
                                               ', ghost cell behind face ', loc(5), ': degenerate or inverted'
          if ( present(ratio) ) write(*,'(A,F0.4)') ' [ERROR]   r = dl_2/dl_1 of the first two cells across the face: ', ratio
        else
          write(*,'(A,I0,A,3(1X,I0),A)') ' [ERROR] Metric tensor: block ', loc(1), ', cell', loc(2:4), &
                                          ': degenerate or inverted'
        end if
      else
        write(*,'(A)') ' [ERROR] Metric tensor: degenerate or inverted cell'
      end if
      write(*,'(A,ES12.4,A,ES12.4,A)') ' [ERROR]   det = ', det, ', |A1| |A2| |A3| = ', scale, ' (equal for an orthogonal cell)'
      write(*,'(A,ES12.4,A)') ' [ERROR]   volume = ', v, ' (the integral of det J, of the sign of det on a cell not turned over)'
      if ( present(det1) ) then
        if ( det1 /= 0d0 ) write(*,'(A,ES12.4)') ' [ERROR]   det of the first cell next to the face: ', det1
        write(*,'(A)') ' [ERROR]   The ghost cells take the edge across the face of the first cell, and their edges along the'
        write(*,'(A)') ' [ERROR]   face are extrapolated: behind a chimera face (102), whose flux takes its metric, the first'
        write(*,'(A)') ' [ERROR]   layer turns over where those edges change sharply across the first cell (on a plane 2-D'
        write(*,'(A)') ' [ERROR]   mesh, growing 3 times or more); behind any face a layer of zero size stops too.'
        write(*,'(A)') ' [ERROR]   Remedy: index lines that change smoothly across the first two node layers at that face.'
      else
        write(*,'(A)') ' [ERROR]   Check the mesh, or the precision of the restart file.'
      end if
      write(*,'(A,3(1X,ES14.6))') ' [ERROR]   corner node N1:', N1 % c
      write(*,'(A,3(1X,ES14.6))') ' [ERROR]   corner node N8:', N8 % c

    else if ( .not. bad ) then
      
      cofactor(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))
      cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      cofactor(1,3) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))
      cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      cofactor(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))
      cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      cofactor(3,1) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))
      cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      cofactor(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))

      M % c = cofactor / det

    endif

    ! Cell length in i,j,k
    do h = 1, 3
      dl % c(h) = sqrt( A(h,1)**2 + A(h,2)**2 + A(h,3)**2 )
    enddo

    vol = abs ( v )

  end subroutine Compute_Metric_Tensor


  pure function Cross3 ( u, v ) result ( w )
    real(R8), intent(in) :: u(3), v(3)
    real(R8)             :: w(3)
    w = [ u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1) ]
  end function Cross3


  !> Metrics of the gc ghost cells behind face Fm of cell (Im,Jm,Km) of block b, from nodes extrapolated linearly along
  !> the index lines, N(-g) = 2 N(-g+1) - N(-g+2) (MOSE cc4caf1; the cubic 3N - 3N + N of 237c9c2 is gone): every ghost
  !> layer has the edge across the face of the first cell, so its size across the face never vanishes nor turns over;
  !> its edges along the face are extrapolated too, and turn a layer over where they change sharply across the first
  !> cell (on a plane 2-D mesh, growing 3 times or more turns the first layer, 5/3 times the second). Behind a chimera
  !> face (chimera true) the first layer must keep the sign and the size of the first cell (Compute_Metric_Tensor
  !> stops otherwise): its metric is the ghost metric of the flux of chimera faces (BC_Connection takes bc % Mg(1)), the
  !> one flux that reads the metric of an extrapolated ghost cell. Behind the other faces no flux reads it (the walls take
  !> the metric of the cell, symmetry copies it, a connection takes the cell across): only a zero size stops there. The
  !> first layer is centred at 2 x_f - x_1 (x_f the mean of the nodes of the face, x_1 that of the cell): the position of
  !> the ghost value of a wall, 2 Tw - T. The 2-D axisymmetric wedge rotates the cell (exact).
  subroutine BC_Extrapolate_Metrics( Im, Jm, Km, Fm, blk, Mg, dlg, volg, b, chimera )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout) :: blk
    type(FUSS_tensor_3D_type), intent(out)    :: Mg(2)
    type(FUSS_vector_3D_type), intent(out)    :: dlg(2)
    real(R8), intent(out)                :: volg(2)
    integer, intent(in)                  :: b     ! the block, for the stop message of Compute_Metric_Tensor
    logical, intent(in)                  :: chimera   ! a chimera face (102): the sign of the first layer is checked
    ! Local
    integer :: g, loc(5), d, x(3)
    real(R8) :: d1, rr
    type(FUSS_vector_3D_type), dimension(-gc:gc):: N1, N2, N3, N4, N5, N6, N7, N8

    ! det of the first cell, 1/det(M), whose sign and size its ghost cells must keep (Compute_Metric_Tensor)
    associate ( Mc => blk % M(Im,Jm,Km) % c )
      d1 = 1d0 / ( Mc(1,1) * ( Mc(2,2)*Mc(3,3) - Mc(2,3)*Mc(3,2) ) - Mc(1,2) * ( Mc(2,1)*Mc(3,3) - Mc(2,3)*Mc(3,1) ) &
                 + Mc(1,3) * ( Mc(2,1)*Mc(3,2) - Mc(2,2)*Mc(3,1) ) )
    end associate
    loc = [b, Im, Jm, Km, Fm]
    ! r = dl_2/dl_1 of the first two cells across the face, for the stop message (0 with one cell)
    d  = (Fm + 1) / 2
    rr = 0d0
    if ( blk % dim(d) >= 2 ) rr = blk % dl(Im+guide(Fm,1),Jm+guide(Fm,2),Km+guide(Fm,3)) % c(d) / blk % dl(Im,Jm,Km) % c(d)
    ! Behind a face of zero area (an axis, as Compute_Norm_Area finds it) the ghost cells are the reflection of the
    ! first cell through the axis: det = -d1 by construction. No check of their sign there (d1 = 0)
    x = [Im, Jm, Km]
    if ( mod(Fm,2) == 1 ) x(d) = x(d) - 1
    if ( blk % dir(d) % f(x(1),x(2),x(3)) % a == 0d0 ) d1 = 0d0
    ! Behind a face that is not a chimera face no flux reads the metric of the first layer: no check of its sign
    ! either, as for the second layer (docs/theory/numerics.md, Ghost Cells)
    if ( .not. chimera ) d1 = 0d0

    select case (Fm)
      case(1)
        do g = 0, 1
          N1(g)%c = blk % node (Im+g-1,Jm-1,Km-1) % c
          N2(g)%c = blk % node (Im+g-1,Jm-1,Km  ) % c
          N3(g)%c = blk % node (Im+g-1,Jm  ,Km-1) % c
          N4(g)%c = blk % node (Im+g-1,Jm  ,Km  ) % c
        enddo
        do g = -1, -gc, -1
          N1(g)%c = 2d0*N1(g+1)%c - N1(g+2)%c
          N2(g)%c = 2d0*N2(g+1)%c - N2(g+2)%c
          N3(g)%c = 2d0*N3(g+1)%c - N3(g+2)%c
          N4(g)%c = 2d0*N4(g+1)%c - N4(g+2)%c
          call Compute_Metric_Tensor ( N1(g), N2(g), N3(g), N4(g), N1(g+1), N2(g+1), N3(g+1), N4(g+1), &
                                       blk % M(Im+g, Jm, Km), blk % dl(Im+g, Jm, Km), blk % vol(Im+g, Jm, Km), &
                                       det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
          Mg(-g)   = blk % M(Im+g, Jm, Km)
          dlg(-g)  = blk % dl(Im+g, Jm, Km)
          volg(-g) = blk % vol(Im+g, Jm, Km)
        enddo
        
      case(2)
        do g = 0, -1, -1
          N5(g)%c = blk % node (Im+g,Jm-1,Km-1) % c
          N6(g)%c = blk % node (Im+g,Jm-1,Km  ) % c
          N7(g)%c = blk % node (Im+g,Jm  ,Km-1) % c
          N8(g)%c = blk % node (Im+g,Jm  ,Km  ) % c
        enddo
        do g = 1, gc
          N5(g)%c = 2d0*N5(g-1)%c - N5(g-2)%c
          N6(g)%c = 2d0*N6(g-1)%c - N6(g-2)%c
          N7(g)%c = 2d0*N7(g-1)%c - N7(g-2)%c
          N8(g)%c = 2d0*N8(g-1)%c - N8(g-2)%c
          call Compute_Metric_Tensor ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                       blk % M(Im+g, Jm, Km), blk % dl(Im+g, Jm, Km), blk % vol(Im+g, Jm, Km), &
                                       det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
          Mg(g)   = blk % M(Im+g, Jm, Km)
          dlg(g)  = blk % dl(Im+g, Jm, Km)
          volg(g) = blk % vol(Im+g, Jm, Km)
        enddo

      case(3)
        do g = 0, 1
          N1(g)%c = blk % node (Im-1,Jm+g-1,Km-1) % c
          N2(g)%c = blk % node (Im-1,Jm+g-1,Km  ) % c
          N5(g)%c = blk % node (Im  ,Jm+g-1,Km-1) % c
          N6(g)%c = blk % node (Im  ,Jm+g-1,Km  ) % c
        enddo
        do g = -1, -gc, -1
          N1(g)%c = 2d0*N1(g+1)%c - N1(g+2)%c
          N2(g)%c = 2d0*N2(g+1)%c - N2(g+2)%c
          N5(g)%c = 2d0*N5(g+1)%c - N5(g+2)%c
          N6(g)%c = 2d0*N6(g+1)%c - N6(g+2)%c
          call Compute_Metric_Tensor ( N1(g), N2(g), N1(g+1), N2(g+1), N5(g), N6(g), N5(g+1), N6(g+1), &
                                       blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km), &
                                       det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
          Mg(-g)   = blk % M(Im, Jm+g, Km)
          dlg(-g)  = blk % dl(Im, Jm+g, Km)
          volg(-g) = blk % vol(Im, Jm+g, Km)
        enddo

      case(4)
        do g = 0, -1, -1
          N3(g)%c = blk % node (Im-1,Jm+g,Km-1) % c
          N4(g)%c = blk % node (Im-1,Jm+g,Km  ) % c
          N7(g)%c = blk % node (Im  ,Jm+g,Km-1) % c
          N8(g)%c = blk % node (Im  ,Jm+g,Km  ) % c
        enddo
        do g = 1, gc
          N3(g)%c = 2d0*N3(g-1)%c - N3(g-2)%c
          N4(g)%c = 2d0*N4(g-1)%c - N4(g-2)%c
          N7(g)%c = 2d0*N7(g-1)%c - N7(g-2)%c
          N8(g)%c = 2d0*N8(g-1)%c - N8(g-2)%c
          call Compute_Metric_Tensor ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                       blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km), &
                                       det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
          Mg(g)   = blk % M(Im, Jm+g, Km)
          dlg(g)  = blk % dl(Im, Jm+g, Km)
          volg(g) = blk % vol(Im, Jm+g, Km)
        enddo

      case(5)
        if ( ndir==2 .and. delthe /= 0d0 ) then
          ! 2Dax: extrapolation with a rotation angle delthe (exact).
          do g = 0, 1
            N1(g)%c = blk % node (Im-1,Jm-1,Km+g-1) % c
            N3(g)%c = blk % node (Im-1,Jm  ,Km+g-1) % c
            N5(g)%c = blk % node (Im  ,Jm-1,Km+g-1) % c
            N7(g)%c = blk % node (Im  ,Jm  ,Km+g-1) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c(1) = N1(0) % c(1) 
            N1(g)%c(2) = N1(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N1(g)%c(3) = N1(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N3(g)%c(1) = N3(0) % c(1) 
            N3(g)%c(2) = N3(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N3(g)%c(3) = N3(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N5(g)%c(1) = N5(0) % c(1) 
            N5(g)%c(2) = N5(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N5(g)%c(3) = N5(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N7(g)%c(1) = N7(0) % c(1) 
            N7(g)%c(2) = N7(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N7(g)%c(3) = N7(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            call Compute_Metric_Tensor ( N1(g), N1(g+1), N3(g), N3(g+1), N5(g), N5(g+1), N7(g), N7(g+1), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g), &
                                         det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo         
        else
          do g = 0, 1
            N1(g)%c = blk % node (Im-1,Jm-1,Km+g-1) % c
            N3(g)%c = blk % node (Im-1,Jm  ,Km+g-1) % c
            N5(g)%c = blk % node (Im  ,Jm-1,Km+g-1) % c
            N7(g)%c = blk % node (Im  ,Jm  ,Km+g-1) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 2d0*N1(g+1)%c - N1(g+2)%c
            N3(g)%c = 2d0*N3(g+1)%c - N3(g+2)%c
            N5(g)%c = 2d0*N5(g+1)%c - N5(g+2)%c
            N7(g)%c = 2d0*N7(g+1)%c - N7(g+2)%c
            call Compute_Metric_Tensor ( N1(g), N1(g+1), N3(g), N3(g+1), N5(g), N5(g+1), N7(g), N7(g+1), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g), &
                                         det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo
        endif

      case(6)
        if ( ndir==2 .and. delthe /= 0d0 ) then
          ! 2Dax: extrapolation with a rotation angle delthe (exact).
          do g = 0, -1, -1
            N2(g)%c = blk % node (Im-1,Jm-1,Km+g) % c
            N4(g)%c = blk % node (Im-1,Jm  ,Km+g) % c
            N6(g)%c = blk % node (Im  ,Jm-1,Km+g) % c
            N8(g)%c = blk % node (Im  ,Jm  ,Km+g) % c
          enddo
          do g = 1, gc
            N2(g)%c(1) = N2(0) % c(1) 
            N2(g)%c(2) = N2(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N2(g)%c(3) = N2(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N4(g)%c(1) = N4(0) % c(1) 
            N4(g)%c(2) = N4(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N4(g)%c(3) = N4(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N6(g)%c(1) = N6(0) % c(1) 
            N6(g)%c(2) = N6(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N6(g)%c(3) = N6(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N8(g)%c(1) = N8(0) % c(1) 
            N8(g)%c(2) = N8(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N8(g)%c(3) = N8(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g), &
                                         det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo         
        else
          do g = 0, -1, -1
            N2(g)%c = blk % node (Im-1,Jm-1,Km+g) % c
            N4(g)%c = blk % node (Im-1,Jm  ,Km+g) % c
            N6(g)%c = blk % node (Im  ,Jm-1,Km+g) % c
            N8(g)%c = blk % node (Im  ,Jm  ,Km+g) % c
          enddo
          do g = 1, gc
            N2(g)%c = 2d0*N2(g-1)%c - N2(g-2)%c
            N4(g)%c = 2d0*N4(g-1)%c - N4(g-2)%c
            N6(g)%c = 2d0*N6(g-1)%c - N6(g-2)%c
            N8(g)%c = 2d0*N8(g-1)%c - N8(g-2)%c
            call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g), &
                                         det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo
        endif

      end select

  end subroutine BC_Extrapolate_Metrics


  !> Layer g of the ghost cells of a symmetry face (300): a copy of the cell g - 1 behind the face, in M, dl and vol of
  !> the block and in the copies of the record, Mg(g), dlg(g), volg(g). In a direction with fewer cells than g that cell
  !> is a ghost cell of the opposite face: Setup_Metrics fills the layers one at a time, every face of a layer before
  !> the next (MOSE f63f9af; layers 1 and 2 of faces 5 and 6 of a one-cell-thick mesh raced in one parallel loop).
  subroutine BC_Symmetry_Metrics( Im, Jm, Km, Fm, blk, Mg, dlg, volg, g )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm, g
    type(FUSS_block_type), intent(inout)     :: blk
    type(FUSS_tensor_3D_type), intent(inout) :: Mg(2)
    type(FUSS_vector_3D_type), intent(inout) :: dlg(2)
    real(R8), intent(inout)                  :: volg(2)

    select case (Fm)
      case(1)
        blk % M(Im-g,Jm,Km)   = blk % M(Im+g-1,Jm,Km)
        blk % dl(Im-g,Jm,Km)  = blk % dl(Im+g-1,Jm,Km)
        blk % vol(Im-g,Jm,Km) = blk % vol(Im+g-1,Jm,Km)
        Mg(g)   = blk % M(Im-g,Jm,Km)
        dlg(g)  = blk % dl(Im-g,Jm,Km)
        volg(g) = blk % vol(Im-g,Jm,Km)
      case(2)
        blk % M(Im+g,Jm,Km)   = blk % M(Im-g+1,Jm,Km)
        blk % dl(Im+g,Jm,Km)  = blk % dl(Im-g+1,Jm,Km)
        blk % vol(Im+g,Jm,Km) = blk % vol(Im-g+1,Jm,Km)
        Mg(g)   = blk % M(Im+g,Jm,Km)
        dlg(g)  = blk % dl(Im+g,Jm,Km)
        volg(g) = blk % vol(Im+g,Jm,Km)
      case(3)
        blk % M(Im,Jm-g,Km)   = blk % M(Im,Jm+g-1,Km)
        blk % dl(Im,Jm-g,Km)  = blk % dl(Im,Jm+g-1,Km)
        blk % vol(Im,Jm-g,Km) = blk % vol(Im,Jm+g-1,Km)
        Mg(g)   = blk % M(Im,Jm-g,Km)
        dlg(g)  = blk % dl(Im,Jm-g,Km)
        volg(g) = blk % vol(Im,Jm-g,Km)
      case(4)
        blk % M(Im,Jm+g,Km)   = blk % M(Im,Jm-g+1,Km)
        blk % dl(Im,Jm+g,Km)  = blk % dl(Im,Jm-g+1,Km)
        blk % vol(Im,Jm+g,Km) = blk % vol(Im,Jm-g+1,Km)
        Mg(g)   = blk % M(Im,Jm+g,Km)
        dlg(g)  = blk % dl(Im,Jm+g,Km)
        volg(g) = blk % vol(Im,Jm+g,Km)
      case(5)
        blk % M(Im,Jm,Km-g)   = blk % M(Im,Jm,Km+g-1)
        blk % dl(Im,Jm,Km-g)  = blk % dl(Im,Jm,Km+g-1)
        blk % vol(Im,Jm,Km-g) = blk % vol(Im,Jm,Km+g-1)
        Mg(g)   = blk % M(Im,Jm,Km-g)
        dlg(g)  = blk % dl(Im,Jm,Km-g)
        volg(g) = blk % vol(Im,Jm,Km-g)
      case(6)
        blk % M(Im,Jm,Km+g)   = blk % M(Im,Jm,Km-g+1)
        blk % dl(Im,Jm,Km+g)  = blk % dl(Im,Jm,Km-g+1)
        blk % vol(Im,Jm,Km+g) = blk % vol(Im,Jm,Km-g+1)
        Mg(g)   = blk % M(Im,Jm,Km+g)
        dlg(g)  = blk % dl(Im,Jm,Km+g)
        volg(g) = blk % vol(Im,Jm,Km+g)
    end select

  end subroutine BC_Symmetry_Metrics


  subroutine BC_Connect_Metrics ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, d11s, d12s, d21s, d22s, Mg, dlg, volg )
    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(FUSS_block_type), intent(in)    :: blks
    type(FUSS_block_type), intent(inout) :: blkm
    type(FUSS_tensor_3D_type), intent(out) :: Mg(2)
    type(FUSS_vector_3D_type), intent(out) :: dlg(2)
    real(R8), intent(out)                  :: volg(2)
    ! Local
    integer      :: g, Is1, Js1, Ks1
    integer      :: II, JJ, KK, III, JJJ, KKK      
    integer      :: guidem(6,3), guides(6,3), guidem2(6,3), guides2(6,3)
    type(FUSS_vector_3D_type), dimension(0:gc) :: N1, N2, N3, N4, N5, N6, N7, N8
    
    ! ---------------------------------------------------------------------------------------------
    ! Preliminary definitions
    guidem(1,1)=-2 ; guidem(1,2)= 0 ; guidem(1,3)= 0 ; guidem(2,1)=1 ; guidem(2,2)=0 ; guidem(2,3)=0
    guidem(3,1)= 0 ; guidem(3,2)=-2 ; guidem(3,3)= 0 ; guidem(4,1)=0 ; guidem(4,2)=1 ; guidem(4,3)=0
    guidem(5,1)= 0 ; guidem(5,2)=0  ; guidem(5,3)=-2 ; guidem(6,1)=0 ; guidem(6,2)=0 ; guidem(6,3)=1

    guidem2(1,1)=-3 ; guidem2(1,2)= 0 ; guidem2(1,3)= 0 ; guidem2(2,1)=2 ; guidem2(2,2)=0 ; guidem2(2,3)=0
    guidem2(3,1)= 0 ; guidem2(3,2)=-3 ; guidem2(3,3)= 0 ; guidem2(4,1)=0 ; guidem2(4,2)=2 ; guidem2(4,3)=0
    guidem2(5,1)= 0 ; guidem2(5,2)= 0 ; guidem2(5,3)=-3 ; guidem2(6,1)=0 ; guidem2(6,2)=0 ; guidem2(6,3)=2

    guides(1,1)= 0 ; guides(1,2)= 0 ; guides(1,3)= 0 ; guides(2,1)=-1 ; guides(2,2)= 0 ; guides(2,3)= 0
    guides(3,1)= 0 ; guides(3,2)= 0 ; guides(3,3)= 0 ; guides(4,1)= 0 ; guides(4,2)=-1 ; guides(4,3)= 0
    guides(5,1)= 0 ; guides(5,2)= 0 ; guides(5,3)= 0 ; guides(6,1)= 0 ; guides(6,2)= 0 ; guides(6,3)=-1

    guides2(1,1)= 1 ; guides2(1,2)= 0 ; guides2(1,3)= 0 ; guides2(2,1)=-2 ; guides2(2,2)= 0 ; guides2(2,3)= 0
    guides2(3,1)= 0 ; guides2(3,2)= 1 ; guides2(3,3)= 0 ; guides2(4,1)= 0 ; guides2(4,2)=-2 ; guides2(4,3)= 0
    guides2(5,1)= 0 ; guides2(5,2)= 0 ; guides2(5,3)= 1 ; guides2(6,1)= 0 ; guides2(6,2)= 0 ; guides2(6,3)=-2
    ! ---------------------------------------------------------------------------------------------

    select case (Fm)
      case(1)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
      case(2)
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(3)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
      case(4)
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(5)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
      case(6)
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
    end select

    do g = 1, gc

      ! connected cell node index
      Is1 = Is + (gc-g)*guides(Fs,1)+(g-1)*guides2(Fs,1)
      Js1 = Js + (gc-g)*guides(Fs,2)+(g-1)*guides2(Fs,2)
      Ks1 = Ks + (gc-g)*guides(Fs,3)+(g-1)*guides2(Fs,3)

      ! compute the indexes in the connected block corresponding to the ghost nodes to update in the local block 
      select case(Fs)
        case(1:2) ! i-faces
          II=Is1
          JJ=(2*Js1+d11s+d21s)/2
          KK=(2*Ks1+d12s+d22s)/2
        case(3:4) ! j-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=Js1
          KK=(2*Ks1+d12s+d22s)/2
        case(5:6) ! k-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=(2*Js1+d12s+d22s)/2
          KK=Ks1
      end select
      if (Fm == 1) then
          N4(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 2) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 3) then
          N6(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 4) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 5) then
          N7(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 6) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      endif

      ! update the remaining 3 indexes (only necessary in the corners, but for good measure)
      ! node i1-1,i2 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s
          KKK=KK-d12s
        case(3:4)
          III=II-d11s
          JJJ=JJ
          KKK=KK-d12s
        case(5:6)
          III=II-d11s
          JJJ=JJ-d12s
          KKK=KK
      end select
      if (Fm == 1) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! node i1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d21s
          KKK=KK-d22s
        case(3:4)
          III=II-d21s
          JJJ=JJ
          KKK=KK-d22s
        case(5:6)
          III=II-d21s
          JJJ=JJ-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! node i1-1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s-d21s
          KKK=KK-d12s-d22s
        case(3:4)
          III=II-d11s-d21s
          JJJ=JJ
          KKK=KK-d12s-d22s
        case(5:6)
          III=II-d11s-d21s
          JJJ=JJ-d12s-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! compute metric variables
      select case(Fm)
        case(1)
          call Compute_Metric_Tensor ( N1(g), N2(g), N3(g), N4(g), N1(g-1), N2(g-1), N3(g-1), N4(g-1), &
                                       blkm%M(Im-g, Jm, Km), blkm%dl(Im-g, Jm, Km), blkm%vol(Im-g, Jm, Km))
          Mg(g)   = blkm % M(Im-g, Jm, Km)
          dlg(g)  = blkm % dl(Im-g, Jm, Km)
          volg(g) = blkm % vol(Im-g, Jm, Km)
        case(2)
          call Compute_Metric_Tensor ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                       blkm%M(Im+g, Jm, Km), blkm%dl(Im+g, Jm, Km), blkm%vol(Im+g, Jm, Km))
          Mg(g)   = blkm % M(Im+g, Jm, Km)
          dlg(g)  = blkm % dl(Im+g, Jm, Km)
          volg(g) = blkm % vol(Im+g, Jm, Km)
        case(3)
          call Compute_Metric_Tensor ( N1(g), N2(g), N1(g-1), N2(g-1), N5(g), N6(g), N5(g-1), N6(g-1), &
                                       blkm%M(Im, Jm-g, Km), blkm%dl(Im, Jm-g, Km), blkm%vol(Im, Jm-g, Km))
          Mg(g)   = blkm % M(Im, Jm-g, Km)
          dlg(g)  = blkm % dl(Im, Jm-g, Km)
          volg(g) = blkm % vol(Im, Jm-g, Km)
        case(4)
          call Compute_Metric_Tensor ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                       blkm%M(Im, Jm+g, Km), blkm%dl(Im, Jm+g, Km), blkm%vol(Im, Jm+g, Km))
          Mg(g)   = blkm % M(Im, Jm+g, Km)
          dlg(g)  = blkm % dl(Im, Jm+g, Km)
          volg(g) = blkm % vol(Im, Jm+g, Km)
        case(5)
          call Compute_Metric_Tensor ( N1(g), N1(g-1), N3(g), N3(g-1), N5(g), N5(g-1), N7(g), N7(g-1), &
                                       blkm%M(Im, Jm, Km-g), blkm%dl(Im, Jm, Km-g), blkm%vol(Im, Jm, Km-g))
          Mg(g)   = blkm % M(Im, Jm, Km-g)
          dlg(g)  = blkm % dl(Im, Jm, Km-g)
          volg(g) = blkm % vol(Im, Jm, Km-g)
        case(6)
          call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                       blkm%M(Im, Jm, Km+g), blkm%dl(Im, Jm, Km+g), blkm%vol(Im, Jm, Km+g))
          Mg(g)   = blkm % M(Im, Jm, Km+g)
          dlg(g)  = blkm % dl(Im, Jm, Km+g)
          volg(g) = blkm % vol(Im, Jm, Km+g)
        end select

    enddo ! End loop over ghost cell           

  end subroutine BC_Connect_Metrics

end module FUSS_Lib_Metrics