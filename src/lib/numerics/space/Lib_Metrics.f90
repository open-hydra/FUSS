module FUSS_Lib_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Base_Types_m
  use FUSS_Advanced_Types_m
  use FUSS_Global_m
  use FUSS_Parameters_m

  implicit none
  real(R8), public :: delthe = 0d0   ! grid axisymmetric angle

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
    use FUSS_Mod_MPI, only: mpi_abort_all
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
        call mpi_abort_all ( trim(msg) )
      end if
      if ( nearer_move .and. .not. moved ) then
        write(msg,'(A,I0,A,ES9.2,A)') '[ERROR] Mesh type: block ', b, ' is nearer a plane block than a wedge, but its'// &
          ' k = 1 node plane is the k = 0 plane translated only to ', fit_move, ' of its thickness (1e-2 needed):'// &
          ' its k planes must be parallel, or the coordinates written with more significant digits'
        call mpi_abort_all ( trim(msg) )
      end if
      if ( .not. turned ) cycle
      nturn = nturn + 1
      if ( nturn == 1 ) then
        ang = ang_b
      elseif ( abs ( ang_b - ang ) > ETA * abs ( ang ) ) then
        write(msg,'(A,I0,A,ES12.5,A,ES12.5,A)') '[ERROR] Mesh type: block ', b, ' is a wedge of angle ', ang_b, &
          ' rad, another block one of ', ang, ' rad'
        call mpi_abort_all ( trim(msg) )
      end if
    enddo
    if ( nturn == nb ) then
      delthe = ang
    elseif ( nturn > 0 ) then
      write(msg,'(A,I0,A,I0,A)') '[ERROR] Mesh type: ', nturn, ' of the ', nb, &
        ' blocks of the 2-D mesh are axisymmetric wedges (k = 1 plane turned about the x axis), the others are not'
      call mpi_abort_all ( trim(msg) )
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


  subroutine Compute_Metric_Tensor( N1, N2, N3, N4, N5, N6, N7, N8, M, dl, vol )
    implicit none
    type(FUSS_vector_3D_type), intent(in)  :: N1, N2, N3, N4, N5, N6, N7, N8
    type(FUSS_tensor_3D_type), intent(out) :: M
    type(FUSS_vector_3D_type), intent(out) :: dl
    real(R8), intent(out)             :: vol
    ! Local
    integer  :: h
    real(R8) :: det, A(3,3), cofactor(3,3)
    real(R8) :: vx(8), vy(8), vz(8)

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

    if ( abs(det) == 0d0 ) then
      ! Degenerate cell (zero volume): a mesh error, or ghost cells extrapolated from a
      ! wall-thin first cell whose coordinates were read back in single precision.
      ! error stop, not mpi_abort_all: this can run inside OpenMP threads (MPI is FUNNELED).
      write(*,'(A)') ' [ERROR] Metric tensor det=0: degenerate cell. Check the mesh, or the precision of the restart file'
      write(*,'(A,3(1X,ES14.6))') ' [ERROR]   corner node N1:', N1 % c
      write(*,'(A,3(1X,ES14.6))') ' [ERROR]   corner node N8:', N8 % c
      error stop

    else
      
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

    ! Cell volume computation
    vx(1) = N1%c(1)
    vy(1) = N1%c(2)
    vz(1) = N1%c(3)

    vx(2) = N5%c(1)
    vy(2) = N5%c(2)
    vz(2) = N5%c(3)

    vx(3) = N3%c(1)
    vy(3) = N3%c(2)
    vz(3) = N3%c(3)

    vx(4) = N7%c(1)
    vy(4) = N7%c(2)
    vz(4) = N7%c(3)
  
    vx(5) = N2%c(1)
    vy(5) = N2%c(2)
    vz(5) = N2%c(3)

    vx(6) = N6%c(1)
    vy(6) = N6%c(2)
    vz(6) = N6%c(3)

    vx(7) = N4%c(1)
    vy(7) = N4%c(2)
    vz(7) = N4%c(3)
  
    vx(8) = N8%c(1)
    vy(8) = N8%c(2)
    vz(8) = N8%c(3)
  
    vol = tvol(vx,vy,vz,1,2,3,5) + tvol(vx,vy,vz,2,4,3,8) &
        + tvol(vx,vy,vz,5,8,6,2) + tvol(vx,vy,vz,5,7,8,3) &
        + tvol(vx,vy,vz,5,8,2,3)
  
    if( vol <= 0d0 ) then
      write(*,'(A90)') ' [ERROR] Negative volume'
      stop
    endif

    contains
  
      pure function tvol(vx, vy, vz, i1, i2, i3, i4) result(volume)
        implicit none
        real(R8), intent(in) :: vx(8), vy(8), vz(8)
        integer, intent(in)  :: i1, i2, i3, i4
        real(R8)             :: volume
  
        volume = abs(((vx(i2)-vx(i1))* &
          ((vy(i3)-vy(i1))*(vz(i4)-vz(i1))-(vy(i4)-vy(i1))*(vz(i3)-vz(i1)))+ &
                                  (vy(i2)-vy(i1))* &
          ((vx(i4)-vx(i1))*(vz(i3)-vz(i1))-(vx(i3)-vx(i1))*(vz(i4)-vz(i1)))+ &
                                  (vz(i2)-vz(i1))* &
          ((vx(i3)-vx(i1))*(vy(i4)-vy(i1))-(vx(i4)-vx(i1))*(vy(i3)-vy(i1)))) &
          /6.d0)
  
      end function tvol
  
  end subroutine Compute_Metric_Tensor


  subroutine BC_Extrapolate_Metrics( Im, Jm, Km, Fm, blk, Mg, dlg, volg )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout) :: blk
    type(FUSS_tensor_3D_type), intent(out)    :: Mg(2)
    type(FUSS_vector_3D_type), intent(out)    :: dlg(2)
    real(R8), intent(out)                :: volg(2)
    ! Local
    integer :: g
    type(FUSS_vector_3D_type), dimension(-gc:gc):: N1, N2, N3, N4, N5, N6, N7, N8

    select case (Fm)
      case(1)
        do g = 0, gc
          N1(g)%c = blk % node (Im+g-1,Jm-1,Km-1) % c
          N2(g)%c = blk % node (Im+g-1,Jm-1,Km  ) % c
          N3(g)%c = blk % node (Im+g-1,Jm  ,Km-1) % c
          N4(g)%c = blk % node (Im+g-1,Jm  ,Km  ) % c
        enddo
        do g = -1, -gc, -1
          N1(g)%c = 3d0*N1(g+1)%c - 3d0*N1(g+2)%c + N1(g+3)%c
          N2(g)%c = 3d0*N2(g+1)%c - 3d0*N2(g+2)%c + N2(g+3)%c
          N3(g)%c = 3d0*N3(g+1)%c - 3d0*N3(g+2)%c + N3(g+3)%c
          N4(g)%c = 3d0*N4(g+1)%c - 3d0*N4(g+2)%c + N4(g+3)%c
          call Compute_Metric_Tensor ( N1(g), N2(g), N3(g), N4(g), N1(g+1), N2(g+1), N3(g+1), N4(g+1), &
                                       blk % M(Im+g, Jm, Km), blk % dl(Im+g, Jm, Km), blk % vol(Im+g, Jm, Km))
          Mg(-g)   = blk % M(Im+g, Jm, Km)
          dlg(-g)  = blk % dl(Im+g, Jm, Km)
          volg(-g) = blk % vol(Im+g, Jm, Km)
        enddo
        
      case(2)
        do g = 0, -gc, -1
          N5(g)%c = blk % node (Im+g,Jm-1,Km-1) % c
          N6(g)%c = blk % node (Im+g,Jm-1,Km  ) % c
          N7(g)%c = blk % node (Im+g,Jm  ,Km-1) % c
          N8(g)%c = blk % node (Im+g,Jm  ,Km  ) % c
        enddo
        do g = 1, gc
          N5(g)%c = 3d0*N5(g-1)%c - 3d0*N5(g-2)%c + N5(g-3)%c
          N6(g)%c = 3d0*N6(g-1)%c - 3d0*N6(g-2)%c + N6(g-3)%c
          N7(g)%c = 3d0*N7(g-1)%c - 3d0*N7(g-2)%c + N7(g-3)%c
          N8(g)%c = 3d0*N8(g-1)%c - 3d0*N8(g-2)%c + N8(g-3)%c
          call Compute_Metric_Tensor ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                       blk % M(Im+g, Jm, Km), blk % dl(Im+g, Jm, Km), blk % vol(Im+g, Jm, Km))
          Mg(g)   = blk % M(Im+g, Jm, Km)
          dlg(g)  = blk % dl(Im+g, Jm, Km)
          volg(g) = blk % vol(Im+g, Jm, Km)
        enddo

      case(3)
        if ( blk % dim(2) > 2 ) then
          do g = 0, gc
            N1(g)%c = blk % node (Im-1,Jm+g-1,Km-1) % c
            N2(g)%c = blk % node (Im-1,Jm+g-1,Km  ) % c
            N5(g)%c = blk % node (Im  ,Jm+g-1,Km-1) % c
            N6(g)%c = blk % node (Im  ,Jm+g-1,Km  ) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 3d0*N1(g+1)%c - 3d0*N1(g+2)%c + N1(g+3)%c
            N2(g)%c = 3d0*N2(g+1)%c - 3d0*N2(g+2)%c + N2(g+3)%c
            N5(g)%c = 3d0*N5(g+1)%c - 3d0*N5(g+2)%c + N5(g+3)%c
            N6(g)%c = 3d0*N6(g+1)%c - 3d0*N6(g+2)%c + N6(g+3)%c
            call Compute_Metric_Tensor ( N1(g), N2(g), N1(g+1), N2(g+1), N5(g), N6(g), N5(g+1), N6(g+1), &
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(-g)   = blk % M(Im, Jm+g, Km)
            dlg(-g)  = blk % dl(Im, Jm+g, Km)
            volg(-g) = blk % vol(Im, Jm+g, Km)
          enddo
        else
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
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(-g)   = blk % M(Im, Jm+g, Km)
            dlg(-g)  = blk % dl(Im, Jm+g, Km)
            volg(-g) = blk % vol(Im, Jm+g, Km)
          enddo
        endif

      case(4)
        if ( blk % dim(2) > 2 ) then
          do g = 0, -gc, -1
            N3(g)%c = blk % node (Im-1,Jm+g,Km-1) % c
            N4(g)%c = blk % node (Im-1,Jm+g,Km  ) % c
            N7(g)%c = blk % node (Im  ,Jm+g,Km-1) % c
            N8(g)%c = blk % node (Im  ,Jm+g,Km  ) % c
          enddo
          do g = 1, gc
            N3(g)%c = 3d0*N3(g-1)%c - 3d0*N3(g-2)%c + N3(g-3)%c
            N4(g)%c = 3d0*N4(g-1)%c - 3d0*N4(g-2)%c + N4(g-3)%c
            N7(g)%c = 3d0*N7(g-1)%c - 3d0*N7(g-2)%c + N7(g-3)%c
            N8(g)%c = 3d0*N8(g-1)%c - 3d0*N8(g-2)%c + N8(g-3)%c
            call Compute_Metric_Tensor ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(g)   = blk % M(Im, Jm+g, Km)
            dlg(g)  = blk % dl(Im, Jm+g, Km)
            volg(g) = blk % vol(Im, Jm+g, Km)
          enddo
        else
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
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(g)   = blk % M(Im, Jm+g, Km)
            dlg(g)  = blk % dl(Im, Jm+g, Km)
            volg(g) = blk % vol(Im, Jm+g, Km)
          enddo
        endif

      case(5)
        if ( ndir==2 .and. delthe==0d0 .or. ndir==1 ) then
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
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo
        elseif (ndir==2 .and. delthe/=0d0) then
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
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo         
        elseif (ndir==3) then
          do g = 0, gc
            N1(g)%c = blk % node (Im-1,Jm-1,Km+g-1) % c
            N3(g)%c = blk % node (Im-1,Jm  ,Km+g-1) % c
            N5(g)%c = blk % node (Im  ,Jm-1,Km+g-1) % c
            N7(g)%c = blk % node (Im  ,Jm  ,Km+g-1) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 3d0*N1(g+1)%c - 3d0*N1(g+2)%c + N1(g+3)%c
            N3(g)%c = 3d0*N3(g+1)%c - 3d0*N3(g+2)%c + N3(g+3)%c
            N5(g)%c = 3d0*N5(g+1)%c - 3d0*N5(g+2)%c + N5(g+3)%c
            N7(g)%c = 3d0*N7(g+1)%c - 3d0*N7(g+2)%c + N7(g+3)%c
            call Compute_Metric_Tensor ( N1(g), N1(g+1), N3(g), N3(g+1), N5(g), N5(g+1), N7(g), N7(g+1), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo
        endif

      case(6)
        if ( ndir==2 .and. delthe==0d0 .or. ndir==1 ) then
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
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo
        elseif (ndir==2 .and. delthe/=0d0) then
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
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo         
        elseif (ndir==3) then
          do g = 0, -gc, -1
            N2(g)%c = blk % node (Im-1,Jm-1,Km+g) % c
            N4(g)%c = blk % node (Im-1,Jm  ,Km+g) % c
            N6(g)%c = blk % node (Im  ,Jm-1,Km+g) % c
            N8(g)%c = blk % node (Im  ,Jm  ,Km+g) % c
          enddo
          do g = 1, gc
            N2(g)%c = 3d0*N2(g-1)%c - 3d0*N2(g-2)%c + N2(g-3)%c
            N4(g)%c = 3d0*N4(g-1)%c - 3d0*N4(g-2)%c + N4(g-3)%c
            N6(g)%c = 3d0*N6(g-1)%c - 3d0*N6(g-2)%c + N6(g-3)%c
            N8(g)%c = 3d0*N8(g-1)%c - 3d0*N8(g-2)%c + N8(g-3)%c
            call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo
        endif

      end select

  end subroutine BC_Extrapolate_Metrics


  subroutine BC_Symmetry_Metrics( Im, Jm, Km, Fm, blk, Mg, dlg, volg )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(FUSS_block_type), intent(inout)   :: blk
    type(FUSS_tensor_3D_type), intent(out) :: Mg(2)
    type(FUSS_vector_3D_type), intent(out) :: dlg(2)
    real(R8), intent(out)                  :: volg(2)
    ! Local
    integer :: g

    do g = 1, gc
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
    enddo

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