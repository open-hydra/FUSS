!> Periodic faces (201) on ours' moving mesh (plan 14 S9d, Part B).
!>
!> A periodic pair (Lib_Periodic, upstream) makes the domain one period of the ring, or row, of its copies: face f of
!> block b takes as ghost cells the cells inside face fs of block bs, moved by the rigid motion x -> R x + t that
!> Setup_Periodic fits once, at set-up. On a moving mesh the copies still fit together, and the ghost geometry that
!> Update_Mesh rebuilds every step with that motion (Setup_BC_Metrics, Setup_Face_Metrics) is still the geometry of the
!> neighbour copy, exactly while the motion takes image nodes to image nodes: u(x_f) = R u(x_fs) for every node x_f of
!> face f and its image x_fs on face fs. Then x_f + u(x_f) = R (x_fs + u(x_fs)) + t: the pair keeps its set-up motion,
!> and the swept volumes of the two faces are images, so the remap carries the same volume through both. A motion that
!> breaks it describes no periodic domain, even when the two faces stay images by another rigid motion (a translation d
!> across the axis of a rotational pair: the faces are images by (R, t + (I - R) d), but the copies of the sector would
!> move by R^m d, not by d, and the swept volumes of the two faces differ); refitting (R, t) would hide it.
!>  - Check_Periodic_Translation (set-up, law = translation): R vel = vel for every pair, to 2 MISFIT_MAX of |vel|.
!>  - Check_Periodic_Faces (every step, every law, after exchange_node_planes): the misfit of every periodic record of
!>    this rank's blocks against the pair's motion of the set-up, by Setup_Periodic's own measure (Record_Points: the
!>    centres of the two faces and the points half a cell across them); it stops on GROWTH: beyond the pair's set-up
!>    misfit plus the set-up tolerance (MISFIT_MAX half cells). The ranks stop together.
!> The ablating wall (law = nodeshift) keeps the images through Setup_Recession_Seams (Mod_Recession).
!> The bound 2 MISFIT_MAX of the set-up checks: an invariant direction or tensor seen through a rotation fitted within
!> the pair's tolerance (an angle error of at most MISFIT_MAX half cell over a face of at least a half cell) moves by at
!> most that fraction.
module FUSS_Mod_Periodic_Checks
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Check_Periodic_Translation, Check_Periodic_Faces

contains

  !> Set-up, law = translation, on every rank (the pairs and the velocity are the same on every rank): a velocity that
  !> the rotation of a pair does not leave alone moves the copies of the domain apart. Stops every rank.
  subroutine Check_Periodic_Translation ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_mesh_motion
    use FUSS_Lib_Periodic,   only: MISFIT_MAX
    use FUSS_Mod_MPI,        only: mpi_stop_all
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    integer  :: g
    real(R8) :: v(3), dv
    character(len=600) :: msg

    if ( trim(obj_mesh_motion%law) /= 'translation' ) return
    if ( .not. allocated(domain%per) ) return
    v = obj_mesh_motion%vel
    if ( norm2(v) == 0d0 ) return
    do g = 1, size(domain%per)
      dv = norm2 ( matmul ( domain%per(g)%R, v ) - v ) / norm2(v)
      if ( dv > 2d0 * MISFIT_MAX ) then
        write(msg,'(A,4(I0,A),ES10.3,A)') '[ERROR] periodic faces (201): block ', domain%per(g)%b, ' face ', &
          domain%per(g)%f, ' <- block ', domain%per(g)%bs, ' face ', domain%per(g)%fs, &
          ': law = translation moves the mesh across the axis of the pair (|R vel - vel| / |vel| = ', dv, &
          '): the copies of a periodic domain would move by R vel, not by vel. Translate along the axis of the '// &
          'pair, or mesh the full ring'
        call mpi_stop_all ( trim(msg) )
      end if
    end do

  end subroutine Check_Periodic_Translation


  !> Every step of a moving mesh, after exchange_node_planes (the node planes of the other face of a pair come to this
  !> rank with the planes of that face), before the ghost metrics are rebuilt: the misfit of every periodic record of
  !> this rank's blocks against the pair's motion of the set-up, as Setup_Periodic measures it (Record_Points, nodes
  !> only: the centres of face f and face fs and the points half a cell across them) and with its tolerance
  !> (MISFIT_MAX of the pair's half cell of the set-up), on growth: a record stops when its misfit passes the pair's
  !> set-up misfit (per(g)%misfit, Setup_Periodic: the pair's maximum over its records, so every record gets the worst
  !> record's slack, intended) plus that tolerance; a pair accepted near the bound keeps the full tolerance for its
  !> motion. The worst record beyond it stops every rank (mpi_stop_if_any: one reduction per step, only with periodic
  !> pairs). A restart refits the pair (Setup_Periodic), so each restart re-baselines the misfit (open item
  !> O-PERIODIC-RESTART). Like the set-up measure it does not see an in-plane node error whose sign alternates between
  !> adjacent node lines (it cancels in every face centre and keeps every face in its plane); no smooth motion makes one.
  subroutine Check_Periodic_Faces ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_mesh_motion
    use FUSS_Lib_Periodic,   only: Periodic_Move, Record_Points, MISFIT_MAX
    use FUSS_Mod_MPI,        only: is_local_block, mpi_stop_if_any
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    integer  :: i, g, n, iw
    real(R8) :: p(3,2), q(3,2), h, dev, worst
    character(len=800) :: detail

    if ( .not. allocated(domain%per) ) return
    if ( size(domain%per) == 0 ) return
    worst = 0d0 ; iw = 0
    do i = 1, domain%nbound
      g = domain%bc(i)%per
      if ( g <= 0 ) cycle
      if ( .not. is_local_block( domain%bc(i)%b ) ) cycle
      call Record_Points ( domain, i, p, q, h )
      dev = 0d0
      do n = 1, 2
        dev = max ( dev, norm2 ( Periodic_Move ( domain%per(g), p(:,n) ) - q(:,n) ) )
      enddo
      if ( dev > domain%per(g)%misfit + MISFIT_MAX * domain%per(g)%hmin .and. dev / domain%per(g)%hmin > worst ) then
        worst = dev / domain%per(g)%hmin
        iw = i
      end if
    enddo
    detail = ''
    if ( iw > 0 ) then
      g = domain%bc(iw)%per
      write(detail,'(A,5(I0,A),ES10.3,A,ES10.3,A,ES10.3,A)') ' [ERROR] periodic faces (201): bc.txt record ', iw, &
        ', block ', domain%bc(iw)%b, ' face ', domain%bc(iw)%f, ' <- block ', domain%bc(iw)%bs, ' face ', &
        domain%bc(iw)%fs, ': after the mesh motion (law = '//trim(obj_mesh_motion%law)//') the faces lie ', &
        worst * domain%per(g)%hmin, ' m from images by the motion of the set-up (set-up misfit ', domain%per(g)%misfit, &
        ' m, tolerance ', MISFIT_MAX * domain%per(g)%hmin, ' m on top of it): the motion does not take image nodes to '// &
        'image nodes, '// &
        'so the problem is not periodic any more (translation: along the axis of the pair; prescribed: a '// &
        'displacement whose values at a node and at its image are images)'
    end if
    call mpi_stop_if_any ( iw > 0, 'periodic faces (201) no longer images after the mesh motion (see above)', detail )

  end subroutine Check_Periodic_Faces

end module FUSS_Mod_Periodic_Checks
