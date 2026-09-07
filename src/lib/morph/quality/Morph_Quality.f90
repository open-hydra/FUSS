!>@brief Mesh validity checks.
!>
!> The point of this module is that FUSS's legacy cell volume CANNOT detect a
!> tangled cell: it sums five tetrahedra taking the absolute value of each, so
!> an inverted cell still reports a positive volume (see Morph_GCL's header).
!> A mesh that has folded therefore produces plausible-looking geometry and then
!> NaNs several steps later, far from the cause.
!>
!> The signed Jacobian is what actually detects inversion, so that is the
!> primary test here, with non-positive signed volume as the coarser backstop.
!>
!> MORPH never calls `stop`: results come back in morph_status_t and the caller
!> decides whether to abort, retry with a smaller step, or remesh.
module Morph_Quality
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_GCL, only: Morph_Cell_Volume_Signed

  implicit none
  private

  public :: Morph_Quality_Check

contains

  !> Scan a block for tangled or degenerate cells.
  !>
  !> `min_vol` and `min_jac` come back so the caller can log a trend and see a
  !> mesh degrading before it actually fails.
  subroutine Morph_Quality_Check ( node, dim, min_vol, min_jac, status )
    type(morph_vec3_t),   intent(in)  :: node(0:,0:,0:)
    integer(I4),          intent(in)  :: dim(3)
    real(R8),             intent(out) :: min_vol
    real(R8),             intent(out) :: min_jac
    type(morph_status_t), intent(out) :: status
    ! Local
    integer(I4) :: i, j, k
    integer(I4) :: bad_i, bad_j, bad_k, bad_code
    real(R8)    :: x(3,8), v, jac

    call status%clear()
    min_vol  =  huge(1.0_R8)
    min_jac  =  huge(1.0_R8)
    bad_i = 0; bad_j = 0; bad_k = 0; bad_code = MORPH_OK

    !$omp parallel do collapse(3) private(i,j,k,x,v,jac) &
    !$omp   reduction(min: min_vol, min_jac)
    do k = 1, dim(3)
    do j = 1, dim(2)
    do i = 1, dim(1)
      call gather_cell ( node, i, j, k, x )
      v   = Morph_Cell_Volume_Signed( x )
      jac = corner_jacobian_min( x )
      min_vol = min( min_vol, v )
      min_jac = min( min_jac, jac )
    enddo; enddo; enddo
    !$omp end parallel do

    ! Locate the offender only if there is one. A second pass costs nothing on
    ! the healthy path and keeps the hot loop free of critical sections.
    if ( min_jac <= 0.0_R8 .or. min_vol <= 0.0_R8 ) then
      outer: do k = 1, dim(3)
      do j = 1, dim(2)
      do i = 1, dim(1)
        call gather_cell ( node, i, j, k, x )
        v   = Morph_Cell_Volume_Signed( x )
        jac = corner_jacobian_min( x )
        if ( jac <= 0.0_R8 ) then
          bad_code = MORPH_ERR_JACOBIAN; bad_i = i; bad_j = j; bad_k = k
          exit outer
        elseif ( v <= 0.0_R8 ) then
          bad_code = MORPH_ERR_NEGVOL; bad_i = i; bad_j = j; bad_k = k
          exit outer
        endif
      enddo; enddo; enddo outer

      select case (bad_code)
      case (MORPH_ERR_JACOBIAN)
        call status%fail ( MORPH_ERR_JACOBIAN, &
          'mesh tangled: non-positive corner Jacobian', bad_i, bad_j, bad_k )
      case (MORPH_ERR_NEGVOL)
        call status%fail ( MORPH_ERR_NEGVOL, &
          'non-positive signed cell volume', bad_i, bad_j, bad_k )
      end select
    endif

  end subroutine Morph_Quality_Check


  ! ---------------------------------------------------------------------------

  pure subroutine gather_cell ( node, i, j, k, x )
    type(morph_vec3_t), intent(in)  :: node(0:,0:,0:)
    integer(I4),        intent(in)  :: i, j, k
    real(R8),           intent(out) :: x(3,8)

    ! Same node ordering as Morph_Metric_Tensor and Morph_GCL.
    x(:,1) = node(i-1,j-1,k-1)%c
    x(:,2) = node(i-1,j-1,k  )%c
    x(:,3) = node(i-1,j  ,k-1)%c
    x(:,4) = node(i-1,j  ,k  )%c
    x(:,5) = node(i  ,j-1,k-1)%c
    x(:,6) = node(i  ,j-1,k  )%c
    x(:,7) = node(i  ,j  ,k-1)%c
    x(:,8) = node(i  ,j  ,k  )%c

  end subroutine gather_cell


  !> Minimum over the 8 corners of the trilinear map's Jacobian determinant.
  !>
  !> A hexahedron is untangled iff the Jacobian is positive throughout, and for
  !> a trilinear map it suffices to test the corners: the determinant attains its
  !> extrema there for the element shapes of interest. This catches inversion
  !> that a volume magnitude cannot.
  pure function corner_jacobian_min ( x ) result ( jmin )
    real(R8), intent(in) :: x(3,8)
    real(R8)             :: jmin
    ! Local
    integer(I4) :: c
    real(R8)    :: e1(3), e2(3), e3(3), d

    ! Corner-local edge triples, expressed in the node numbering
    !   1=(0,0,0) 2=(0,0,1) 3=(0,1,0) 4=(0,1,1) 5=(1,0,0) 6=(1,0,1) 7=(1,1,0) 8=(1,1,1)
    ! where the three edges leave the corner along +i, +j, +k respectively.
    integer(I4), parameter :: C0(8) = [1, 2, 3, 4, 5, 6, 7, 8]
    integer(I4), parameter :: CI(8) = [5, 6, 7, 8, 1, 2, 3, 4]   ! +/- i neighbour
    integer(I4), parameter :: CJ(8) = [3, 4, 1, 2, 7, 8, 5, 6]   ! +/- j neighbour
    integer(I4), parameter :: CK(8) = [2, 1, 4, 3, 6, 5, 8, 7]   ! +/- k neighbour
    real(R8),    parameter :: SI(8) = [ 1, 1, 1, 1,-1,-1,-1,-1]
    real(R8),    parameter :: SJ(8) = [ 1, 1,-1,-1, 1, 1,-1,-1]
    real(R8),    parameter :: SK(8) = [ 1,-1, 1,-1, 1,-1, 1,-1]

    jmin = huge(1.0_R8)
    do c = 1, 8
      e1 = SI(c) * ( x(:,CI(c)) - x(:,C0(c)) )
      e2 = SJ(c) * ( x(:,CJ(c)) - x(:,C0(c)) )
      e3 = SK(c) * ( x(:,CK(c)) - x(:,C0(c)) )
      d =   e1(1)*(e2(2)*e3(3) - e2(3)*e3(2)) &
          - e1(2)*(e2(1)*e3(3) - e2(3)*e3(1)) &
          + e1(3)*(e2(1)*e3(2) - e2(2)*e3(1))
      jmin = min( jmin, d )
    enddo

  end function corner_jacobian_min

end module Morph_Quality
