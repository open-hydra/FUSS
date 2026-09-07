!>@brief The ONLY file in FUSS that names both FUSS and MORPH types.
!>
!> MORPH is a standalone geometry component (src/lib/morph) that must never
!> depend on FUSS, so that extracting it into its own hydra solver is a
!> directory move rather than a rewrite (plan 09). This adapter is the seam.
!> Keeping it to one file is checked mechanically by test/run_morph_guards.sh.
!>
!> WHY THERE IS A COPY HERE
!> ------------------------
!> Plan 09 assumed that giving MORPH structurally identical types would make the
!> adapter a copy-free pass-through. It does not: Fortran treats structurally
!> identical derived types as DISTINCT types, so a FUSS_vector_3D_type array
!> cannot be passed to a morph_vec3_t dummy. The options were to alias the
!> storage with c_f_pointer, or to copy. This copies, because:
!>
!>   * it is obviously correct, with no aliasing or contiguity assumptions;
!>   * it is O(nodes) per call, the same order as the metric computation it
!>     feeds, so it cannot be worse than a constant factor;
!>   * plan 05 section 4.2 says to profile the per-step metric refresh before
!>     optimising it, and this is part of that same measurement.
!>
!> If profiling later shows the copy matters, the fix is local to this file
!> (alias via c_f_pointer) and changes no interface.
!>
!> BIT-IDENTITY
!> ------------
!> Only cells 1..dim are copied back, which is exactly the range the original
!> Setup_Metrics filled -- the ghost entries of blk%M / blk%dl / blk%vol were
!> never written by it (ghost-cell metrics live separately in bc%Mg/dlg/volg).
!> Leaving them untouched is what lets the zero-motion regression reproduce the
!> Phase 0 results bit-identically.
module FUSS_Adapter_Morph
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Base_Types_m
  use FUSS_Advanced_Types_m
  use Morph_API

  implicit none
  private

  public :: Adapter_Mesh_Type
  public :: Adapter_Block_Metrics
  public :: Adapter_Free

  !> One geometry bundle per block, kept across calls so the per-step path does
  !> not churn the allocator. Indexed by block number.
  type(morph_geom_t), allocatable, save :: geom_(:)

contains

  !> Classify the mesh (1-D / 2-D / 2-D-axisymmetric / 3-D) using block 1.
  !> Returns what FUSS's globals `ndir` and `delthe` used to be set to directly.
  subroutine Adapter_Mesh_Type ( blk, ndir_out, delthe_out )
    type(FUSS_block_type), intent(in)  :: blk
    integer(I4),           intent(out) :: ndir_out
    real(R8),              intent(out) :: delthe_out
    ! Local
    type(morph_vec3_t), allocatable :: node(:,:,:)

    call pack_nodes ( blk, node )
    call Morph_Mesh_Type ( node, blk%dim, ndir_out, delthe_out )

  end subroutine Adapter_Mesh_Type


  !> Compute cell volume, metric tensor, cell lengths, face normals and areas
  !> for one block, writing them into the block's own arrays.
  subroutine Adapter_Block_Metrics ( blk, nb, b, ok, message, ci, cj, ck )
    type(FUSS_block_type), intent(inout) :: blk
    integer(I4),           intent(in)    :: nb        !> number of blocks
    integer(I4),           intent(in)    :: b         !> this block's index
    logical,               intent(out)   :: ok
    character(len=*),      intent(out)   :: message
    integer(I4),           intent(out)   :: ci, cj, ck
    ! Local
    type(morph_vec3_t), allocatable :: node(:,:,:)
    type(morph_status_t) :: status
    integer(I4) :: i, j, k, d

    if ( .not. allocated(geom_) ) allocate( geom_(1:nb) )

    call pack_nodes ( blk, node )
    call Morph_Metrics_Block_Wrapper ( node, blk%dim, geom_(b), status )

    ok = status%ok()
    message = status%message
    ci = status%i; cj = status%j; ck = status%k

    ! ---- Copy back, interior cells only (see the header note on bit-identity).
    do k = 1, blk%dim(3)
    do j = 1, blk%dim(2)
    do i = 1, blk%dim(1)
      blk%vol(i,j,k)   = geom_(b)%vol(i,j,k)
      blk%M(i,j,k)%c   = geom_(b)%M(i,j,k)%c
      blk%dl(i,j,k)%c  = geom_(b)%dl(i,j,k)%c
    enddo; enddo; enddo

    ! Face normals and areas. Bounds match FUSS's allocation exactly:
    !   dir(1)%f(0:ni, nj, nk)  dir(2)%f(ni, 0:nj, nk)  dir(3)%f(ni, nj, 0:nk)
    do d = 1, 3
      do k = lbound(blk%dir(d)%f,3), ubound(blk%dir(d)%f,3)
      do j = lbound(blk%dir(d)%f,2), ubound(blk%dir(d)%f,2)
      do i = lbound(blk%dir(d)%f,1), ubound(blk%dir(d)%f,1)
        blk%dir(d)%f(i,j,k)%N = geom_(b)%dir(d)%f(i,j,k)%n
        blk%dir(d)%f(i,j,k)%A = geom_(b)%dir(d)%f(i,j,k)%A
      enddo; enddo; enddo
    enddo

  end subroutine Adapter_Block_Metrics


  subroutine Adapter_Free ()
    if ( allocated(geom_) ) deallocate( geom_ )
  end subroutine Adapter_Free


  ! ---------------------------------------------------------------------------

  !> FUSS_vector_3D_type -> morph_vec3_t. See the header for why this copies.
  subroutine pack_nodes ( blk, node )
    type(FUSS_block_type),           intent(in)  :: blk
    type(morph_vec3_t), allocatable, intent(out) :: node(:,:,:)
    ! Local
    integer(I4) :: i, j, k

    allocate( node( 0:blk%dim(1), 0:blk%dim(2), 0:blk%dim(3) ) )

    do k = 0, blk%dim(3)
    do j = 0, blk%dim(2)
    do i = 0, blk%dim(1)
      node(i,j,k)%c = blk%node(i,j,k)%c
    enddo; enddo; enddo

  end subroutine pack_nodes


  !> Thin indirection so the `use` of Morph_Metrics_Block stays inside the
  !> adapter rather than leaking the non-facade module name into FUSS.
  subroutine Morph_Metrics_Block_Wrapper ( node, dim, geom, status )
    use Morph_API, only: Morph_Init
    type(morph_vec3_t),   intent(in)    :: node(0:,0:,0:)
    integer(I4),          intent(in)    :: dim(3)
    type(morph_geom_t),   intent(inout) :: geom
    type(morph_status_t), intent(out)   :: status

    call Morph_Init ( node, dim, geom, status )

  end subroutine Morph_Metrics_Block_Wrapper

end module FUSS_Adapter_Morph
