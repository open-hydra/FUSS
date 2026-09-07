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
    use FUSS_Lib_Metrics, only: delthe, BC_Connect_Metrics, BC_Symmetry_Metrics, &
                                BC_Extrapolate_Metrics
    use FUSS_Adapter_Morph, only: Adapter_Mesh_Type, Adapter_Block_Metrics
    use FUSS_Mod_MPI, only: is_local_block
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i
    integer :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
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
        stop
      endif
    enddo

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

  end subroutine Setup_Metrics

end module FUSS_Mod_Metrics
