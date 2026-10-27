!>@brief MPI infrastructure module for FUSS.
!> Provides MPI environment management, block partitioning, and
!> global reduction wrappers. When compiled without USE_MPI, all
!> routines are no-ops or identity operations (serial fallback).
module FUSS_Mod_MPI
#ifdef USE_MPI
  use mpi
#endif
  use iso_fortran_env, only: R8 => real64, I4 => int32, output_unit, error_unit

  implicit none
  private

  ! --- Public state ---
  integer, public :: mpi_rank_ = 0        !< This process rank (0 in serial)
  integer, public :: mpi_size_ = 1        !< Total number of MPI processes
  logical, public :: mpi_is_root = .true.  !< True on rank 0

  ! --- Block-to-rank mapping ---
  integer, allocatable, public :: block_owner(:)     !< block_owner(b) = rank owning block b
  integer, allocatable, public :: local_block_ids(:) !< Global block IDs owned by this rank
  integer, public :: n_local_blocks = 0

  ! --- Public procedures ---
  public :: mpi_init_env, mpi_finalize_env
  public :: is_local_block
  public :: partition_blocks
  public :: mpi_allreduce_sum_r8, mpi_allreduce_min_r8, mpi_allreduce_max_r8
  public :: mpi_gather_r8
  public :: mpi_allreduce_sum_r8_array
  public :: mpi_reduce_sum_r8, mpi_reduce_sum_r8_array
  public :: mpi_allreduce_norm2
  public :: mpi_bcast_logical, mpi_bcast_integer, mpi_bcast_chars
  public :: mpi_abort_all, mpi_stop_all, mpi_stop_if_any, stop_from_thread
  public :: check_mpi_error

contains


  !> Initialize MPI environment with thread support. Call early in main program.
  !> Uses MPI_THREAD_FUNNELED: only master thread makes MPI calls.
  subroutine mpi_init_env()
#ifdef USE_MPI
    integer :: ierr, provided
    logical :: already
    
    call MPI_Initialized(already, ierr)
    if (.not. already) then
      call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, ierr)
      call check_mpi_error(ierr)
      if (provided < MPI_THREAD_FUNNELED) then
        write(*,'(A)') ' ERROR: MPI does not support the required threading level (MPI_THREAD_FUNNELED)'
        call flush_output_units()
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      end if
    end if
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank_, ierr)
    call check_mpi_error(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size_, ierr)
    call check_mpi_error(ierr)
    mpi_is_root = (mpi_rank_ == 0)
#else
    mpi_rank_ = 0
    mpi_size_ = 1
    mpi_is_root = .true.
#endif
  end subroutine mpi_init_env


  !> Finalize MPI environment. Call at end of main program.
  subroutine mpi_finalize_env()
#ifdef USE_MPI
    integer :: ierr
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
#endif
  end subroutine mpi_finalize_env


  !> Check if block b is owned by this rank.
  logical function is_local_block(b)
    integer, intent(in) :: b

    if (.not. allocated(block_owner)) then
      is_local_block = .true.  ! Serial mode or not yet partitioned
      return
    end if
    is_local_block = (block_owner(b) == mpi_rank_)
  end function is_local_block


  !> Partition nb blocks across MPI ranks: largest block first, each to the least loaded rank (LPT).
  !> blk_ncells(b) = total number of cells in block b.
  subroutine partition_blocks(nb, blk_ncells)
    integer, intent(in) :: nb
    integer, intent(in) :: blk_ncells(nb)
    ! Local
    integer :: b, i, j, r, nloc, tmp
    integer, allocatable :: rank_load(:), order(:)

    if (allocated(block_owner)) deallocate(block_owner)
    if (allocated(rank_load)) deallocate(rank_load)
    allocate(block_owner(nb))
    allocate(rank_load(0:mpi_size_-1))
    allocate(order(nb))
    rank_load = 0

    ! Blocks by descending cell count (insertion sort; ties keep mesh order)
    do b = 1, nb
      order(b) = b
    end do
    do i = 2, nb
      tmp = order(i)
      j = i - 1
      do while (j >= 1)
        if (blk_ncells(order(j)) >= blk_ncells(tmp)) exit
        order(j+1) = order(j)
        j = j - 1
      end do
      order(j+1) = tmp
    end do

    do i = 1, nb
      b = order(i)
      r = minloc(rank_load, dim=1) - 1   ! rank with minimum load (0-indexed)
      block_owner(b) = r
      rank_load(r) = rank_load(r) + blk_ncells(b)
    end do
    deallocate(order)

    ! Build local block list
    nloc = count(block_owner == mpi_rank_)
    n_local_blocks = nloc
    if (allocated(local_block_ids)) deallocate(local_block_ids)
    allocate(local_block_ids(nloc))
    nloc = 0
    do b = 1, nb
      if (block_owner(b) == mpi_rank_) then
        nloc = nloc + 1
        local_block_ids(nloc) = b
      end if
    end do

    call report_partition_balance(nb, blk_ncells, rank_load)

    deallocate(rank_load)
  end subroutine partition_blocks


  !> Report how evenly the blocks landed on the ranks. In a coupled run FUSS
  !> partitions the solid independently of MOSE's fluid partition, so this is
  !> the only place the solid's balance becomes visible.
  subroutine report_partition_balance(nb, blk_ncells, rank_load)
    integer, intent(in) :: nb
    integer, intent(in) :: blk_ncells(nb)
    integer, intent(in) :: rank_load(0:)
    ! Local
    real(R8) :: ideal, eff

    if (mpi_rank_ /= 0) return
    if (mpi_size_ <= 1) return

    ideal = real(sum(blk_ncells), R8) / real(mpi_size_, R8)
    eff   = ideal / real(maxval(rank_load), R8) * 100.0_R8

    write(*,'(A,I0,A,I0,A,F5.1,A)') '  FUSS MPI partition: ', nb, ' blocks over ', &
      mpi_size_, ' ranks, balance ', eff, '% of ideal'

    if (mpi_size_ > nb) then
      write(*,'(A,I0,A,I0,A)') '  WARNING: blocks are indivisible - only ', nb, &
        ' of ', mpi_size_, ' ranks have work, the rest idle. Reduce ranks or split the mesh.'
    else if (real(maxval(blk_ncells), R8) > ideal) then
      write(*,'(A,I0,A)') '  NOTE: balance is capped by the largest block (', &
        maxval(blk_ncells), ' cells); only a finer mesh split can improve it.'
    end if
  end subroutine report_partition_balance


  !> MPI_ALLREDUCE with MPI_SUM for a scalar real(R8).
  subroutine mpi_allreduce_sum_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_sum_r8


  !> MPI_ALLREDUCE with MPI_MIN for a scalar real(R8).
  subroutine mpi_allreduce_min_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_min_r8


  !> MPI_ALLREDUCE with MPI_MAX for a scalar real(R8).
  subroutine mpi_allreduce_max_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_max_r8


  !> MPI_GATHER of one real(R8) per rank onto root. `arr` must be at least
  !> mpi_size_ long on root; off root its contents are undefined.
  subroutine mpi_gather_r8(local_val, arr)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: arr(:)
#ifdef USE_MPI
    integer :: ierr
    call MPI_GATHER(local_val, 1, MPI_DOUBLE_PRECISION, &
                    arr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#else
    arr(1) = local_val
#endif
  end subroutine mpi_gather_r8


  !> MPI_ALLREDUCE with MPI_SUM for an array of real(R8): arr is overwritten with the global result.
  !> The send buffer is a copy of arr, not MPI_IN_PLACE.
  !> With gfortran the MPI sentinels in COMMON blocks (MPI_IN_PLACE, MPI_BOTTOM, MPI_STATUS_IGNORE,
  !> MPI_STATUSES_IGNORE, MPI_ERRCODES_IGNORE, MPI_ARGV_NULL, MPI_ARGVS_NULL) are not recognised while
  !> -fno-underscoring is set (cmake/SetFortranFlags.cmake): do not use them.
  subroutine mpi_allreduce_sum_r8_array(arr, n)
    integer, intent(in)     :: n
    real(R8), intent(inout) :: arr(n)
#ifdef USE_MPI
    integer :: ierr
    real(R8), allocatable :: send(:)
    send = arr
    call MPI_ALLREDUCE(send, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_allreduce_sum_r8_array


  !> Distributed L2 norm: local sum of squares + MPI_ALLREDUCE SUM + sqrt.
  function mpi_allreduce_norm2(x, n) result(global_norm)
    integer, intent(in) :: n
    real(R8), intent(in) :: x(n)
    real(R8) :: global_norm
    real(R8) :: local_ss

    local_ss = dot_product(x, x)
#ifdef USE_MPI
    call mpi_allreduce_sum_r8(local_ss, global_norm)
    global_norm = sqrt(global_norm)
#else
    global_norm = sqrt(local_ss)
#endif
  end function mpi_allreduce_norm2


  !> MPI_REDUCE with MPI_SUM for a scalar real(R8) to root rank 0: val holds the sum on the root and is
  !> left as it is elsewhere, where the receive buffer is not significant. Every rank sends a copy of val:
  !> no buffer passed twice (an OUT argument cannot be aliased with another one) and no MPI_IN_PLACE.
  !> With gfortran the MPI sentinels in COMMON blocks (MPI_IN_PLACE, MPI_BOTTOM, MPI_STATUS_IGNORE,
  !> MPI_STATUSES_IGNORE, MPI_ERRCODES_IGNORE, MPI_ARGV_NULL, MPI_ARGVS_NULL) are not recognised while
  !> -fno-underscoring is set (cmake/SetFortranFlags.cmake): do not use them.
  subroutine mpi_reduce_sum_r8(val)
    real(R8), intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    real(R8) :: send
    send = val
    call MPI_REDUCE(send, val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_reduce_sum_r8


  !> MPI_REDUCE with MPI_SUM for an array of real(R8) to root rank 0.
  !> Result is only valid on root. More efficient than ALLREDUCE when
  !> only root needs the result (e.g. for residual output). arr is left as it is elsewhere, where the
  !> receive buffer is not significant. Every rank sends a copy of arr: no buffer passed twice (an OUT
  !> argument cannot be aliased with another one) and no MPI_IN_PLACE.
  !> With gfortran the MPI sentinels in COMMON blocks (MPI_IN_PLACE, MPI_BOTTOM, MPI_STATUS_IGNORE,
  !> MPI_STATUSES_IGNORE, MPI_ERRCODES_IGNORE, MPI_ARGV_NULL, MPI_ARGVS_NULL) are not recognised while
  !> -fno-underscoring is set (cmake/SetFortranFlags.cmake): do not use them.
  subroutine mpi_reduce_sum_r8_array(arr, n)
    integer, intent(in)     :: n
    real(R8), intent(inout) :: arr(n)
#ifdef USE_MPI
    integer :: ierr
    real(R8), allocatable :: send(:)
    send = arr
    call MPI_REDUCE(send, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_reduce_sum_r8_array


  !> Broadcast a logical value from root (rank 0) to all ranks.
  subroutine mpi_bcast_logical(val)
    logical, intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    call MPI_BCAST(val, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_bcast_logical


  !> Broadcast an integer value from root (rank 0) to all ranks.
  subroutine mpi_bcast_integer(val)
    integer, intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    call MPI_BCAST(val, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_bcast_integer


  !> Broadcast a character string from root (rank 0) to all ranks: its length first, then its characters, so that every
  !> rank passes MPI_BCAST the count of root (MPI-4.1, section 6.4: the type signature must be the same on every
  !> process). The other ranks receive the string allocated to that length.
  subroutine mpi_bcast_chars(text)
    character(len=:), allocatable, intent(inout) :: text
    integer :: n
#ifdef USE_MPI
    integer :: ierr
#endif
    n = 0
    if (mpi_is_root .and. allocated(text)) n = len(text)
    call mpi_bcast_integer(n)
    if (.not. mpi_is_root) then
      if (allocated(text)) deallocate(text)
      allocate(character(len=n) :: text)
    end if
#ifdef USE_MPI
    if (n > 0) then
      call MPI_BCAST(text, n, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
    end if
#endif
  end subroutine mpi_bcast_chars


  !> Check MPI return code. Abort all ranks on error.
  subroutine check_mpi_error(ierr)
    integer, intent(in) :: ierr
#ifdef USE_MPI
    integer :: abort_ierr
    if (ierr /= 0) then
      write(*,'(A,I4,A,I6)') ' [RANK ', mpi_rank_, '] MPI error code: ', ierr
      call flush_output_units()
      call MPI_ABORT(MPI_COMM_WORLD, 1, abort_ierr)
    end if
#endif
  end subroutine check_mpi_error


  !> MPI-aware abort: prints message, writes out the output of this process (flush_output_units), then aborts all
  !> ranks. For an error that one rank finds alone, outside OpenMP regions, where no collective call can be reached
  !> (mpi_stop_all and mpi_stop_if_any are the stops of the errors that every rank, or the ranks together, can decide).
  subroutine mpi_abort_all(message)
    character(len=*), intent(in) :: message
#ifdef USE_MPI
    integer :: ierr
    write(*,'(A,I4,A,A)') ' [RANK ', mpi_rank_, '] ABORT: ', trim(message)
    call flush_output_units()
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#else
    write(*,'(A,A)') ' ABORT: ', trim(message)
    call flush_output_units()
    error stop 1
#endif
  end subroutine mpi_abort_all


  !> Write out what this process still holds for output_unit and error_unit, before MPI_ABORT. MPI_ABORT makes only a
  !> "best attempt" to abort the processes (MPI-4.1, section 11.4.2) and need not end this one through the end of the
  !> Fortran program, where the run-time library writes out its units: a run without this call can lose the message,
  !> what it wrote before and the lines on error_unit. It writes out the units of this process only: what the other
  !> ranks still hold when the run is aborted is not written out by this call.
  subroutine flush_output_units()
    flush(output_unit)
    flush(error_unit)
  end subroutine flush_output_units


  !> Stop of every rank on an error that every rank has found, a decision taken on all ranks (the outcome of reading
  !> input.ini, which rank 0 broadcasts; a file that every rank reads, such as bc.txt): each rank prints the message and
  !> flushes its output, the ranks wait for each other, so that what any rank printed before (rank 0 the explanation of
  !> the error) is out before any rank ends, then the ranks finalize MPI and end with exit code 1. MPI_ABORT
  !> (mpi_abort_all) makes only a "best attempt" to abort the processes (MPI-4.1, section 11.4.2): with Intel MPI
  !> 2021.12 the launcher at times ends before writing out what the ranks printed (measured), and the first rank to
  !> call it can stop another one before that one prints. A program that ends normally calls MPI_FINALIZE in every
  !> process, which can then give an exit code (MPI-4.1, section 11.2.2).
  !> Precondition: every rank must reach this call on the same error, from the same data (files that every rank reads at
  !> set-up, the outcome of a broadcast), outside the parallel regions of OpenMP. A decision that can differ between ranks,
  !> such as one on the initial condition, which each rank reads itself (Check_IC in Wrap_Setup), must be reduced first
  !> (mpi_stop_if_any); else the ranks that do not stop wait in their next collective call and the run hangs.
  subroutine mpi_stop_all(message)
    character(len=*), intent(in) :: message
#ifdef USE_MPI
    integer :: ierr
    write(*,'(A,I4,A,A)') ' [RANK ', mpi_rank_, '] STOP: ', trim(message)
    flush(output_unit)
    ! MPI_BARRIER returns at a rank only after every rank has entered it (MPI-4.1, section 6.3)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
#else
    write(*,'(A,A)') ' STOP: ', trim(message)
    flush(output_unit)
#endif
    error stop 1
  end subroutine mpi_stop_all

  !> Stop of the ranks on an error that only some of them can find (the initial condition or the restart solution that
  !> each rank reads, a file that a rank opens for its blocks): every rank calls this routine, found true on the ranks
  !> that found the error. The largest of the flag decides; the lowest rank that found the error writes detail (with its
  !> rank when the run has more than one), then every rank stops together with mpi_stop_all(message). Without an error it
  !> returns after one MPI_ALLREDUCE. Not for OpenMP regions.
  subroutine mpi_stop_if_any(found, message, detail)
    logical,          intent(in) :: found     !< The verdict of this rank
    character(len=*), intent(in) :: message   !< The line of the stop, as mpi_stop_all gives it
    character(len=*), intent(in) :: detail    !< What this rank found, written by the first rank that found it
    ! Local
    real(R8) :: flag, flag_all, first

    flag = 0.0_R8
    if (found) flag = 1.0_R8
    call mpi_allreduce_max_r8(flag, flag_all)
    if (flag_all <= 0.0_R8) return
    call mpi_allreduce_min_r8(merge(real(mpi_rank_, R8), huge(1.0_R8), found), first)
    if (found .and. nint(first) == mpi_rank_) then
      if (mpi_size_ > 1) then
        write(*,'(A,I0,A)') trim(detail)//' (rank ', mpi_rank_, ')'
      else
        write(*,'(A)') trim(detail)
      end if
    end if
    call mpi_stop_all(message)

  end subroutine mpi_stop_if_any


  !> Stop of the process from inside an OpenMP parallel region (the time step state check), where no MPI call is allowed
  !> (MPI is FUNNELED): writes the message, writes out the output (output_unit), and ends this process with exit code 1;
  !> the launcher ends the other ranks.
  subroutine stop_from_thread(message)
    character(len=*), intent(in) :: message

    write(*,'(A,I4,A,A)') ' [RANK ', mpi_rank_, '] STOP: ', trim(message)
    flush(output_unit)
    stop 1

  end subroutine stop_from_thread


end module FUSS_Mod_MPI
