module FUSS_Backend_INI
  use iso_fortran_env, only: R8 => real64
  implicit none
  private
  public :: Load_Ini, Scan_Ini, Check_Unknown_Options, Reload_Runtime_Ini

  !> Parameters that MAY be changed while a run is in progress.
  !>
  !> Everything not on this list is immutable, which is the safe default: a new
  !> parameter added later is immutable until someone deliberately decides it is
  !> not. The test is not "does it look harmless" but "is it READ LIVE" -- these
  !> seven are fetched from obj_sim_param / obj_io every step or every write, so
  !> assigning them mid-run actually changes behaviour.
  !>
  !> `iter-threshold` deliberately is NOT here, and the reason is the whole
  !> point of this list. It is copied once into domain%itermax at
  !> Wrap_Setup.f90:48 (and per level in Mod_Multigrid.f90:53), and the stopping
  !> test reads domain%itermax, not obj_sim_param%iter_threshold. Re-assigning
  !> it mid-run would therefore do exactly nothing, with no message -- the
  !> failure mode this whole check exists to remove. Treated as immutable, an
  !> attempt to change it is at least reported. Making it genuinely mutable
  !> means re-deriving itermax, which needs the domain and so does not belong in
  !> the config layer.
  !>
  !> This was not reasoned out in advance: the list first included
  !> `iter-threshold`, and test/run_runtime_ini_gate.sh hung because lowering it
  !> mid-run never stopped the run.
  character(len=*), parameter :: RUNTIME_MUTABLE(7) = [ character(len=16) :: &
    'res-threshold ', 'time-threshold',                                      &
    'sol-diter     ', 'sol-dtime     ', 'res-diter     ', 'shell-diter   ',  &
    'ini-diter     ' ]

  !> What input.ini said for each registered parameter the last time it was
  !> read. Compared against on every re-read, so that a mid-run edit is detected
  !> once, when it happens.
  character(len=256), allocatable, save :: ini_seen(:)

contains

  !>@brief Re-read input.ini mid-run, applying ONLY the parameters that may
  !>       legitimately change, and reporting any attempt to change the rest.
  !>
  !> WHAT WAS WRONG
  !> --------------
  !> Read_Inifile_Runtime used to call Load_Ini, which re-assigns EVERY
  !> registered parameter from the file, and then computed Validate_Registry()
  !> into a variable it never looked at. Three consequences:
  !>
  !>   * an out-of-range value introduced mid-run was accepted in silence;
  !>   * parameters that cannot change meaning mid-run -- `newrun`,
  !>     `integration-variables`, `law` -- were re-assignable, while the state
  !>     derived from them at setup (obj_mesh_motion%enabled, and so on) was not
  !>     re-derived;
  !>   * Check_Mesh_Motion_Compatibility runs only at setup, so editing
  !>     `irs = true` into a running moving-mesh case reached exactly the
  !>     configuration that is refused at startup.
  !>
  !> THE POLICY, AS DECIDED
  !> ----------------------
  !> A change to an immutable parameter is IGNORED, not fatal -- the run
  !> continues on the value it started with -- and is recorded. The caller
  !> writes the record to a file that is created only if something was actually
  !> ignored, so the file's existence is itself the signal.
  !>
  !> Note what this buys beyond the warning: because `irs`, `law` and
  !> `integration-variables` are never re-assigned, the guarded combinations
  !> cannot be reached mid-run at all, so the setup-time guards do not need to
  !> re-run. Do not "simplify" this into a plain re-read later without also
  !> re-running Check_Mesh_Motion_Compatibility.
  subroutine Reload_Runtime_Ini ( fini, iter, report )
    use Finer,               only: file_ini
    use FUSS_Input_Registry, only: reg
    implicit none
    type(file_ini),   intent(in)  :: fini
    integer,          intent(in)  :: iter
    character(len=*), intent(out) :: report
    ! Local
    integer             :: i, error
    character(len=256)  :: now
    character(len=256)  :: entry_

    report = ''
    if ( .not. allocated(ini_seen) ) return     ! Load_Ini has not run yet

    do i = 1, reg%size

      now = Raw_Ini_Text(fini, i, error)
      if ( now == ini_seen(i) ) cycle           ! the file did not change here

      if ( Is_Runtime_Mutable(reg%params(i)%name) ) then
        error = 1
        if (associated(reg%params(i)%value%i)) then
          call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%i, error=error)
        else if (associated(reg%params(i)%value%r)) then
          call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%r, error=error)
        else if (associated(reg%params(i)%value%l)) then
          call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%l, error=error)
        else if (associated(reg%params(i)%value%s)) then
          call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%s, error=error)
        endif
      else
        write(entry_,'(A,I0,A)') 'iteration ', iter, ':  ['//                  &
          trim(reg%params(i)%section)//'] '//trim(reg%params(i)%name)//        &
          '  ignored "'//trim(adjustl(now))//'", kept "'//                     &
          trim(adjustl(ini_seen(i)))//'"'
        if ( len_trim(report) + len_trim(entry_) + 1 < len(report) ) then
          report = trim(report)//trim(entry_)//new_line('a')
        endif
      endif

      ! Record what the file says NOW either way, so a change is reported once
      ! when it is made rather than on every re-read for the rest of the run.
      ini_seen(i) = now
    enddo

  end subroutine Reload_Runtime_Ini


  !> The option's text exactly as it stands in the file, or '<absent>'.
  !>
  !> Everything is compared as raw text, and the comparison is FILE-TO-FILE --
  !> against what the file said last time, not against the live variable. That
  !> matters: several registered variables are overwritten after setup with
  !> values that never came from the file at all. `integration-variables` is
  !> read as 'cons' and then rewritten to 'Conservative' by
  !> Assign_Integration_Variables once the procedure pointer is set. Comparing
  !> the file against the live variable therefore reported a spurious change on
  !> every single re-read -- a 4.2 MB report file from a run nobody had edited.
  !>
  !> The same fact is why this routine does NOT re-run Validate_Registry. After
  !> setup the registry's targets no longer hold raw ini values, so
  !> 'Conservative' fails its own 'cons , prim' allowed-list. The validation
  !> result was historically discarded, which is exactly what hid this.
  function Raw_Ini_Text ( fini, i, error ) result ( txt )
    use Finer,               only: file_ini
    use FUSS_Input_Registry, only: reg
    implicit none
    type(file_ini), intent(in)  :: fini
    integer,        intent(in)  :: i
    integer,        intent(out) :: error
    character(len=256)          :: txt
    character(len=256)          :: buf

    buf = ''
    call fini%get(reg%params(i)%section, reg%params(i)%name, val=buf, error=error)
    if ( error /= 0 ) then
      txt = '<absent>'
    else
      txt = adjustl(buf)
    endif

  end function Raw_Ini_Text


  !> Fill the file-state snapshot. Called at the end of Load_Ini, i.e. once the
  !> setup read has happened and before anything has had a chance to overwrite
  !> a registered variable with a derived value.
  subroutine Snapshot_Ini_Text ( fini )
    use Finer,               only: file_ini
    use FUSS_Input_Registry, only: reg
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: i, error

    if ( allocated(ini_seen) ) deallocate( ini_seen )
    allocate( ini_seen(reg%size) )
    do i = 1, reg%size
      ini_seen(i) = Raw_Ini_Text(fini, i, error)
    enddo

  end subroutine Snapshot_Ini_Text


  logical function Is_Runtime_Mutable ( name ) result ( ok )
    implicit none
    character(len=*), intent(in) :: name
    integer :: k
    ok = .false.
    do k = 1, size(RUNTIME_MUTABLE)
      if ( trim(name) == trim(RUNTIME_MUTABLE(k)) ) then
        ok = .true. ; return
      endif
    enddo
  end function Is_Runtime_Mutable


  !>@brief Report input.ini content that FUSS will silently not apply.
  !>
  !> Load_Ini pulls values BY NAME for every registered parameter, so a key
  !> nobody registered is simply never looked at, and Validate_Registry cannot
  !> see it either -- it iterates the registry, not the file. The consequence is
  !> that `vnn` misspelt as `vnnn` does not fail: it silently reverts to the
  !> default and the run produces a plausible answer at the wrong time step.
  !>
  !> Not hypothetical. Two of the seven guard demonstrations in
  !> test/run_guard_demos.sh did not fire the first time they were written,
  !> because their parameters were being written into a section the registry
  !> never reads. Both runs completed happily.
  !>
  !> SCOPE. Only sections belonging to this solver are checked -- those prefixed
  !> `<codename>-` or `HYDRA-`. A FUSS input.ini also carries [GRIB-*], [GPB-*],
  !> [ICB-*], [BCB-*] and arbitrarily named per-BC blocks such as [Twall-hot],
  !> which belong to ATLAS and which FUSS never reads. Checking those would
  !> reject every existing case.
  subroutine Check_Unknown_Options ( fini, out )
    use Finer,               only: file_ini
    use FUSS_Input_Registry, only: reg
    use FUSS_Global_m,       only: codename
    implicit none
    type(file_ini),   intent(in)  :: fini
    character(len=*), intent(out) :: out
    ! Local
    character(len=:), allocatable :: sections(:), pairs(:)
    character(len=256) :: sec, opt
    integer :: s, i, nbad
    logical :: known

    out   = ''
    nbad  = 0

    call fini%get_sections_list ( sections )
    if ( .not. allocated(sections) ) return

    ! A REPEATED section header is silently ignored -- everything in the second
    ! occurrence is dropped. Measured: adding a second [FUSS-Numerics] carrying
    ! `vnn = 0.123` left dt at the value implied by `vnn = 0.5` in the first.
    ! The option-loop below cannot see it either, because it enumerates the
    ! first occurrence only, so a duplicate would otherwise slip past BOTH
    ! checks.
    do s = 1, size(sections)
      sec = adjustl(sections(s))
      if ( .not. Section_Is_Ours(sec) ) cycle
      do i = 1, s - 1
        if ( Same(sections(i), sec) ) then
          nbad = nbad + 1
          if ( nbad == 1 ) then
            out = '[ERROR] repeated section in input.ini: ['//trim(sec)//']'
          else
            out = trim(out)//' ; repeated section ['//trim(sec)//']'
          endif
          exit
        endif
      enddo
    enddo

    do s = 1, size(sections)
      sec = adjustl(sections(s))
      if ( .not. Section_Is_Ours(sec) ) cycle

      ! Each section's loop must be run to completion: FiNeR keeps the cursor in
      ! saved state, so abandoning one part-way would corrupt the next.
      do while ( fini%loop( section_name=trim(sec), option_pairs=pairs ) )
        opt = adjustl(pairs(1))

        known = .false.
        do i = 1, reg%size
          if ( Same(reg%params(i)%section, sec) .and. Same(reg%params(i)%name, opt) ) then
            known = .true.
            exit
          endif
        enddo

        ! `levels` is read directly by Scan_Ini, before the registry exists --
        ! the registry cannot contain it, because how many multigrid levels to
        ! register is exactly what it answers.
        if ( .not. known ) then
          if ( Same(opt,'levels') .and. index(sec,'-Multigrid') > 0 ) known = .true.
        endif

        if ( .not. known ) then
          nbad = nbad + 1
          if ( nbad == 1 ) then
            out = '[ERROR] unknown option in input.ini: ['//trim(sec)//'] '//trim(opt)
          else
            out = trim(out)//' ; ['//trim(sec)//'] '//trim(opt)
          endif
        endif
      enddo
    enddo

  contains

    !> Does this section belong to this solver? Anything else in the file is
    !> another tool's business.
    logical function Section_Is_Ours ( name ) result ( ours )
      character(len=*), intent(in) :: name
      ours = ( index(name, trim(codename)//'-') == 1 ) .or. ( index(name, 'HYDRA-') == 1 )
    end function Section_Is_Ours

    !> Trimmed, case-insensitive comparison. Case-insensitive on purpose: a
    !> rejection over capitalisation would be a worse failure than the one this
    !> is meant to catch.
    logical function Same ( a, b ) result ( eq )
      character(len=*), intent(in) :: a, b
      character(len=256) :: x, y
      integer :: k, c
      x = adjustl(a); y = adjustl(b)
      do k = 1, len_trim(x)
        c = iachar(x(k:k)); if (c >= 65 .and. c <= 90) x(k:k) = achar(c+32)
      enddo
      do k = 1, len_trim(y)
        c = iachar(y(k:k)); if (c >= 65 .and. c <= 90) y(k:k) = achar(c+32)
      enddo
      eq = ( trim(x) == trim(y) )
    end function Same

  end subroutine Check_Unknown_Options


  subroutine Load_Ini(fini)
    use FUSS_Input_Registry, only: reg
    use Finer, only: file_ini
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: i, error

    do i = 1, reg%size
      error = 1

      ! Each registry entry has exactly one associated typed pointer.
      if (associated(reg%params(i)%value%i)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%i, error=error)
      else if (associated(reg%params(i)%value%r)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%r, error=error)
      else if (associated(reg%params(i)%value%l)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%l, error=error)
      else if (associated(reg%params(i)%value%s)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%s, error=error)
      else if (associated(reg%params(i)%value%iarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%iarr, error=error)
      else if (associated(reg%params(i)%value%rarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%rarr, error=error)
      end if

      if (error == 0) reg%params(i)%is_set = .true.
    end do

    ! Record the file's own text for every parameter, for Reload_Runtime_Ini to
    ! compare against later. Must happen HERE, before anything downstream
    ! overwrites a registered variable with a derived value.
    call Snapshot_Ini_Text ( fini )

  end subroutine


  subroutine Scan_Ini(fini, nprobes, probes_name, nmgl)
    use Finer, only: file_ini
    use FUSS_Config_Types_m, only: obj_sim_param
    use FUSS_Global_m
    use IR_precision
    implicit none
    type(file_ini), intent(in) :: fini
    integer, intent(out) :: nprobes, nmgl
    character(len=16), allocatable, intent(out) :: probes_name(:)
    ! Local
    character(len=llen) :: wholestring
    integer :: error, i
  
    ! Read the MGL
    nmgl = 1
    if (.not. obj_sim_param%HYDRA_MG) then
      call fini%get(section_name=trim(codename)//'-Multigrid', option_name='levels', val=nmgl, error=error)
    else
      call fini%get(section_name='HYDRA-Multigrid', option_name='levels', val=nmgl, error=error)
    endif

    ! Count the probes
    nprobes = 0
    do 
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,nprobes+1)), val=wholestring, error=error)
      if (error/=0) exit
      nprobes = nprobes+1
    enddo

    allocate(probes_name(nprobes))
    do i = 1, nprobes
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,i)), val=probes_name(i), error=error)
    enddo

  end subroutine Scan_Ini

end module FUSS_Backend_INI