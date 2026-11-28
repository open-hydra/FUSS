module FUSS_Backend_INI
  use iso_fortran_env, only: R8 => real64
  implicit none
  private
  public :: Load_Ini, Load_Param, Scan_Ini, Snapshot_Ini_Text, Check_Unknown_Options, Check_Value_Lengths, Reload_Runtime_Ini

  !> Parameters that MAY be changed while a run is in progress.
  !>
  !> Everything not on this list is immutable, which is the safe default: a new
  !> parameter added later is immutable until someone deliberately decides it is
  !> not. The test is not "does it look harmless" but "is it READ LIVE" -- these
  !> eight are fetched from obj_sim_param / obj_io every step or every write, so
  !> assigning them mid-run actually changes behaviour.
  !>
  !> iter-threshold joined the list with plan 14 S6: FUSS_postprocess re-derives
  !> domain%itermax from it after every reading (upstream f3401c5). With more than
  !> one multigrid level it has no effect (each level stops at its level<m>-iter,
  !> which stays immutable).
  character(len=*), parameter :: RUNTIME_MUTABLE(8) = [ character(len=16) :: &
    'res-threshold ', 'time-threshold', 'iter-threshold',                    &
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
  !>
  !> Plan 14 S6 (D8-A2): called only for a reading that passed the whole
  !> validation (Read_Inifile_Runtime); the reading itself comes from rank 0 on
  !> every rank (Load_Input_File).
  subroutine Reload_Runtime_Ini ( fini, iter, report )
    use Finer,               only: file_ini
    use FUSS_Input_Registry, only: reg, Validate_Param
    implicit none
    type(file_ini),   intent(inout)  :: fini
    integer,          intent(in)  :: iter
    character(len=*), intent(out) :: report
    ! Local
    integer             :: i, error
    character(len=256)  :: now
    character(len=256)  :: entry_
    character(len=1024) :: vout
    integer             :: old_i
    real(R8)            :: old_r
    logical             :: old_l
    character(len=1024) :: old_s

    report = ''
    if ( .not. allocated(ini_seen) ) return     ! Load_Ini has not run yet

    do i = 1, reg%size

      now = Raw_Ini_Text(fini, i, error)
      if ( now == ini_seen(i) ) cycle           ! the file did not change here

      if ( Is_Runtime_Mutable(reg%params(i)%name) ) then
        ! Apply, then validate the ONE parameter that changed against its own
        ! rule -- the same rule it passed at setup. If it fails, put the old
        ! value back and report; a bad edit never takes effect. This is
        ! per-parameter on purpose: the whole registry cannot be re-validated
        ! after setup, see Raw_Ini_Text below.
        if      (associated(reg%params(i)%value%i)) then ; old_i = reg%params(i)%value%i
        else if (associated(reg%params(i)%value%r)) then ; old_r = reg%params(i)%value%r
        else if (associated(reg%params(i)%value%l)) then ; old_l = reg%params(i)%value%l
        else if (associated(reg%params(i)%value%s)) then ; old_s = reg%params(i)%value%s
        endif

        call Load_Param ( fini, i )

        vout = ''
        call Validate_Param(i, vout)
        if ( len_trim(vout) > 0 ) then
          if      (associated(reg%params(i)%value%i)) then ; reg%params(i)%value%i = old_i
          else if (associated(reg%params(i)%value%r)) then ; reg%params(i)%value%r = old_r
          else if (associated(reg%params(i)%value%l)) then ; reg%params(i)%value%l = old_l
          else if (associated(reg%params(i)%value%s)) then ; reg%params(i)%value%s = old_s
          endif
          write(entry_,'(A,I0,A)') 'iteration ', iter, ':  '//trim(vout)//      &
            '  -- rejected "'//trim(adjustl(now))//'", kept "'//                 &
            trim(adjustl(ini_seen(i)))//'"'
          if ( len_trim(report) + len_trim(entry_) + 1 < len(report) ) then
            report = trim(report)//trim(entry_)//new_line('a')
          endif
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
  !> The comparison stays file-to-file (a registered variable may hold a derived
  !> value after set-up).
  function Raw_Ini_Text ( fini, i, error ) result ( txt )
    use Finer,               only: file_ini
    use FUSS_Input_Registry, only: reg
    implicit none
    type(file_ini), intent(inout)  :: fini
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


  !> Fill the file-state snapshot. Called by Read_Inifile right after the set-up
  !> Load_Ini (not inside Load_Ini, which every reading during the run calls too),
  !> before anything has had a chance to overwrite a registered variable with a
  !> derived value.
  subroutine Snapshot_Ini_Text ( fini )
    use Finer,               only: file_ini
    use FUSS_Input_Registry, only: reg
    implicit none
    type(file_ini), intent(inout) :: fini
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
  !> `<codename>-`, and [HYDRA-Multigrid], where FUSS reads its multigrid levels
  !> when hydra drives it (Register_Multigrid_Levels). A FUSS input.ini also carries
  !> [GRIB-*], [GPB-*], [ICB-*], [BCB-*] and arbitrarily named per-BC blocks such as
  !> [Twall-hot], which belong to ATLAS and which FUSS never reads. Checking those
  !> would reject every existing case. The other [HYDRA-*] sections ([HYDRA-Parameters],
  !> [HYDRA-IO]) are hydra's: until 2026-10-14 every HYDRA- section counted as FUSS's,
  !> so FUSS refused each key of a hydra case's input.ini (hydra-MF stopped at setup).
  subroutine Check_Unknown_Options ( fini, out )
    use Finer,               only: file_ini
    use FUSS_Input_Registry, only: reg
    use FUSS_Global_m,       only: codename
    implicit none
    type(file_ini),   intent(inout)  :: fini
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

      ! Each section's loop must be run to completion: FiNeR keeps the loop cursor
      ! in the section object (v2.1.0+, hence intent(inout) on fini; before, in
      ! saved module state), so abandoning one part-way would corrupt the next.
      do while ( fini%loop( section_name=trim(sec), option_pairs=pairs ) )
        opt = adjustl(pairs(1))
        ! "key =" without a value: FiNeR lists it but gives no value; Warn_Unread_Keys
        ! reports it with its line (upstream 189ba21): a warning, not a stop (plan 14 D-S6-3)
        if ( len_trim(pairs(2)) == 0 ) cycle

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
      ours = ( index(name, trim(codename)//'-') == 1 ) .or. Same(name, 'HYDRA-Multigrid')
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


  !> A string value longer than the variable it is read into is cut by FiNeR
  !> without a message. Measured 2026-10-15: a 21-character probe name became 16
  !> characters, the probe's own section was then looked up under the cut name,
  !> and the run stopped later on an unrelated "unknown probe variable". Every
  !> string key is checked; its limit is the length of its target.
  subroutine Check_Value_Lengths ( fini, out )
    use Finer,               only: file_ini
    use FUSS_Input_Registry, only: reg
    implicit none
    type(file_ini),   intent(inout) :: fini
    character(len=*), intent(out)   :: out
    ! Local
    character(len=4096) :: buf
    character(len=16)   :: sn, sl
    integer :: i, error, n, limit, nbad

    out  = ''
    nbad = 0
    do i = 1, reg%size
      if ( .not. associated(reg%params(i)%value%s) ) cycle
      buf = ''
      call fini%get(reg%params(i)%section, reg%params(i)%name, val=buf, error=error)
      if ( error /= 0 ) cycle
      n     = len_trim(adjustl(buf))
      limit = len(reg%params(i)%value%s)
      if ( n <= limit ) cycle
      nbad = nbad + 1
      write(sn,'(I0)') n
      write(sl,'(I0)') limit
      if ( nbad == 1 ) then
        out = '[ERROR] value too long in input.ini: ['//trim(reg%params(i)%section)//'] '// &
              trim(reg%params(i)%name)//' has '//trim(sn)//' characters, the limit is '//trim(sl)
      else
        out = trim(out)//' ; ['//trim(reg%params(i)%section)//'] '//trim(reg%params(i)%name)// &
              ' has '//trim(sn)//', limit '//trim(sl)
      endif
    enddo

  end subroutine Check_Value_Lengths


  subroutine Load_Ini(fini)
    use FUSS_Input_Registry, only: reg
    use Finer, only: file_ini
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: i

    do i = 1, reg%size
      call Load_Param(fini, i)
    end do

  end subroutine


  !> One key of the registry from input.ini (upstream 189ba21): its value, is_set, unread (a value that does not read
  !> as its type, an integer written with a period or an exponent) and the text as written
  subroutine Load_Param(fini, i)
    use FUSS_Input_Registry, only: reg
    use Finer, only: file_ini, err_option_vals
    implicit none
    type(file_ini), intent(in) :: fini
    integer,        intent(in) :: i
    integer :: error, error_text
    character(len=:), allocatable :: text

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

    ! A key of the file whose value the reader cannot take as the type of its parameter (a value that does not
    ! convert, a list with more values than the parameter holds, no value: FiNeR gives err_option_vals, and leaves
    ! the parameter as it was) is an error of the file, which Validate_Registry gives with the value as written
    reg%params(i)%unread = ( error == err_option_vals )
    if (allocated(reg%params(i)%text)) deallocate(reg%params(i)%text)
    call fini%get_string(section_name=reg%params(i)%section, option_name=reg%params(i)%name, val=text, error=error_text)
    if (error_text == 0) reg%params(i)%text = text
    ! An integer, or a list of integers, written with a period or an exponent is not one (the input form of an integer is
    ! a signed digit string): gfortran refuses it, ifx takes it and cuts the fraction (2.7 reads 2), so the text decides,
    ! in every build
    if (error == 0 .and. error_text == 0 .and. &
        (associated(reg%params(i)%value%i) .or. associated(reg%params(i)%value%iarr))) then
      if (scan(text, '.eEdDqQ') > 0) reg%params(i)%unread = .true.
    end if

  end subroutine Load_Param


  subroutine Scan_Ini(fini, nprobes, probes_name, nmgl)
    use Finer, only: file_ini, err_option_vals
    use FUSS_Input_Registry, only: reg
    use FUSS_Config_Types_m, only: obj_sim_param, obj_pyrolysis, obj_darcy
    use FUSS_Global_m
    use IR_precision
    implicit none
    type(file_ini), intent(inout) :: fini
    integer, intent(out) :: nprobes, nmgl
    character(len=16), allocatable, intent(out) :: probes_name(:)
    ! Local
    character(len=llen) :: wholestring
    integer :: error, error_text, i
    character(len=:), allocatable :: section, text
  
    ! Read the MGL
    nmgl = 1
    if (.not. obj_sim_param%HYDRA_MG) then
      section = trim(codename)//'-Multigrid'
    else
      section = 'HYDRA-Multigrid'
    endif
    call fini%get(section_name=section, option_name='levels', val=nmgl, error=error)
    ! A value that is not an integer, or no value: the run would go on with one grid level, without a word; an integer
    ! written with a period or an exponent is not one either (gfortran refuses it, ifx takes it and cuts the fraction).
    ! levels is no key of the registry: Validate_Registry gives this error first (Check_Input), as those of its keys
    text = ''
    call fini%get_string(section_name=section, option_name='levels', val=text, error=error_text)
    if (error == err_option_vals .or. (error == 0 .and. error_text == 0 .and. scan(text, '.eEdDqQ') > 0)) then
      if (len_trim(text) > 0) then
        reg%load_error = '[ERROR] input.ini, ['//section//'] levels = '//trim(text)//': the value is not an integer'
      else
        reg%load_error = '[ERROR] input.ini, ['//section//'] levels =: the value is empty'
      endif
    endif

    ! Pyrolysis (plan 06 section 3.5): the feature switch is the PRESENCE of
    ! its section. Decided here, before the registry exists, so that
    ! Register_Pyrolysis and every later setup step see it; the mid-run re-read
    ! never calls Scan_Ini, so adding or removing the section mid-run does nothing.
    obj_pyrolysis%enabled = fini%has_section( section_name=trim(codename)//'-Pyrolysis' )
    ! Plan 07: the Darcy section is only READ under gas-closure = type2; its
    ! presence is recorded so that Check_Pyrolysis_Compatibility can refuse a
    ! section nothing would apply.
    obj_darcy%section_present = fini%has_section( section_name=trim(codename)//'-Darcy' )

    ! Count the probes: probe1, probe2, ... with no gaps
    nprobes = 0
    do 
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,nprobes+1)), val=wholestring, error=error)
      if (error/=0) exit
      nprobes = nprobes+1
    enddo
    ! A probe numbered past a gap is an error of the file: Validate_Registry gives it (Check_Input)
    do i = nprobes+2, nprobes+100
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,i)), val=wholestring, error=error)
      if (error == 0) then
        reg%load_error = '[ERROR] input.ini, ['//trim(codename)//'-Probes] probe'//trim(str(.true.,i))//' found but probe'// &
                         trim(str(.true.,nprobes+1))//' is missing: the probes must be numbered without gaps'
        exit
      endif
    enddo

    allocate(probes_name(nprobes))
    do i = 1, nprobes
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,i)), val=probes_name(i), error=error)
    enddo

  end subroutine Scan_Ini

end module FUSS_Backend_INI