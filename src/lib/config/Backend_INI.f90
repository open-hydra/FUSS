module FUSS_Backend_INI
  implicit none
  private
  public :: Load_Ini, Scan_Ini, Check_Unknown_Options

contains

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