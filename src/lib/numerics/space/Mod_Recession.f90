!> The ablating wall's per-step and per-output bookkeeping (plan 08 sections
!> 2.1, 2.3 R3, 3.2): the env values of every 306 face at the step's time, the
!> restart state of the faces, and the two files the faces write.
!>
!> ENV VALUES. The env table of a 306 face (t p_w rhoeUeCH0 h_r q_rad chem
!> hconv T_edge) is evaluated ONCE per step, after the clock has advanced to
!> t^(n+1) and before the type-2 pressure solve, which reads the face's p_w as
!> its Dirichlet value: a refresh inside the BC routine (stage time) would hand
!> the backward-Euler p^(n+1) the previous step's p_w, and p-wall on step 1.
!> Every RK stage runs at t^(n+1), so the stored vector is what a stage-time
!> evaluation would give. Also called once at setup, after the clock is set.
!>
!> RESTART STATE. Besides the node positions the solution file carries, a
!> face's balance has state: the recession so far, the wall temperature (the
!> Newton start), the recession rate (it moves the mesh one step later), the
!> blowing correction (the Picard start) and the B'_pg of the last chemistry-
!> on solve (the fractional-chem rule). OUTPUT/<prefix>bc306-state.dat holds
!> them, rewritten at every solution output with 17 significant digits (an
!> exact binary64 round trip: the converged T_w depends on its start at the
!> seb-tol level) and read back when newrun = false.
!>
!> SEB FILE. OUTPUT/<prefix>seb.txt: one row per 306 face per solution output
!> with every term of the last converged balance (plan 08 section 3.2).
module FUSS_Mod_Recession
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Ablation_Update_Env, Setup_Ablation, Write_Ablation_Outputs

  logical, save :: seb_file_started = .false.

contains

  !> Env values of every 306 face at domain%time (serial; outside the
  !> parallel region).
  subroutine Ablation_Update_Env ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_ablation
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    ! Local
    integer  :: i
    real(R8) :: v(7)

    if ( .not. obj_ablation%enabled ) return
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 306 ) cycle
      v = domain%bc(i)%env%eval( domain%time )
      domain%bc(i)%pw        = v(1)
      domain%bc(i)%rhoeUeCH0 = v(2)
      domain%bc(i)%hr        = v(3)
      domain%bc(i)%qrad_env  = v(4)
      domain%bc(i)%chem      = v(5)
      domain%bc(i)%hconv_env = v(6)
      domain%bc(i)%Tedge     = v(7)
    enddo

  end subroutine Ablation_Update_Env


  !> Once at setup, after Setup_BC: on a restart, the faces' state from
  !> OUTPUT/<prefix>bc306-state.dat (every 306 face must have its row, and
  !> every row its face).
  subroutine Setup_Ablation ( domain, newrun )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_ablation
    use FUSS_Global_m,       only: FUSS_phase_prefix
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    logical,                intent(in)    :: newrun
    ! Local
    character(len=512) :: file, line
    integer  :: u, ios, ln, i, b, ii, jj, kk, f, nrows
    real(R8) :: v(5)
    logical, allocatable :: seen(:)

    if ( .not. obj_ablation%enabled ) return
    if ( newrun ) return

    file = 'OUTPUT/'//trim(FUSS_phase_prefix)//'bc306-state.dat'
    open( newunit=u, file=trim(file), status='old', action='read', iostat=ios )
    if ( ios /= 0 ) then
      write(*,'(A)') ' [ERROR] restart with ablating faces (306): '//trim(file)//' is missing'
      write(*,'(A)') '         it holds each face''s recession, wall temperature, recession rate and Picard start'
      error stop 1
    endif
    allocate( seen(domain%nbound) )
    seen = .false.
    ln = 0 ; nrows = 0
    do
      read( u, '(A)', iostat=ios ) line
      if ( ios /= 0 ) exit
      ln = ln + 1
      line = adjustl(line)
      if ( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
      read( line, *, iostat=ios ) b, ii, jj, kk, f, v
      if ( ios /= 0 ) call Fail( 'expected "b i j k f s_total Tw_prev sdot beta Bpg_lag"' )
      nrows = nrows + 1
      do i = 1, domain%nbound
        if ( domain%bc(i)%type /= 306 ) cycle
        if ( domain%bc(i)%b == b .and. domain%bc(i)%i == ii .and. domain%bc(i)%j == jj .and. &
             domain%bc(i)%k == kk .and. domain%bc(i)%f == f ) exit
      enddo
      if ( i > domain%nbound ) call Fail( 'no 306 face of bc.txt matches this row' )
      if ( seen(i) ) call Fail( 'a second row for the same face' )
      seen(i) = .true.
      domain%bc(i)%s_total = v(1)
      domain%bc(i)%Tw_prev = v(2)
      domain%bc(i)%sdot    = v(3)
      domain%bc(i)%beta    = v(4)
      domain%bc(i)%Bpg_lag = v(5)
      domain%bc(i)%Tw_face = v(2)
    enddo
    close( u )
    do i = 1, domain%nbound
      if ( domain%bc(i)%type == 306 .and. .not. seen(i) ) then
        write(*,'(A,I0,A)') ' [ERROR] '//trim(file)//': no row for the 306 face of bc.txt record ', i, &
          ' (written by a run with another bc.txt?)'
        error stop 1
      endif
    enddo

  contains

    subroutine Fail ( reason )
      character(len=*), intent(in) :: reason
      write(*,'(A,I0,A)') ' [ERROR] '//trim(file)//' line ', ln, ': '//reason
      error stop 1
    end subroutine Fail

  end subroutine Setup_Ablation


  !> At every solution output (root): the seb.txt rows and the restart state.
  subroutine Write_Ablation_Outputs ( domain )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_ablation, obj_sim_param
    use FUSS_Global_m,       only: FUSS_phase_prefix
    implicit none
    type(FUSS_domain_type), intent(in) :: domain
    ! Local
    character(len=512) :: file
    integer :: u, ios, i

    if ( .not. obj_ablation%enabled ) return

    ! seb.txt: appended, header on a new run's first write
    file = 'OUTPUT/'//trim(FUSS_phase_prefix)//'seb.txt'
    if ( .not. seb_file_started ) then
      if ( obj_sim_param%newrun ) then
        open( newunit=u, file=trim(file), status='replace', action='write', iostat=ios )
        if ( ios == 0 ) write(u,'(A)') '# time b i j k f T_w T_c beta Bpg Bch h_w rhoeUeCH0 h_r mdot_pg mdot_ch iter '// &
          'nclamp_p nclamp_Bg nclamp_T nclamp_rc  cond conv adv_pg adv_ch alpha_qrad eps_Tinf4 m_eps_Tw4 F   (SI)'
      else
        open( newunit=u, file=trim(file), status='unknown', position='append', action='write', iostat=ios )
      endif
      seb_file_started = .true.
    else
      open( newunit=u, file=trim(file), status='old', position='append', action='write', iostat=ios )
    endif
    if ( ios /= 0 ) then
      write(*,'(A)') ' [ERROR] ablating wall: cannot open '//trim(file)
      error stop 1
    endif
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 306 ) cycle
      associate ( bc => domain%bc(i) )
        write(u,'(ES24.16,5I8,10ES24.16,5I8,8ES24.16)') domain%time, bc%b, bc%i, bc%j, bc%k, bc%f, &
          bc%Tw_prev, domain%blk(bc%b)%T(bc%i,bc%j,bc%k), bc%beta, bc%Bpg, bc%Bch, bc%hw, bc%rhoeUeCH0, bc%hr, &
          bc%mdot_pg, bc%mdot_ch, bc%seb_iter, bc%nclamp, bc%seb
      end associate
    enddo
    close(u)

    ! bc306-state.dat: the restart state, rewritten at every output
    file = 'OUTPUT/'//trim(FUSS_phase_prefix)//'bc306-state.dat'
    open( newunit=u, file=trim(file), status='replace', action='write', iostat=ios )
    if ( ios /= 0 ) then
      write(*,'(A)') ' [ERROR] ablating wall: cannot open '//trim(file)
      error stop 1
    endif
    write(u,'(A)') '# b i j k f s_total Tw_prev sdot beta Bpg_lag'
    write(u,'(A,ES24.16)') '# min-thickness ', 0.0d0
    do i = 1, domain%nbound
      if ( domain%bc(i)%type /= 306 ) cycle
      associate ( bc => domain%bc(i) )
        write(u,'(5I8,5ES24.16)') bc%b, bc%i, bc%j, bc%k, bc%f, bc%s_total, bc%Tw_prev, bc%sdot, bc%beta, bc%Bpg_lag
      end associate
    enddo
    close(u)

  end subroutine Write_Ablation_Outputs

end module FUSS_Mod_Recession
