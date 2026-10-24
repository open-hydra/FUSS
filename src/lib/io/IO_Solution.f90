module FUSS_IO_Solution
  use FUSS_Mod_MPI, only: fuss_abort

  implicit none

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(r_solution_if), pointer, public :: Read_IC
  procedure(w_solution_if), pointer, public :: Write_solution

  !> Abstract interface relative to the finite-rate reactions source procedure
  abstract interface
    subroutine r_solution_if ( IOfield )
      use Lib_ORION_data
      implicit none
      type(ORION_data), intent(inout) :: IOfield
    end subroutine r_solution_if

    subroutine w_solution_if ( domain, IOfield, file )
      use FUSS_Advanced_Types_m
      use FUSS_Config_Types_m, only: obj_io
      use FUSS_Global_m
      use FUSS_Parameters_m
      use IR_Precision
      use Lib_ORION_data
      use Lib_VTK
      use Lib_Tecplot
      use FUSS_Lib_Solid
      implicit none
      type(FUSS_domain_type), intent(inout) :: domain
      type(ORION_data), intent(inout)    :: IOfield
      character(llen), intent(in)        :: file
    end subroutine w_solution_if
  end interface

contains

  subroutine Setup_Input_Solution()
    use FUSS_Config_Types_m, only: obj_sim_param, obj_io
    use FUSS_Global_m
    use FUSS_Parameters_m
    use IR_Precision
    implicit none
    logical         :: present
    integer         :: i
    character(llen) :: try
    character(6)    :: extension

    if (obj_sim_param%newrun) then
      if (index(obj_io%ic_format,'vtk')>0) then
        extension = '.vtm'
        Read_IC  => Read_vtk_tec
      else
        if (index(obj_io%ic_format,'ascii')>0) then
          extension = '.tec'
        else
          extension = '.szplt'
        endif
        Read_IC  => Read_vtk_tec
      endif
      obj_io%nameinit   = 'INPUT/'//trim(FUSS_phase_prefix)//'ic'//extension
      inquire(file=obj_io%nameinit, exist=present)

    else

      if (index(obj_io%sol_format,'vtk')>0) then
        extension = '.vtm'
        Read_IC => Read_vtk_tec
      else
        if (index(obj_io%sol_format,'ascii')>0) then
          extension = '.tec'
        else
          extension = '.szplt'
        endif
        Read_IC => Read_vtk_tec
      end if
      obj_io%nameinit = 'OUTPUT/'//trim(FUSS_phase_prefix)//'field'//extension
      inquire(file=obj_io%nameinit, exist=present)

      if (.not.present) then
        i = 0
        do
          i = i+1
          try = 'OUTPUT/'//trim(FUSS_phase_prefix)//'field'//trim(str(.true.,i))//extension
          inquire(file=try,exist=present)
          if (present) then
            obj_io%nameinit = try
          else
            exit
          endif
        enddo
      endif

    endif

    ! The volumetric source is an INPUT in both cases, so its name is set for
    ! both. It used to be assigned only on the newrun branch, which left it
    ! undefined on a restart: Read_vtk_tec then opened a garbage filename, the
    ! read failed, and qvol was silently set to zero for the whole restarted
    ! run. Nothing reported it -- a restart of a case with a volumetric source
    ! simply continued without the source.
    obj_io%namesource = 'INPUT/'//trim(FUSS_phase_prefix)//'st'//extension

  end subroutine Setup_Input_Solution

  subroutine Read_vtk_tec ( IOfield )
    use FUSS_Config_Types_m, only: obj_sim_param, obj_io, obj_pyrolysis, obj_darcy
    use FUSS_Parameters_m
    use FUSS_Mod_MPI, only: mpi_is_root
    use Lib_ORION_data
    use Lib_VTK
    use Lib_Tecplot
    use strings, only: parse
    implicit none
    type(ORION_data), intent(inout) :: IOfield
    ! Local
    type(ORION_data) :: IOinit, IOsource
    integer          :: error_init, error_source, error_dim, b, nblocks, nchi, nv, r, ndar
    character(clen)  :: format(2)
    character(len=32) :: expected

    error_init   = 0
    error_source = 0
    error_dim    = 0

    if ( obj_sim_param%newrun ) then
      call parse(obj_io%ic_format,' ', format)
    else
      call parse(obj_io%sol_format,' ', format)
    endif

    obj_io%IOtime = 0.d0
    
    select case(trim(format(1)))
    case('tecplot')
      IOfield%tec%format  = trim(format(2))
      IOinit%tec%format   = trim(format(2))
      IOsource%tec%format = trim(format(2))
      error_init   = tec_read_structured_multiblock(orion=IOinit,  filename=trim(obj_io%nameinit))
      error_source = tec_read_structured_multiblock(orion=IOsource,filename=trim(obj_io%namesource))
    case('vtk')
      IOfield%vtk%format  = trim(format(2))
      IOinit%vtk%format   = trim(format(2))
      IOsource%vtk%format = trim(format(2))
      error_init   = vtk_read_structured_multiblock(orion=IOinit,  vtmpath=obj_io%nameinit(1:len(trim(obj_io%nameinit))-4),vtspath='INPUT/vtk/field',time=obj_io%IOtime)
      error_source = vtk_read_structured_multiblock(orion=IOsource,vtmpath=obj_io%namesource(1:len(trim(obj_io%namesource))-4),vtspath='INPUT/vtk/field',time=obj_io%IOtime)
    end select

    if (error_init/=0) then
      obj_io%error_message = "[ERROR] reading input file "//trim(obj_io%nameinit)
      return
    endif

    ! Allocate IOfield (plus one slot per reaction when pyrolysis is on, plus
    ! the four Darcy state slots under type 2: p_pg, T_prev, eta_prev, eta_st)
    nchi = 0
    if ( obj_pyrolysis%enabled ) nchi = obj_pyrolysis%nsub
    ndar = 0
    if ( obj_darcy%enabled ) ndar = 4
    nblocks = size(IOinit%block)
    allocate(IOfield%block(1:nblocks))
    do b = 1, nblocks
      allocate(IOfield%block(b)%mesh(1:3,0:IOinit%block(b)%Ni,0:IOinit%block(b)%Nj,0:IOinit%block(b)%Nk))
      allocate(IOfield%block(b)%vars(1:3+nchi+ndar,1:IOinit%block(b)%Ni,1:IOinit%block(b)%Nj,1:IOinit%block(b)%Nk))
    enddo

    ! Check dimensions IOinit and IOsource are compatible
    if (error_source == 0) then
      if (size(IOsource%block) /= size(IOinit%block)) error_dim = 1
      do b = 1, nblocks
        if (IOsource%block(b)%Ni /= IOinit%block(b)%Ni) error_dim = 1
        if (IOsource%block(b)%Nj /= IOinit%block(b)%Nj) error_dim = 1
        if (IOsource%block(b)%Nk /= IOinit%block(b)%Nk) error_dim = 1
      enddo
    endif
    if (error_dim == 1) then
      obj_io%error_message = "[ERROR] source file dimensions not compatible with init file"
      return
    endif

    ! Fill IOfield
    IOfield%solutiontime = IOinit%solutiontime
    do b = 1, nblocks
      IOfield%block(b)%Ni = IOinit%block(b)%Ni
      IOfield%block(b)%Nj = IOinit%block(b)%Nj
      IOfield%block(b)%Nk = IOinit%block(b)%Nk
      IOfield%block(b)%mesh = IOinit%block(b)%mesh
      ! Variable layout differs between the two input files, because ORION puts
      ! only the CELL-CENTRED variables in %vars (x,y,z go to %mesh):
      !   the ic.tec file : "x" "y" "z" T matID -> vars(1)=T, vars(2)=matID
      !   the st.tec file : "x" "y" "z" qvol   -> vars(1)=qvol  [extent (1:1)]
      ! (note: no globs in these comments -- the build uses -cpp, and a
      !  slash-star sequence would open a C block comment)
      ! so the source is read from index 1, not 3. Indexing IOsource with 3 is
      ! out of bounds; in a RELEASE build it silently linearises to
      ! vars(1,i+2,j,k) -- the source shifted by two cells, plus a read past the
      ! end of the array in the last two cells.
      IOfield%block(b)%vars(1:2,:,:,:) = IOinit%block(b)%vars(1:2,:,:,:)
      if (error_source == 0) then
        IOfield%block(b)%vars(3,:,:,:) = IOsource%block(b)%vars(1,:,:,:)
      else
        IOfield%block(b)%vars(3,:,:,:) = 0.0d0
      endif

      ! Pyrolysis state (plan 06 section 3.6), by column POSITION. A new run's
      ! initial field never carries it (virgin). A restart file written with
      ! pyrolysis on carries the seven legacy columns, chi_1..chi_nsub (one
      ! reaction: chi), rho_s and mdot_pg (nv = 9 + nsub) -- the names are
      ! checked when the reader kept them; a legacy restart file (nv = 7)
      ! resets to virgin with a notice; anything else is refused, because a
      ! file written with a different pyrolysis.dat must not be read
      ! positionally. A type-2 file (plan 07 section 3.8) adds p_pg,
      ! u_pg_1..3, eta_prev, eta_st and T_prev (nv = 16 + nsub).
      if ( obj_pyrolysis%enabled ) then
        if ( obj_sim_param%newrun ) then
          IOfield%block(b)%vars(4:3+nchi,:,:,:) = 0.0d0
        else
          nv = size(IOinit%block(b)%vars, 1)
          if ( nv == 9 + nchi .or. nv == 16 + nchi ) then
            do r = 1, nchi
              if ( allocated(IOinit%varnames) ) then
                if ( size(IOinit%varnames) >= 10 + r ) then
                  ! the writer names a one-component field without a suffix
                  ! (one reaction: 'chi'), several as chi_1..chi_nsub
                  if ( nchi == 1 ) then
                    expected = 'chi'
                  else
                    write(expected,'(A,I0)') 'chi_', r
                  endif
                  if ( trim(adjustl(IOinit%varnames(10+r))) /= trim(expected) ) then
                    write(obj_io%error_message,'(A,I0,A)') "[ERROR] pyrolysis restart: "//trim(obj_io%nameinit)// &
                      " column ", 10+r, " is '"//trim(adjustl(IOinit%varnames(10+r)))//"', expected '"//trim(expected)//"'"
                    return
                  endif
                endif
              endif
              IOfield%block(b)%vars(3+r,:,:,:) = IOinit%block(b)%vars(7+r,:,:,:)
            enddo
          elseif ( nv == 7 ) then
            IOfield%block(b)%vars(4:3+nchi,:,:,:) = 0.0d0
            if ( mpi_is_root .and. b == 1 ) &
              write(*,'(A)') ' [NOTICE] restart file carries no chi_* columns: material state reset to virgin'
          else
            write(obj_io%error_message,'(A,I0,A,I0,A,I0)') "[ERROR] pyrolysis restart: "//trim(obj_io%nameinit)//" has ", nv, &
              " variables, expected 7, ", 9 + nchi, " or ", 16 + nchi
            return
          endif
        endif
      endif

      ! Darcy state (plan 07 section 3.8): on a restart from a type-2 file the
      ! columns p_pg (10 + nsub), eta_prev (14 + nsub), eta_st (15 + nsub) and
      ! T_prev (16 + nsub), names checked; a new run, or a type-1 / legacy
      ! file, starts from p-init, T_prev = T and the equation of state
      ! (negative eta sentinels, resolved by Darcy_Init_Cells).
      if ( obj_darcy%enabled ) then
        nv = size(IOinit%block(b)%vars, 1)
        if ( .not. obj_sim_param%newrun .and. nv == 16 + nchi ) then
          if ( allocated(IOinit%varnames) ) then
            if ( size(IOinit%varnames) >= 19 + nchi ) then
              if ( trim(adjustl(IOinit%varnames(13+nchi))) /= 'p_pg'     .or. &
                   trim(adjustl(IOinit%varnames(17+nchi))) /= 'eta_prev' .or. &
                   trim(adjustl(IOinit%varnames(18+nchi))) /= 'eta_st'   .or. &
                   trim(adjustl(IOinit%varnames(19+nchi))) /= 'T_prev' ) then
                write(obj_io%error_message,'(A)') "[ERROR] Darcy restart: "//trim(obj_io%nameinit)// &
                  " does not carry p_pg, eta_prev, eta_st, T_prev at the type-2 column positions"
                return
              endif
            endif
          endif
          IOfield%block(b)%vars(4+nchi,:,:,:) = IOinit%block(b)%vars(10+nchi,:,:,:)
          IOfield%block(b)%vars(5+nchi,:,:,:) = IOinit%block(b)%vars(16+nchi,:,:,:)
          IOfield%block(b)%vars(6+nchi,:,:,:) = IOinit%block(b)%vars(14+nchi,:,:,:)
          IOfield%block(b)%vars(7+nchi,:,:,:) = IOinit%block(b)%vars(15+nchi,:,:,:)
        else
          IOfield%block(b)%vars(4+nchi,:,:,:) = obj_darcy%p_init
          IOfield%block(b)%vars(5+nchi,:,:,:) = IOfield%block(b)%vars(1,:,:,:)
          IOfield%block(b)%vars(6+nchi,:,:,:) = -1.0d0
          IOfield%block(b)%vars(7+nchi,:,:,:) = -1.0d0
          if ( .not. obj_sim_param%newrun .and. mpi_is_root .and. b == 1 ) &
            write(*,'(A)') ' [NOTICE] restart file carries no Darcy state: gas pressure reset to p-init, gas storage to the equation of state'
        endif
      endif
    enddo

  end subroutine Read_vtk_tec


  !> Output setup
  subroutine Setup_Output_Solution ( IOfield )
    use IR_Precision
    use FUSS_Config_Types_m, only: obj_multigrid, obj_io, obj_pyrolysis
    use FUSS_Global_m
    use FUSS_Parameters_m
    use FUSS_Field_List_m,   only: Field_List_Count, Field_List_Get, field_entry_t
    use Lib_ORION_data
    use strings, only: parse
    implicit none
    type(ORION_data), intent(inout) :: IOfield(obj_multigrid%MGL)
    ! Local
    integer :: b, m, fl, c
    type(field_entry_t) :: e
    character(len=40)   :: colname

    ! IO Variables specification
    obj_io%Ovarnames=' "T" "matID" "qvol" "k" "rho" "cs" "h" '
    obj_io%Onvar = 1 + 6

    ! Pyrolysis (plan 06 section 3.6): the field-list rows with output = .true.
    ! are appended after the seven legacy columns, in list order, ncomp
    ! columns each (chi_1..chi_nsub, rho_s, mdot_pg); "h" carries E and
    ! "k"/"rho"/"cs" the caches. The legacy string stays byte-identical
    ! otherwise, so every existing post-processing script keeps working.
    if ( obj_pyrolysis%enabled ) then
      do fl = 1, Field_List_Count()
        e = Field_List_Get( fl )
        if ( .not. ( e%managed .and. e%output ) ) cycle
        if ( e%ncomp == 1 ) then
          obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "'//trim(e%name)//'" '
          obj_io%Onvar = obj_io%Onvar + 1
        else
          do c = 1, e%ncomp
            write(colname,'(A,A,I0)') trim(e%name), '_', c
            obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "'//trim(colname)//'" '
            obj_io%Onvar = obj_io%Onvar + 1
          enddo
        endif
      enddo
    endif

    do m = 1, obj_multigrid%MGL
      IOfield(m)%vtk%node = .false.
      IOfield(m)%tec%node = .false.
      IOfield(m)%tec%bc = .false.
      if (index(obj_io%sol_format,'ascii')>0) then
        IOfield(m)%tec%extension = '.tec'
        IOfield(m)%tec%format = 'ascii'
        IOfield(m)%vtk%format = 'ascii'
      else
        IOfield(m)%tec%extension = '.szplt'
        IOfield(m)%tec%format = 'binary'
        IOfield(m)%vtk%format = 'binary'
      endif
    end do

    ! Concretize the sol subroutine
    Write_solution => Write_vtk_tec

    do m = 1, obj_multigrid%MGL
      do b = 1, Size(IOfield(m)%block)
        IOfield(m)%block(b)%name = 'Block'//trim(str(.true.,b))
      enddo
    enddo

    ! Reallocate IOfield vars if the number of variables read into the backup file is different from the solution one
    ! The reallocation is performed after reading the ICs and before the first solution is written!
    do m = 1, obj_multigrid%MGL
      if ( obj_io%Onvar /= Size(IOfield(m)%block(1)%vars,1) ) then
        do b = 1, Size ( IOfield(m)%block )
          IOfield(m)%block(b)%name = 'Block'//trim(str(.true.,b))
          deallocate ( IOfield(m)%block(b)%vars )
          allocate( IOfield(m)%block(b)%vars(1:obj_io%Onvar, &
                                             1:IOfield(m)%block(b)%Ni, &
                                             1:IOfield(m)%block(b)%Nj, &
                                             1:IOfield(m)%block(b)%Nk) )
        enddo
      endif
    enddo

  end subroutine Setup_Output_Solution

  subroutine Write_vtk_tec ( domain, IOfield, file )
    use FUSS_Advanced_Types_m
    use FUSS_Config_Types_m, only: obj_io, obj_mesh_motion, obj_pyrolysis, obj_darcy, obj_ablation
    use FUSS_Global_m
    use FUSS_Parameters_m
    use IR_Precision
    use Lib_ORION_data
    use Lib_VTK
    use Lib_Tecplot
    use FUSS_Lib_Solid
    use FUSS_Lib_Pyrolysis,  only: Solid_State, Gas_State, Energy_Of, Energy_Of2
    use FUSS_Field_List_m,   only: Field_List_Count, Field_List_Get, field_entry_t
#ifdef USE_MPI
    use FUSS_Field_List_m,   only: FLD_TO
    use FUSS_Mod_GhostExchange, only: gather_field_to_root
#endif
    use FUSS_Mod_Pyrolysis,  only: Write_Energy_Budget
    use FUSS_Mod_Recession,  only: Write_Ablation_Outputs, Write_Interface_Flux
    use FUSS_Mod_MPI, only: mpi_is_root
    use FUSS_Mod_GhostExchange, only: gather_T_to_root, gather_qvol_to_root, gather_matID_to_root, mpi_io_barrier, &
                                      gather_node_to_root
    implicit none
    type(FUSS_domain_type), intent(inout) :: domain
    type(ORION_data), intent(inout)       :: IOfield
    character(llen), intent(in)           :: file
    ! Local
    character(len=llen) :: path
    character(len=llen) :: localpath_vtk
    integer             :: E_IO, b, i, j, k, m, fl, c, col
    real(8)             :: kappa, rho, cp, h
    real(8)             :: rho_s, cp_s, k_s, tau, eps_pg, W, h_pg, mu, G
    type(field_entry_t) :: e

    ! Gather P from all ranks to root (collective — all ranks must call)
    call gather_T_to_root(domain)
    call gather_matID_to_root(domain)
    call gather_qvol_to_root(domain)
    ! a moving mesh: root writes every block's CURRENT nodes (plan 11 stage 2); this
    ! file is also the restart, which would otherwise carry the setup geometry
    if ( obj_mesh_motion%enabled ) call gather_node_to_root(domain)
#ifdef USE_MPI
    ! pyrolysis (plan 11 stage 3): the field-list output rows -- chi, rho_s, mdot_pg; type 2
    ! p_pg, u_pg, eta_prev, eta_st, T_prev -- and, under type 2, TO (the T_prev column once a
    ! step has run). Root then evaluates the auxiliary columns k, rho, cs, h of every block
    ! from the gathered (T, chi, p_pg, eta_st); the k_s / cp_s caches are not gathered.
    if ( obj_pyrolysis%enabled ) then
      do fl = 1, Field_List_Count()
        e = Field_List_Get( fl )
        if ( e%managed .and. e%output ) call gather_field_to_root( domain, e%code, e%ncomp )
      enddo
      if ( obj_darcy%enabled ) call gather_field_to_root( domain, FLD_TO, 1 )
    endif
#endif

    if (mpi_is_root) then
      path = 'OUTPUT/'

      ! Update IOfield variables with domain primitives and other variables
      do b = 1, size(IOfield%block)

        ! Write the CURRENT node positions, not the ones read at setup.
        !
        ! With a moving mesh the solution would otherwise be plotted, probed and
        ! post-processed against the initial geometry. That is not merely a
        ! cosmetic problem: any external check that pairs the written field with
        ! written cell volumes -- an energy-conservation audit, for instance --
        ! silently compares moved-mesh values against initial-mesh volumes and
        ! reports an imbalance that is entirely an artefact of the output.
        ! Restarts need this too (plan 05 section 4.5).
        if ( obj_mesh_motion%enabled ) then
          do k = 0, IOfield%block(b)%Nk ; do j = 0, IOfield%block(b)%Nj ; do i = 0, IOfield%block(b)%Ni
            IOfield%block(b)%mesh(1,i,j,k) = domain%blk(b)%node(i,j,k)%c(1)
            IOfield%block(b)%mesh(2,i,j,k) = domain%blk(b)%node(i,j,k)%c(2)
            IOfield%block(b)%mesh(3,i,j,k) = domain%blk(b)%node(i,j,k)%c(3)
          enddo; enddo; enddo
        endif

        IOfield%block(b)%vars(1,:,:,:) = domain%blk(b)%T(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        IOfield%block(b)%vars(2,:,:,:) = domain%blk(b)%matID(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        IOfield%block(b)%vars(3,:,:,:) = domain%blk(b)%qvol(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        ! Auxiliary variables
        if ( obj_pyrolysis%enabled ) then
          ! plan 06 section 3.6: k_s, total density rho_s + eps_pg rho_pg, cp_s,
          ! the integrated energy E in the "h" slot, then the list's output rows
          ! type 2 (plan 07): the "rho" slot carries rho_s + G (the stored gas)
          ! and "h" the type-2 E with the cell's p_pg and G = eta_st p_pg
          if ( obj_darcy%enabled ) then
            do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
              m = nint( domain%blk(b)%matID(i,j,k) )
              G = domain%blk(b)%eta_st(i,j,k) * domain%blk(b)%p_pg(i,j,k)
              call Solid_State( m, domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), rho_s, cp_s, k_s, tau, eps_pg )
              IOfield%block(b)%vars(4,i,j,k) = k_s
              IOfield%block(b)%vars(5,i,j,k) = rho_s + G
              IOfield%block(b)%vars(6,i,j,k) = cp_s
              IOfield%block(b)%vars(7,i,j,k) = Energy_Of2( m, domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), &
                                                           domain%blk(b)%p_pg(i,j,k), G )
            enddo; enddo; enddo
          else
          do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
            m = nint( domain%blk(b)%matID(i,j,k) )
            call Solid_State( m, domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), rho_s, cp_s, k_s, tau, eps_pg )
            call Gas_State( m, obj_pyrolysis%p_gas, domain%blk(b)%T(i,j,k), W, h_pg, mu, [b, i, j, k] )
            IOfield%block(b)%vars(4,i,j,k) = k_s
            IOfield%block(b)%vars(5,i,j,k) = rho_s + eps_pg * obj_pyrolysis%p_gas * W / ( R_univ * domain%blk(b)%T(i,j,k) )
            IOfield%block(b)%vars(6,i,j,k) = cp_s
            IOfield%block(b)%vars(7,i,j,k) = Energy_Of( m, domain%blk(b)%T(i,j,k), domain%blk(b)%chi(:,i,j,k), obj_pyrolysis%p_gas )
          enddo; enddo; enddo
          endif
          col = 7
          do fl = 1, Field_List_Count()
            e = Field_List_Get( fl )
            if ( .not. ( e%managed .and. e%output ) ) cycle
            select case ( trim(e%name) )
              case ('chi')
                do c = 1, e%ncomp
                  col = col + 1
                  IOfield%block(b)%vars(col,:,:,:) = domain%blk(b)%chi(c,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
                enddo
              case ('rho_s')
                col = col + 1
                IOfield%block(b)%vars(col,:,:,:) = domain%blk(b)%rho_s(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
              case ('mdot_pg')
                col = col + 1
                IOfield%block(b)%vars(col,:,:,:) = domain%blk(b)%mdot_pg(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
              ! Darcy closure (plan 07 section 3.8)
              case ('p_pg')
                col = col + 1
                IOfield%block(b)%vars(col,:,:,:) = domain%blk(b)%p_pg(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
              case ('u_pg')
                do c = 1, e%ncomp
                  col = col + 1
                  IOfield%block(b)%vars(col,:,:,:) = domain%blk(b)%u_pg(c,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
                enddo
              case ('eta_prev')
                col = col + 1
                IOfield%block(b)%vars(col,:,:,:) = domain%blk(b)%eta_prev(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
              case ('eta_st')
                col = col + 1
                IOfield%block(b)%vars(col,:,:,:) = domain%blk(b)%eta_st(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
              case ('T_prev')
                ! The temperature one step before this file's T: TO at the end of a
                ! step. The in-memory T_prev is one level older there (the T^(n-1) of
                ! the step just taken), and a restart's first step keeps the file's
                ! value as its T^(n-1) (Copy_State skips iter 1), so the file carries
                ! TO. Before the first step TO is unset: the setup value is written.
                col = col + 1
                if ( domain%iter >= 1 ) then
                  IOfield%block(b)%vars(col,:,:,:) = domain%blk(b)%TO(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
                else
                  IOfield%block(b)%vars(col,:,:,:) = domain%blk(b)%T_prev(1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
                endif
              case default
                write(*,'(A)') ' [ERROR] Write_vtk_tec: no handler for the output field '//trim(e%name)
                call fuss_abort ( 1 )
            end select
          enddo
        else
          do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
            call co_Kappa( domain%blk(b)%matID(i,j,k), domain%blk(b)%T(i,j,k), kappa )
            call co_Rho( domain%blk(b)%matID(i,j,k), domain%blk(b)%T(i,j,k), rho )
            call co_Cp( domain%blk(b)%matID(i,j,k), domain%blk(b)%T(i,j,k), cp )
            call co_H( domain%blk(b)%matID(i,j,k), domain%blk(b)%T(i,j,k), h )
            IOfield%block(b)%vars(4,i,j,k) = kappa
            IOfield%block(b)%vars(5,i,j,k) = rho
            IOfield%block(b)%vars(6,i,j,k) = cp
            IOfield%block(b)%vars(7,i,j,k) = h
          enddo; enddo; enddo
        endif
      enddo
    
      ! Write the IOfield accordingly to the solution format
      if (index(obj_io%sol_format,'vtk')>0) then
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        E_IO = vtk_write_structured_multiblock(orion=IOfield,vtspath=trim(localpath_vtk)//trim(file), &
                                                             vtmpath=trim(path)//trim(file),varnames=obj_io%Ovarnames,time=domain%time)
      else
        E_IO = tec_write_structured_multiblock(orion=IOfield,varnames=obj_io%Ovarnames,filename=trim(path)//trim(file)//trim(IOfield%tec%extension))
      end if

      ! Ablating walls (plan 08 sections 2.3 R3, 3.2): seb.txt and the restart state
      if ( obj_ablation%enabled ) call Write_Ablation_Outputs( domain )
      call Write_Interface_Flux( domain )
    end if

    ! Energy budget line at every solution write (plan 06 section 3.6), on every rank:
    ! the block sums travel to root, which writes (plan 11 stage 3)
    if ( obj_pyrolysis%enabled ) call Write_Energy_Budget( domain )

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine Write_vtk_tec

end module FUSS_IO_Solution