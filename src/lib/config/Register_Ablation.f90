!> Registry entries of the [FUSS-Ablation] section (plan 08 section 2.10).
!>
!> The section tunes the ablating wall BC 306: the surface energy balance is
!> solved per face by a safeguarded Newton on the wall temperature with a
!> Picard update of the blowing correction, and the B' table lookup clamps
!> (PATO's behaviour) unless told to stop. The keys are registered
!> unconditionally, as Register_Darcy does, so an absent section leaves every
!> one of them at its default. Every key is immutable mid-run (none is in
!> RUNTIME_MUTABLE). The per-face parameters (lambda, T_inf, the B' table, the
!> env table) are bc.txt's, not this section's.
module FUSS_Register_Ablation
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use FUSS_Config_Types_m
  use FUSS_Input_Registry

  implicit none
  private
  public :: Register_Ablation

contains

  subroutine Register_Ablation ()
    use FUSS_Parameters_m
    implicit none
    character(len=llen) :: section

    section = trim(codename)//'-Ablation'

    call reg%add( trim(section), 'seb-tol', obj_ablation%seb_tol, '1.0e-8', &
                  'Stop test of the surface energy balance of BC 306: |dT_w| <= seb-tol*T_w, |F| <= seb-tol*max|term| '// &
                  'and |d(beta)| <= seb-tol', &
                  '> 0', .false. )

    call reg%add( trim(section), 'seb-max-iter', obj_ablation%seb_max_iter, '50', &
                  'Iteration cap of the surface energy balance; reached -> error stop 3', &
                  '> 0', .false. )

    call reg%add( trim(section), 'no-char', obj_ablation%no_char, '.false.', &
                  'TEST ONLY: B''_ch = mdot_ch = 0 on every 306 face; h_w is still read from the B'' table', &
                  'logical', .false. )

    call reg%add( trim(section), 'bprime-out-of-range', obj_ablation%bprime_out_of_range, 'clamp', &
                  'A B'' table query outside the table: clamp (PATO''s behaviour, counted in seb.txt) or abort '// &
                  '(error stop 3 at the first p or B''g outside)', &
                  'clamp ,  abort', .false. )

    call reg%add( trim(section), 'min-thickness', obj_ablation%min_thickness, '0.0', &
                  'law = nodeshift: a receding line shorter than this stops the run (error stop 2) [m]; 0 = 10 x the '// &
                  'smallest initial cell on the receding lines, at most half the shortest line, resolved at a new '// &
                  'run''s setup and kept by a restart', &
                  '>= 0', .false. )

    call reg%add( trim(section), 'sdot-override', obj_ablation%sdot_override, '-1.0', &
                  'TEST ONLY: the recession rate of every 306 face [m/s] (kinematic recession); < 0 = off', &
                  'real', .false. )

    call reg%add( trim(section), 'dump-interface-flux', obj_ablation%dump_interface_flux, '.false.', &
                  'TEST ONLY: write the flux of every 101 face to OUTPUT/<prefix>interface-flux.txt at each output', &
                  'logical', .false. )

    call reg%add( trim(section), 'perturb-interface-node', obj_ablation%perturb_interface_node, '.false.', &
                  'TEST ONLY: at step 2 move the first shared 101 node of the lowest receding block by 1e-12 m in x '// &
                  '(the interface-node assertion must stop the run)', &
                  'logical', .false. )

    call reg%add( trim(section), 'smoothing-width', obj_ablation%smoothing_width, '5', &
                  'law = nodeshift: nodes of the stencil that averages the surface recession rate along each surface '// &
                  'grid line (thesis Eq. 5.1, truncated at a line end; odd; 1 = no smoothing)', &
                  '>= 1', .false. )

    call reg%add( trim(section), 'coupled-env-file', obj_ablation%coupled_env_file, 'none', &
                  'TEST ONLY: every coupled ablating wall (307) takes its exchange values from this INPUT/ file '// &
                  '(rows t p_w rhoeUeCH h_r q_rad, linear in t) instead of hydra; lines // face b i j k f split it '// &
                  'into one block per 307 face (a replay of hydra''s env-dump); none = from hydra', &
                  '', .false. )

  end subroutine Register_Ablation

end module FUSS_Register_Ablation
