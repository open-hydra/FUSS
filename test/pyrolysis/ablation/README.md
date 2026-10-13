# Ablating wall (BC 306) cases — plan 08

Plan 08 builds the ablating wall in four stages; stage A (2026-10-13) is the wall on a static mesh. Its cases:

| case | what it checks | gate | runtime (-p 1) |
|---|---|---|---|
| `atc2off_306_static`, `atc2off_306_static_t2` (`.slow`) | 3.1 (a): a 306 face with the chemistry off and zero optics hands the face to the BC 303 kernel. Every output except the wall files must be byte-identical to the unmodified 303 case, under type 1 and under type 2. | `run_physics_gates.sh -c ablation_306_static` | 77 s |
| `zero_blowing` | 3.1 (b): with λ = 0, β = 1 bitwise and the convective term is `ρ_e u_e C_H0 (h_r − h_w)` to 1e-15 | `-c ablation_zero_blowing` | 1 s |
| `seb_closure` | 3.2: the surface energy balance closes on every row (chemistry on, the ramp, chemistry off), and β sits at the Picard's closed-form fixed point; it is also the 306 case of `run_restart_gate.sh` | `-c ablation_seb_closure` | 1 s |
| `ATC21` (`.slow`) | 3.3: Ablation Test Case 2.1 against PATO (static mesh) | `-c ablation_atc21` | 45 s |
| `demo_base`, `demo_base_pyro` | the bases of the ablating-wall guard demos (`run_guard_demos.sh`) | — | 1 s |

Generators:
- `atc2x306` and `atc2off --bc306-chemoff` of `../tools/make_case.py` (the command is in each case's README);
- `../tools/make_abl_demo_base.py` for the demo bases.

Shared verify code is `../tools/abl_tools.py` (seb.txt, env tables, §3.2 closure, PATO readers, PATO's front rule).

Other material in this directory:
- `tables/`: the equilibrium B′ table and its provenance;
- `edit_pato_atc2x.py`: the PATO oracle edits for ATC 2.1, 2.2 and 2.3.

The face writes two files of its own:
- `OUTPUT/seb.txt`: one row per 306 face at every solution output, with every term of the balance;
- `OUTPUT/bc306-state.dat`: the restart state, 17 digits.

`wall.tec` has 7 columns when a 306 face exists: `Tw qw sdot mdot_pg mdot_ch Bprime recession`.
