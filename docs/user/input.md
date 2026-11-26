# Input File

FUSS reads its settings from `input.ini` in the directory it runs in, with the FiNeR parser.

- With MPI, rank 0 alone reads the file and passes its text to every rank, at set-up and at every reading during the run. Every rank takes the same values.
- A file that cannot be opened, is empty or has no section stops every rank with exit code 1.

## Syntax

```ini
[FUSS-Parameters]
res-threshold = 1e-8

[FUSS-Numerics]
time-scheme = RK3
vnn = 1.0               ; inline comment after a value
time-accurate = false
integration-variables = cons
irs-beta = 0.5
```

- A section starts with `[name]`; each option is `key = value`.
- Lines starting with `!`, `;` or `#` are comments, and `;` starts an inline comment after a value.
- A comment line may contain `=`.
- A `key = value` line that cannot be read (after a malformed section header, or with no value) is reported at set-up with its line number, for example `[WARNING] input.ini, line 7: key of [section] is not read ...`. The parameter keeps its default.

---

## Sections Read by FUSS

| Section | Content | Keys |
|---------|---------|------|
| `[FUSS-Parameters]` | New run or restart, stop criteria | [registry](registry.md#fuss-parameters) |
| `[FUSS-Numerics]` | Time scheme, VNN, IRS, integration variables | [registry](registry.md#fuss-numerics) |
| `[FUSS-Multigrid]` | Grid levels and iterations per level | [registry](registry.md#fuss-multigrid) |
| `[FUSS-IO]` | Formats and frequencies of the output | [registry](registry.md#fuss-io) |
| `[FUSS-Probes]` and one section per probe | Probes | [registry](registry.md#fuss-probes), [Running FUSS](using.md#probes) |
| `[FUSS-Pyrolysis]` | Pyrolysis of a charring solid (its presence switches the model on) | [registry](registry.md#fuss-pyrolysis), [Pyrolysis](#pyrolysis) |
| `[FUSS-Ablation]` | Surface energy balance of the ablating walls 306, 307 | [registry](registry.md#fuss-ablation), [Ablating wall](#ablating-wall-306-307) |
| `[FUSS-Darcy]` | Gas pressure solver of the type-2 gas closure (read only with `gas-closure = type2`) | [registry](registry.md#fuss-darcy), [Type-2 gas closure](#type-2-gas-closure-darcy) |

The [Parameter Registry](registry.md) lists every key with its default and allowed values. It is generated from the source by `bin/DocGen` (`src/app/docgen.f90`): do not edit it by hand. Keys that come in numbered series appear there once as examples: `probe<n>` with its section (`probe-placeholder`), and `level<m>-iter` for two levels.

One key is not in the registry: `levels` in `[FUSS-Multigrid]`, the number of grid levels (default 1, no multigrid). FUSS reads it before building the registry, to know how many `level<m>-iter` keys to expect. Multigrid is described in [Time Integration](../theory/time-integration.md#multigrid-grid-sequencing).

The other sections of the test-case files are not read by FUSS; they are the input of the pre-processors that write `INPUT/`:

| Section | Pre-processor | Described in |
|---------|---------------|--------------|
| `[GRIB-meshgen]`, `[GRIB-Block<n>]` | GRIB (mesh) | [Running FUSS](using.md#mesh) |
| `[GPB-Phase1]` | ATLAS GPB (`phase.txt`, `properties.dat`) | [Initial Conditions](initial-conditions.md#materials) |
| `[ICB-Block<n>]` | ATLAS ICB (`ic.tec`) | [Initial Conditions](initial-conditions.md#initial-field-inputic) |
| `[STB-Block<n>]` | ATLAS STB (`st.tec`) | [Initial Conditions](initial-conditions.md#heat-source-inputst) |
| `[ATLAS-Parameters]`, `[BCB-Block<n>]`, named BC sections | ATLAS BCB (`bc.txt`) | [Boundary Conditions](boundary-conditions.md#atlas-input) |

---

## Defaults and Validation

- A key that is not set takes its default, except `time-accurate`, which is required: without it the run stops at set-up.
- In a FUSS section an unknown (misspelt) key, a repeated section, an invalid value or a string longer than the variable that holds it is refused at setup with the key named. Check the spelling against the [registry](registry.md).
- Every rule of the registry is checked at set-up: bounds and intervals on numbers (on each value of a list), the allowed names of a text, `true` or `false` for a logical. A value out of its rule stops the run with `[ERROR] <key> must be ...` and exit code 1.
- Text values are case-sensitive: `time-scheme = rk3` is rejected.
- A value that does not read as the type of its key stops the run, naming the file, section, key and value. Examples: `sol-diter = 2.5` or `1e4` (an integer is a plain digit string), `sol-overwrite = yes`, `vnn = abc`, a list with too many values, `vnn = NaN` or `1e999`.
- Probes must be numbered without gaps: `probe3` without `probe2` stops the run.
- IRS is switched on by `irs-beta > 0`.
- With more than one multigrid level `iter-threshold` has no effect: each level stops at its `level<m>-iter`. A run that sets it anyway gets a `WARNING` line under `Input file OK`.

## Changes During a Run

FUSS reads `input.ini` again every `ini-diter` iterations (default 10000); every rank reads the text that rank 0 passes, so every rank takes the same decision. A reading without error changes only the runtime-mutable keys (`res-threshold`, `time-threshold`, `iter-threshold`, `sol-diter`, `sol-dtime`, `res-diter`, `shell-diter`, `ini-diter`: `RUNTIME_MUTABLE` in `src/lib/config/Backend_INI.f90`). Any other changed key keeps its set-up value and is listed in `OUTPUT/<prefix>runtime-ini-ignored.txt`, `vnn`, `level<m>-iter`, `irs-beta` and the multigrid keys among them.

- IRS can be neither switched on nor off during a run (see [Time Integration](../theory/time-integration.md#implicit-residual-smoothing-irs)).
- The values read at runtime must keep to the rules of the registry, as at set-up.

`ini-error` in `[FUSS-IO]` decides what an error in a runtime reading does:

| `ini-error` | Effect |
|---|---|
| `warn` (default) | The whole reading is discarded, including its valid keys. One warning gives the iteration and the error (`Input file, read again NOT TAKEN`); the run goes on with its previous values, and `OUTPUT/<prefix>runtime-ini-ignored.txt` gets the line `reading not taken: <error>`. |
| `stop` | Every rank stops with the message and exit code 1, as at set-up. |

An error of the reading is a file that cannot be read, or a value out of its rule, that does not read as its type or that is not finite. The value of `ini-error` in force when the reading starts decides.

## Pyrolysis

A section `[FUSS-Pyrolysis]` switches on the pyrolysis model of a charring solid: the solid is a set of phases, each decomposing through Arrhenius reactions, the stored energy `E(T, χ)` is the integrated variable and the pyrolysis gas leaves through one face of the (single) block. Without the section the solver runs the legacy conduction model unchanged, bit for bit.

```ini
[FUSS-Pyrolysis]
gas-closure       = type1        ; type1: quasi-steady 1-D outflow; type2: Darcy pressure equation (below)
outflow-direction = +i           ; type1 only: +i, -i, +j, -j, +k, -k, the gas leaves through that end of the block
p-gas             = 101325.0     ; type1 only: [Pa] pressure at which the gas table is evaluated
hbar-guard-report = false        ; report the materials with rho_v = rho_c at setup
```

Requirements enforced at setup (`error stop 1` with a message): `integration-variables = cons`, `time-accurate = true`, no IRS, one grid level, no chimera (102, 104) or coupled (103) face; under `type1` also one block without block connections (101) or periodic faces (201: its outflow line integral does not cross block faces), an `outflow-direction` and a wall BC (301–307) on the outflow end of the block. Under `type2` blocks may be connected (101): the reaction state, the gas pressure and both conductivities cross the connection, and its gas mass and enthalpy flows are conservative (the two blocks evaluate one face flux).

`type1` is exact only for one-dimensional, quasi-steady gas flow with an impermeable back face. The gas crosses an interior face with the enthalpy at the mean temperature of the two cells and leaves the block at the wall temperature of the outflow face (second order in space). `type2` solves the gas momentum (Darcy) and the gas storage instead; see [Type-2 gas closure](#type-2-gas-closure-darcy).

### Material files

Besides `INPUT/<prefix>properties.dat` (still loaded, not used by the pyrolysis path), every material of `phase.txt` needs four more files in `INPUT/`. Until ATLAS can write them, `test/pyrolysis/tools/make_tacot_tables.py --pato-dir <PATO 3.1 root> --out INPUT` writes them from PATO's TACOT data, and `make_case.py` builds complete cases.

| file | content |
|------|---------|
| `<prefix>pyrolysis.dat` | plain text, `//` comments, one `material <name>` block per material in `phase.txt` order, `key value` lines with PATO's `constantProperties` names: `nSolidPhases`, `rhoI[i]`, `epsI[i]`, `nPyroReac[i]`, `eps_g_v`, `eps_g_c`, and per reaction `F[i][j] A[i][j] E[i][j] m[i][j] n[i][j] T[i][j] h[i][j]` (`h` is parsed and ignored: the decomposition heat is in the enthalpy tables). Each key exactly once per block. Type 2 also needs the permeability tensors `K_v` and `K_c` [m²] of the virgin and the charred material, nine reals each, row-major, in the mesh frame (symmetric, positive diagonal; an isotropic `K` is `K 0 0 0 K 0 0 0 K`), blended with the charring degree like the conductivity; `B_v`, `B_c` (same layout) are read when `[FUSS-Darcy] klinkenberg = true`. Type 1 ignores the four keys. |
| `<prefix>virgin.dat` | Tecplot POINT table, one zone per material: `T cp h ki kj kk emissivity absorptivity`, integer `T` rows 1 K apart (TACOT: 256–4000 K); `h` absolute, `h(298.15 K)` = the enthalpy of formation; `ki = kj = kk` (orthotropic conductivity is not supported) |
| `<prefix>char.dat` | same layout and the same `T` range as `virgin.dat` |
| `<prefix>gas.dat` | one zone per (material, pressure level), material-major, `p` ascending: `T p M h mu` (`M` kg/mol, `h` J/kg, `mu` Pa·s), integer `T` rows 1 K apart, the same range in every zone of a material (TACOT: 200–3975 K, five levels) |

Tecplot layout of the three tables: `TITLE = "FUSS pyrolysis table"`, `VARIABLES = "T" ...`, then per zone `ZONE T="<material>"` (gas: `ZONE T="<material> p=<p>"`), `I=<rows>, F=POINT` and the rows; nothing after the last row.

### Output

The solution file gains the columns `chi_1 … chi_nsub` (reaction advancements, read back on restart; `chi` when the material has one reaction), `rho_s` and `mdot_pg` (gas mass flux through the cell's downstream face, kg m⁻² s⁻¹); `k`, `rho`, `cs` carry the blended solid conductivity, the total density and the blended specific heat, and `h` the stored energy `E` [J/m³]. `OUTPUT/<prefix>energy-budget.txt` records, at every solution write, `t E_tot M_s Q_out Q_gas_out M_gas_out Q_src` (the discrete energy and solid-mass balances close to round-off; `test/pyrolysis/energy_audit/verify.py` checks them).

### Type-2 gas closure (Darcy)

`gas-closure = type2` solves the gas pressure `p_pg` of the pyrolysis gas from the gas mass balance with Darcy's law, `∂(ε_g ρ_g)/∂t + ∇·(ρ_g u_g) = π`, `u_g = −(K/(μ ε_g))·∇p` (`ρ_g = p W/(ℜ T)`, `W` and `μ` from `gas.dat`, `K` the blended permeability tensor), once per material step: the conductance `ρ_g K/μ` at the start of the step, the storage coefficient `ε_g W/(ℜT)` extrapolated to its end (an extrapolation below half the start value is replaced by the start value) and the gas production to its middle. The pressure update itself is backward Euler, so a transient pressure field is first order in time (exact pressure-diffusion benchmark, `test/pyrolysis/exact/pressure_transient`); `T` and `χ` keep the second order of the 0-D gates. The stored gas mass `ε_g ρ_g` is carried as a state and corrected at the end of the step so that the gas mass balance closes to round-off; when the correction would empty a cell's stored gas (threshold kinetics), the reactions that ran supply the difference, up to a tenth of the start value. The stored energy becomes `E = H_s(T, χ) + ε_g (ρ_g h_g − p)`, and the gas enthalpy is advected with the Darcy mass flux. The gas leaves through the **open** faces, where `p` is prescribed, with the enthalpy `h_g(p_w, T_w)` at the wall temperature `T_w` that the face's BC computes (the cell temperature on a face without one, 103); every other face is impermeable.

```ini
[FUSS-Pyrolysis]
gas-closure = type2

[FUSS-Darcy]
solver        = tridiag       ; tridiag (line-implicit along the open-face direction) | pseudo (damped Jacobi)
tol           = 1.0e-8        ; sweep stop: max cell mass residual <= tol * the largest flux/source/storage term
p-init        = 101325.0      ; [Pa] uniform initial gas pressure (new run)
p-wall        = 101325.0      ; [Pa] pressure on an open face whose bc.txt line gives none
open-bc-types = 303 305 306   ; BC codes whose faces are open to the gas
quasi-steady  = false         ; true: no gas storage (neither in the mass balance nor in E)
```

All keys and defaults: [Parameter Registry](registry.md#fuss-darcy).

- **Open faces.** By default the convective (303) and radiative-convective (305) faces are open, `p = p_w`; the value line of a 303 or 305 face takes an optional fourth number `p_w` [Pa] (`hconv, qrad, Tref, p_w`; without it `p-wall` applies). Any wall code of `{301, 302, 303, 304, 305, 306, 103}` can be listed in `open-bc-types`. An ablating wall (306, 307) is open whatever the list holds, at the `p_w` of its env table or exchange.
- **Refused at setup:** `[FUSS-Darcy]` without `gas-closure = type2`; `p-gas` with `type2`; `pseudo-vnn` outside (0, 1]; an unknown `open-bc-types` token; `axes-follow-grid = true` (not built); a 102 or 104 (chimera) face; a periodic face (201) whose rotation $R$ does not leave a material's `K_v`, `K_c` (and `B_v`, `B_c` with `klinkenberg = true`) invariant, $\lVert R K R^T - K\rVert > 2\cdot 10^{-3}\lVert K\rVert$ (the tensors are in the global axes and the cell across the pair is the image of the source cell: the copies of a periodic domain would differ; an isotropic tensor, or one transversely isotropic about the axis of the pair, passes); no open face while gas is produced (kinetics on, or `pi-source > 0`); no open face with `quasi-steady = true`; a material without `K_v`/`K_c`; a `p_w ≤ 0` in `bc.txt`.
- **Solver.** `tridiag` solves each grid line along the direction of a block's open faces exactly. A block without an
  open face whose 101 face ends the lines of a line-implicit block continues them.
  - The lines of a block, and the blocks across a connection, are coupled by line Jacobi with the conductance across the
    lines on the diagonal.
  - After each sweep, an additive correction removes the error that is uniform across the lines. Line-implicit blocks
    joined through 101 faces, in a row along their lines or side by side, share one such correction, so a domain cut
    into blocks converges as the uncut one (to round-off). Lines that carry the same data (a one-dimensional case cut
    into identical lines, in one block or several) are therefore solved exactly in one sweep.
  - A two-dimensional block takes more sweeps (hundreds on a 40 × 20 slab).
  - A block whose open faces lie in several directions, or none, falls back to `pseudo`. `max-iter` sweeps without convergence stop the run (exit code 1). `OUTPUT/<prefix>darcy.txt` records per step `iter time block sweeps resid dp_last`.
- **Quasi-steady.** With no stored gas, the difference between the gas the pressure solve expels (production predicted at the mid-step temperature) and the gas the RK stages release waits as a carry, `(eta_st − eta_prev)·p_pg` (restart columns), and joins the next step's source (never below zero): the gas mass is exact with a one-step lag.
- **Output.** The solution file gains `p_pg`, the cell gas velocity `u_pg_1 u_pg_2 u_pg_3` [m/s] and the state of the scheme `eta_prev`, `eta_st`, `T_prev` (read back on restart with `p_pg`; a type-2 run restarted from a file without them, such as a type-1 solution, starts the gas at `p-init` and the equation of state with a NOTICE, and its first step has no extrapolation history); `mdot_pg` is the magnitude of the cell's Darcy mass flux `|ρ_g (K/μ)·∇p|` [kg m⁻² s⁻¹] and `u_pg` the interstitial velocity. `energy-budget.txt` gains an eighth column `M_g`, the stored gas mass. A probe may record `variables = T p_pg`.

## Ablating wall (306, 307)

A face of type 306 solves the surface energy balance of the material response for its wall temperature. The balance
covers conduction, the convection of the boundary layer blocked by the blowing gas, the advection of the pyrolysis gas and
of the char, absorbed and emitted radiation, and the B′ table of the surface chemistry (thesis eqs. 3.46–3.50,
PATO's `BprimeBoundaryConditions`). Type 307 is the same balance on a face coupled to MOSE, with its boundary-layer
values from hydra. For the `bc.txt` lines, the env table and the B′ table, see
[Boundary Conditions](boundary-conditions.md#ablating-wall-306-and-coupled-ablating-wall-307).

```ini
[FUSS-Ablation]
seb-tol             = 1.0e-8   ; stop test of the surface energy balance
seb-max-iter        = 50       ; reached -> the run stops (exit code 1)
bprime-out-of-range = clamp    ; clamp | abort
```

All keys: [Parameter Registry](registry.md#fuss-ablation).

- **Blowing.** The face blows its pyrolysis gas and consumes no char (`B′_c = ṁ_c = 0`, the wall enthalpy still from
  the table).
- **Coupled ablating wall (307).** It runs with `lambda = 0`, its values written by hydra before every
  step. A standalone run reads them from the TEST ONLY `coupled-env-file`, evaluated as hydra's are: once per step at
  the step's end time (and once at setup at the start time), never per Runge–Kutta stage. The file holds
  - one table `t p_w rhoeUeCH h_r q_rad` for every 307 face (linear in `t`, clamped; `p_w > 0` on every row), or
  - one block per 307 face (hydra's `[HYDRA-Ablation] env-dump` replayed): a line `// face b i j k f`, the face's
    `bc.txt` indices as in `bc306-state.dat`, opens the block, and the rows below it are that face's table. Every 307
    face has exactly one block, in any order; other `//` lines are comments. A block's first row may hold `p_w = 0`
    (no exchange values yet: at setup the face keeps `p-wall`, as under hydra); a step before the block's second row
    stops the run.
- **Output.**
  - `wall.tec` carries seven columns in a run with a 306/307 face (two otherwise): `T_w`, the conducted flux `q_w`,
    `ṡ`, `ṁ_pg`, `ṁ_ch`, `B′` and `s_total` (`ṡ` and `s_total` are zero).
  - `OUTPUT/<prefix>seb.txt` holds every term of each face's last balance, at every solution output.
  - `OUTPUT/<prefix>bc306-state.dat` holds each face's restart state (wall temperature, blowing correction,
    lagged B′_pg). It is rewritten at every output, and a restart needs it.
- **Refused at setup:**
  - multigrid;
  - a chemistry-on env row without `[FUSS-Pyrolysis]` or without a B′ table;
  - a 307 face without pyrolysis or without a B′ table, or with a `lambda`;
  - `coupled-env-file` without a 307 face; with `// face` blocks, a 307 face without a block, a block that names no 307
    face or a face already named, a block of fewer than two rows, or a row above the first `// face` line.
- **TEST ONLY keys:**
  - `no-char`: makes `B′_c = ṁ_c = 0`;
  - `coupled-env-file`: 307 standalone (one table, or one block per face);
  - `dump-interface-flux`: writes `OUTPUT/<prefix>interface-flux.txt`, the conduction, gas-mass and gas-enthalpy flows
    through every block-connection face;
  - `perturb-interface-node`: moves one interface node, to show the interface check.
- **Verification** (`test/pyrolysis/ablation/`):
  - Ablation Test Case 2.1 against PATO;
  - the closed-form balance and blowing checks;
  - the 307/306 identity and the per-face blocks of `coupled-env-file`;
  - the restart cases.
