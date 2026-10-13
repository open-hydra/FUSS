# Input File

FUSS is configured through a single **INI-format** file called `input.ini`, located in the case root directory.

The file is organised into parameter blocks — each block controls a specific aspect of the simulation (run control, numerics, mesh generation, materials, initial conditions, boundary conditions, output). Parameters not specified take their default values.

## File Structure

```ini
[SECTION-NAME]
parameter = value
```

Example:

```ini
[FUSS-Parameters]
res-threshold = 1e-8
iter-threshold = 500000

[FUSS-Numerics]
time-scheme = RK3
vnn = 2.0
time-accurate = false
integration-variables = cons
irs = true
irs-beta = 0.5
```

!!! warning "Default values"
    If a parameter is not specified, the **default value** is used silently.

!!! warning "Parameter names are case-sensitive"
    Incorrect names or values are ignored and the default is used instead. Check the [Parameter Registry](registry.md) for the exact spelling.

---

## Sections

| Section | Description | Reference |
|---------|-------------|-----------|
| `[FUSS-Parameters]` | Simulation control: convergence thresholds, restart | [→](registry.md#fuss-parameters) |
| `[FUSS-Numerics]` | Time scheme, VNN, IRS, integration variables | [→](registry.md#fuss-numerics) |
| `[FUSS-Multigrid]` | Multigrid acceleration settings | [→](registry.md#fuss-multigrid) |
| `[FUSS-IO]` | Output formats, frequencies, and overwrite policy | [→](registry.md#fuss-io) |
| `[FUSS-Probes]` | Point probes for time-history recording | [→](registry.md#fuss-probes) |
| `[FUSS-Pyrolysis]` | Pyrolysis of a charring solid (its presence switches the model on) | [→](#pyrolysis) |
| `[FUSS-Darcy]` | Gas pressure solver of the type-2 gas closure (read only with `gas-closure = type2`) | [→](#type-2-gas-closure-darcy) |
| `[GRIB-meshgen]` | Built-in mesh generator global settings | [User Guide](using.md#mesh) |
| `[GRIB-Block*]` | Per-block mesh geometry and resolution | [User Guide](using.md#mesh) |
| `[GPB-Phase1]` | Material property definitions | [Initial Conditions](initial-conditions.md#materials) |
| `[ICB-Block*]` | Per-block initial conditions | [Initial Conditions](initial-conditions.md) |
| `[BCB-Block*]` | Per-block face boundary condition assignments | [Boundary Conditions](boundary-conditions.md) |
| Named BC sections | BC type and parameters (temperature, flux, convection, radiation) | [Boundary Conditions](boundary-conditions.md) |

---

## Parameter Registry

The complete list of all parameters, defaults, and allowed values is in the **[Parameter Registry](registry.md)**. It is generated automatically from the source code at every release — do not edit it manually.

---

## Pyrolysis

A section `[FUSS-Pyrolysis]` switches on the pyrolysis model of a charring solid (phase 2 of the porous-pyrolysis programme, plan 06): the solid is a set of phases, each decomposing through Arrhenius reactions, the stored energy `E(T, χ)` is the integrated variable and the pyrolysis gas leaves through one face of the (single) block. Without the section the solver runs the legacy conduction model unchanged, bit for bit.

```ini
[FUSS-Pyrolysis]
gas-closure       = type1        ; type1: quasi-steady 1-D outflow; type2: Darcy pressure equation (below)
outflow-direction = +i           ; type1 only: +i, -i, +j, -j, +k, -k, the gas leaves through that end of the block
p-gas             = 101325.0     ; type1 only: [Pa] pressure at which the gas table is evaluated
hbar-guard-report = false        ; report the materials with rho_v = rho_c at setup
```

Requirements enforced at setup (`error stop 1` with a message): `integration-variables = cons`, `time-accurate = true`, no IRS, one grid level, `law = static`, one block without 101/102/103/104 faces, OpenMP only; under `type1` also an `outflow-direction` and a wall BC (301–305) on the outflow end of the block.

`type1` is exact only for one-dimensional, quasi-steady gas flow with an impermeable back face, and it takes the gas enthalpy leaving the block at the last cell's temperature. `type2` solves the gas momentum (Darcy) and the gas storage instead; see [Type-2 gas closure](#type-2-gas-closure-darcy).

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

The solution file gains the columns `chi_1 … chi_nsub` (reaction advancements, read back on restart), `rho_s` and `mdot_pg` (gas mass flux through the cell's downstream face, kg m⁻² s⁻¹); `k`, `rho`, `cs` carry the blended solid conductivity, the total density and the blended specific heat, and `h` the stored energy `E` [J/m³]. `OUTPUT/<prefix>energy-budget.txt` records, at every solution write, `t E_tot M_s Q_out Q_gas_out M_gas_out Q_src` (the discrete energy and solid-mass balances close to round-off; `test/pyrolysis/energy_audit/verify.py` checks them).

### Type-2 gas closure (Darcy)

`gas-closure = type2` (plan 07) solves the gas pressure `p_pg` of the pyrolysis gas from the gas mass balance with Darcy's law, `∂(ε_g ρ_g)/∂t + ∇·(ρ_g u_g) = π`, `u_g = −(K/(μ ε_g))·∇p` (`ρ_g = p W/(ℜ T)`, `W` and `μ` from `gas.dat`, `K` the blended permeability tensor), once per material step: the conductance `ρ_g K/μ` at the start of the step, the storage coefficient `ε_g W/(ℜT)` extrapolated to its end (limited: never below half its start value) and the gas production to its middle. The pressure update itself is backward Euler, so a transient pressure field is first order in time (exact pressure-diffusion benchmark, `test/pyrolysis/exact/pressure_transient`); `T` and `χ` keep the second order of the 0-D gates. The stored gas mass `ε_g ρ_g` is carried as a state and corrected at the end of the step so that the gas mass balance closes to round-off; when the correction would take a cell's stored gas below half its start value (threshold kinetics), the reactions that ran supply the difference. The stored energy becomes `E = H_s(T, χ) + ε_g (ρ_g h_g − p)`, and the gas enthalpy is advected with the Darcy mass flux. The gas leaves through the **open** faces, where `p` is prescribed, with the enthalpy `h_g(p_w, T_w)` at the wall temperature `T_w` that the face's BC computes (the cell temperature on a face without one, 103); every other face is impermeable.

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

- **Open faces.** By default the convective (303) and radiative-convective (305) faces are open, `p = p_w`; the value line of a 303 or 305 face takes an optional fourth number `p_w` [Pa] (`hconv, qrad, Tref, p_w`; without it `p-wall` applies). Any wall code of `{301, 302, 303, 304, 305, 306, 103}` can be listed in `open-bc-types`.
- **Refused at setup:** `[FUSS-Darcy]` without `gas-closure = type2`; `p-gas` with `type2`; `pseudo-vnn` outside (0, 1]; an unknown `open-bc-types` token; `axes-follow-grid = true` (not built); a 102 or 104 (chimera) face; no open face while gas is produced (kinetics on, or `pi-source > 0`); no open face with `quasi-steady = true`; a material without `K_v`/`K_c`; a `p_w ≤ 0` in `bc.txt`.
- **Solver.** `tridiag` solves each grid line along the direction of a block's open faces exactly (one sweep per step on a one-dimensional block); a block whose open faces lie in several directions, or none, falls back to `pseudo`. `max-iter` sweeps without convergence stop the run (`error stop 3`). `OUTPUT/<prefix>darcy.txt` records per step `iter time block sweeps resid dp_last`.
- **Quasi-steady.** With no stored gas, the difference between the gas the pressure solve expels (production predicted at the mid-step temperature) and the gas the RK stages release waits as a carry, `(eta_st − eta_prev)·p_pg` (restart columns), and joins the next step's source (never below zero): the gas mass is exact with a one-step lag. Before 2026-10-13 the difference was dropped, and threshold kinetics lost 5–11 % of the gas (Baum & Atreya benchmark, `test/pyrolysis/lit/baum_atreya2007`).
- **Output.** The solution file gains `p_pg`, the cell gas velocity `u_pg_1 u_pg_2 u_pg_3` [m/s] and the state of the scheme `eta_prev`, `eta_st`, `T_prev` (read back on restart with `p_pg`; a type-2 run restarted from a file without them, such as a type-1 solution, starts the gas at `p-init` and the equation of state with a NOTICE, and its first step has no extrapolation history); `mdot_pg` is the magnitude of the cell's Darcy mass flux `|ρ_g (K/μ)·∇p|` [kg m⁻² s⁻¹] and `u_pg` the interstitial velocity. `energy-budget.txt` gains an eighth column `M_g`, the stored gas mass. A probe may record `variables = T p_pg`.
