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
gas-closure       = type1        ; quasi-steady 1-D outflow (type2 = Darcy, plan 07)
outflow-direction = +i           ; +i, -i, +j, -j, +k, -k: the gas leaves through that end of the block
p-gas             = 101325.0     ; [Pa] pressure at which the gas table is evaluated (type1)
hbar-guard-report = false        ; report the materials with rho_v = rho_c at setup
```

Requirements enforced at setup (`error stop 1` with a message): `integration-variables = cons`, `time-accurate = true`, no IRS, one grid level, `law = static`, one block without 101/102/103 faces, a wall BC (301–305) on the outflow end of the block, OpenMP only.

### Material files

Besides `INPUT/<prefix>properties.dat` (still loaded, not used by the pyrolysis path), every material of `phase.txt` needs four more files in `INPUT/`. Until ATLAS can write them, `test/pyrolysis/tools/make_tacot_tables.py --pato-dir <PATO 3.1 root> --out INPUT` writes them from PATO's TACOT data, and `make_case.py` builds complete cases.

| file | content |
|------|---------|
| `<prefix>pyrolysis.dat` | plain text, `//` comments, one `material <name>` block per material in `phase.txt` order, `key value` lines with PATO's `constantProperties` names: `nSolidPhases`, `rhoI[i]`, `epsI[i]`, `nPyroReac[i]`, `eps_g_v`, `eps_g_c`, and per reaction `F[i][j] A[i][j] E[i][j] m[i][j] n[i][j] T[i][j] h[i][j]` (`h` is parsed and ignored: the decomposition heat is in the enthalpy tables). Each key exactly once per block. |
| `<prefix>virgin.dat` | Tecplot POINT table, one zone per material: `T cp h ki kj kk emissivity absorptivity`, integer `T` rows 1 K apart (TACOT: 256–4000 K); `h` absolute, `h(298.15 K)` = the enthalpy of formation; `ki = kj = kk` (orthotropic conductivity is not supported) |
| `<prefix>char.dat` | same layout and the same `T` range as `virgin.dat` |
| `<prefix>gas.dat` | one zone per (material, pressure level), material-major, `p` ascending: `T p M h mu` (`M` kg/mol, `h` J/kg, `mu` Pa·s), integer `T` rows 1 K apart, the same range in every zone of a material (TACOT: 200–3975 K, five levels) |

Tecplot layout of the three tables: `TITLE = "FUSS pyrolysis table"`, `VARIABLES = "T" ...`, then per zone `ZONE T="<material>"` (gas: `ZONE T="<material> p=<p>"`), `I=<rows>, F=POINT` and the rows; nothing after the last row.

### Output

The solution file gains the columns `chi_1 … chi_nsub` (reaction advancements, read back on restart), `rho_s` and `mdot_pg` (gas mass flux through the cell's downstream face, kg m⁻² s⁻¹); `k`, `rho`, `cs` carry the blended solid conductivity, the total density and the blended specific heat, and `h` the stored energy `E` [J/m³]. `OUTPUT/<prefix>energy-budget.txt` records, at every solution write, `t E_tot M_s Q_out Q_gas_out M_gas_out Q_src` (the discrete energy and solid-mass balances close to round-off; `test/pyrolysis/energy_audit/verify.py` checks them).

### Limitation of the type-1 closure

The gas transport of `gas-closure = type1` is exact only for one-dimensional, quasi-steady flow with an impermeable back face; in two- and three-dimensional configurations it is an approximation, because the direction and the velocity of the gas would have to come from a momentum equation (the Darcy closure of plan 07). The enthalpy leaving with the gas is evaluated at the last cell's temperature, not at the wall.
