# Pyrolysis Kinetics and 0-D Decomposition

Seven cases isolate the decomposition kinetics and its energy coupling from conduction: a slab held
or driven uniform in space, so that the solver's update can be checked cell by cell against a
closed-form or numerically-exact ordinary differential equation (ODE). Together they verify the
reaction-advancement propagator, the type-1 (`prim`) and type-2 (Darcy) energy splittings, and the
line integral of the pyrolysis gas flux, independently of the spatial discretisation.

Setup and run: `test/pyrolysis/kinetics_isothermal/README.md`, `test/pyrolysis/adiabatic0d/README.md`,
`test/pyrolysis/heated0d/README.md` (the type-1 cases); the type-2 twins
`test/pyrolysis/darcy/adiabatic0d_type2/` and `darcy/heated0d_type2/` are run and read by the same
`verify.py` as their type-1 source, see `test/pyrolysis/darcy/README.md`.

---

## Material model

Every case decomposes a two-phase material (TACOT) with three pyrolysis reactions in the second
phase, read from `INPUT/pyrolysis.dat`:

| Reaction | $F_r$ | $A_r$ [1/s] | $E_r$ [J/mol] | $m_r$ | $T_r$ threshold [K] |
|---|---:|---:|---:|---:|---:|
| 1 | 0.25 | $1.2\times 10^{4}$ | 71130.89 | 3 | 333.3 |
| 2 | 0.19 | $4.97777\times 10^{8}$ | 169975 | 3 | 555.6 |
| 3 | 0.06 | $4.97777\times 10^{8}$ | 169975 | 3 | 555.6 |

with phase densities $\rho_{I,1}=1600$, $\epsilon_{I,1}=0.1$ (inert), $\rho_{I,2}=1200$,
$\epsilon_{I,2}=0.1$. The reaction advancement obeys
$\mathrm{d}\chi_r/\mathrm{d}t = \kappa_r(T)(1-\chi_r)^{m_r}$ with
$\kappa_r = A_r T^{n_r}\exp(-E_r/RT)$ above $T_r$, $0$ below. The `kinetics_isothermal` cases use
constant-property (`--isothermal`) solid/gas tables ($c_p$, $k$, $h$ linear and equal between gas and
solid); `adiabatic0d`, `heated0d` and their type-2 twins use the full TACOT property tables.

## Problem and reference

| Case | Faces | $T(0)$ | Source | Oracle |
|---|---|---|---|---|
| `kinetics_isothermal/T0700`, `T0900`, `T1200` | both BC 302 at $T_0$ | 700, 900, 1200 K | none | closed-form propagator (below) |
| `adiabatic0d/T1000` | both BC 301 (adiabatic) | 1000 K | none | `tools/ode0d.py`, type-1 variant |
| `heated0d` | both BC 301 | 300 K | $\dot q=3.75\times10^{8}$ W/m³ | `tools/ode0d.py`, type-1 variant |
| `darcy/adiabatic0d_type2` | both BC 301, open at `p-wall` ($K=10^{-7}$ m²) | 1000 K | none | `ode0d.py`, type-2 variant |
| `darcy/heated0d_type2` | both BC 301, open at `p-wall` ($K=10^{-7}$ m²) | 300 K | $\dot q=3.75\times10^{8}$ W/m³ | `ode0d.py`, type-2 variant |

Every case is a 4-cell, 1 mm slab, RK3, `integration-variables = cons` (`vnn = 0.25` for
`kinetics_isothermal`, `adiabatic0d` and `darcy/adiabatic0d_type2`; `heated0d` and
`darcy/heated0d_type2` run at `vnn = 0.003125`). With a spatially uniform state and (for
`kinetics_isothermal`) a spatially uniform $T$ held exactly by the BC 302 faces, every cell follows
the same 0-D trajectory, so the spatial discretisation drops out of the comparison.

**`kinetics_isothermal`** (both faces held at $T_0$): with $h$, $h_{pg}$ and $c_p$ all linear and equal
in $T$ (the `--isothermal` tables), $T$ stays at $T_0$ exactly and each reaction follows the closed
form at fixed $T_0$ ($\kappa_r$ evaluated once):

$$
1-\chi_r(t) = \bigl[1+(m_r-1)\kappa_r t\bigr]^{-1/(m_r-1)} \quad (m_r\neq1), \qquad
1-\chi_r(t) = e^{-\kappa_r t} \quad (m_r=1).
$$

**`adiabatic0d`, `heated0d`** (adiabatic or uniformly heated, both closures): the uniform state obeys
the ODE system of `test/pyrolysis/tools/ode0d.py` (classical RK4, pure Python; `adiabatic0d` and its
type-2 twin integrate at step $2\times10^{-4}$ s, `heated0d` and its type-2 twin at $10^{-4}$ s), built
from the case's own tables:

$$
\chi_r' = \kappa_r(T)(1-\chi_r)^{m_r}, \qquad
\pi_{\rm tot} = \sum_r \epsilon_i\rho_i F_r\chi_r', \qquad
E_{1,T}\,\frac{\mathrm dT}{\mathrm dt} = \dot q - \pi_{\rm tot}\bigl(h_{pg}(T)-\bar h\bigr) \;\;(\text{type 1})
$$

with $\bar h = (\rho_v h_v(T) - \rho_c h_c(T))/(\rho_v-\rho_c)$ and $E_{1,T}$ the virgin/char-blended
slope of $\rho h$. The type-2 (Darcy) variant adds the gas-storage term
$(E_{1,T}+\epsilon\rho_g c_{pg})\,\mathrm dT/\mathrm dt = \dot q - \pi_{\rm tot}(h_{pg}-\bar h) + p\,\mathrm d\epsilon/\mathrm dt$,
$\rho_g = pW/(RT)$; FUSS's own type-1 scheme is what the solver integrates, and the type-1/type-2
difference is itself an exact, printed quantity.

## Results

Measured maxima over every field-output time, from each case's `verify.py`:

| Case | max $\lvert\chi-\chi_{\rm exact}\rvert$ | max $\lvert T-T_0\rvert$ [K] | line-integral error (type 1) |
|---|---:|---:|---:|
| `kinetics_isothermal` (T0700/T0900/T1200) | $\le 10^{-8}$ | $\le 10^{-7}$ | $\lvert\dot m_{pg}(s)-s\,\dot m_{pg}(1)\rvert/(s\,\dot m_{pg}(1))\le10^{-12}$ |

| Case | max $\lvert T-T_{\rm oracle}\rvert$ [K] | max $\lvert\chi-\chi_{\rm oracle}\rvert$ | non-uniformity [K] |
|---|---:|---:|---:|
| `adiabatic0d/T1000` (type 1) | 0.0510 at the case's vnn 0.25 | $1.22\times10^{-3}$ | $1.0\times10^{-12}$ |
| `heated0d` (type 1) | 0.0020 at the case's vnn 0.003125 | $6.67\times10^{-7}$ | $0$ |
| `darcy/adiabatic0d_type2` (type 2) | 0.0166 | $9.19\times10^{-4}$ | $1.8\times10^{-7}$ |
| `darcy/heated0d_type2` (type 2) | 0.0011 | $1.04\times10^{-5}$ | $6.7\times10^{-7}$ |

Both quantities refine at second order in $\mathrm dt$ (the RK3 stages freeze the kinetics/energy
split within a stage; halving `vnn` roughly quarters the error): `adiabatic0d/T1000`'s own `verify.py`
documents 0.0510 K at `vnn = 0.25` (the case above), 0.014 K at 0.125, 0.0022 K at 0.05, 0.0005 K at
0.025; `heated0d`'s documents 15.8 K at `vnn = 0.25`, 3.8 K at 0.125, 0.59 K at 0.05, 0.14 K at 0.025,
0.034 K at 0.0125, 0.008 K at 0.00625 (the case above runs at 0.003125, one refinement further).
Reported but not gated: the exact type-2 minus type-1 shift (adiabatic: $0.0696$ K; heated: $-3.914$ K
at the end) is what each type-2 run must reproduce, and (for `heated0d`, when PATO's own TACOT tables
are present) the shift from using PATO's table enthalpy instead of the case's own ($-0.150$ K).

## Pass criteria

- `kinetics_isothermal`: max $\lvert\chi_r-\chi_{\rm exact}\rvert\le10^{-8}$; max $\lvert T-T_0\rvert\le10^{-7}$ K;
  $\rho_s$ recomputed from $\chi$ agrees with the solved $\rho_s$ to $10^{-10}$ relative; the type-1 line
  integral $\dot m_{pg}(s)=s\,\dot m_{pg}(1)$ holds to $10^{-12}$ relative at every cell along the line.
- `adiabatic0d/T1000`: max $\lvert T-T_{\rm or1}\rvert\le0.10$ K, max $\lvert\chi-\chi_{\rm or1}\rvert\le2\times10^{-3}$,
  spatial non-uniformity $\le10^{-6}$ K (type 2: non-uniformity $\le10^{-5}$ K, since the gas storage
  $-\epsilon p$ in $E$ follows the sub-pascal pressure excursion).
- `heated0d`: max $\lvert T-T_{\rm or1}\rvert\le0.02$ K, max $\lvert\chi-\chi_{\rm or1}\rvert\le10^{-5}$,
  non-uniformity $\le10^{-6}$ K (type 2: $\chi$ tolerance $2\times10^{-5}$, non-uniformity $\le10^{-5}$ K).
- The type-2 cases reuse the type-1 `verify.py` (it detects `gas-closure` from `input.ini` and switches
  to the type-2 oracle and tolerances above).
