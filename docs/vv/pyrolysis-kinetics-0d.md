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
the ODE system of `test/pyrolysis/tools/ode0d.py`, built from the case's own tables:

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

The right-hand side jumps at every integer kelvin (the slopes of the 1 K tables) and at each reaction
threshold, so the oracle integrates with classical RK4 steps (at most $2\times10^{-4}$ s for `adiabatic0d`
and its twin, $10^{-4}$ s for `heated0d` and its twin) that each stay on one table segment and one on/off
state of every reaction: a step that would leave either is shortened to end on it. Halving the step
changes the oracle by at most $5.1\times10^{-11}$ K (adiabatic) and $1.7\times10^{-9}$ K (heated), five or
more orders below FUSS's own error.

## Results

Measured maxima over every field-output time, from each case's `verify.py`:

| Case | max $\lvert\chi-\chi_{\rm exact}\rvert$ | max $\lvert T-T_0\rvert$ [K] | line-integral error (type 1) |
|---|---:|---:|---:|
| `kinetics_isothermal` (T0700/T0900/T1200) | $\le 10^{-8}$ | $\le 10^{-7}$ | $\lvert\dot m_{pg}(s)-s\,\dot m_{pg}(1)\rvert/(s\,\dot m_{pg}(1))\le10^{-12}$ |

Each type-1 case runs on the ladder `vnn`, `vnn`/2, `vnn`/4 (0.25, 0.125, 0.0625 for `adiabatic0d/T1000`;
0.003125, 0.0015625, 0.00078125 for `heated0d`), against the oracle at every field-output time:

| Case | max $\lvert T-T_{\rm or1}\rvert$ [K] | order | max $\lvert\chi-\chi_{\rm or1}\rvert$ | order |
|---|---|---:|---|---:|
| `adiabatic0d/T1000` | $5.0981\times10^{-2}$, $1.3719\times10^{-2}$, $3.4745\times10^{-3}$ | 1.981 | $1.216\times10^{-3}$, $3.427\times10^{-4}$, $8.969\times10^{-5}$ | 1.934 |
| `heated0d` | $2.1361\times10^{-3}$, $5.2196\times10^{-4}$, $1.3087\times10^{-4}$ | 1.996 | $5.838\times10^{-7}$, $1.321\times10^{-7}$, $3.320\times10^{-8}$ | 1.993 |

(order: on the finest pair.) Both refine at second order in $\mathrm dt$, the RK3 stages freezing the
kinetics/energy split within a stage. The spatial non-uniformity stays at round-off: at most
$2.4\times10^{-12}$ K (`adiabatic0d`) and exactly 0 (`heated0d`).

The type-2 twins run against the type-2 oracle, `heated0d_type2` on the same ladder `vnn`, `vnn`/2, `vnn`/4
and `adiabatic0d_type2` on `vnn`/4, `vnn`/8, `vnn`/16 (0.0625, 0.03125, 0.015625):

| Case | max $\lvert T-T_{\rm or2}\rvert$ [K] | order | max $\lvert\chi-\chi_{\rm or2}\rvert$ | order |
|---|---|---:|---|---:|
| `darcy/adiabatic0d_type2` | $8.5622\times10^{-3}$, $2.4896\times10^{-3}$, $6.1842\times10^{-4}$ | 2.009 | $4.351\times10^{-5}$, $1.144\times10^{-5}$, $2.914\times10^{-6}$ | 1.973 |
| `darcy/heated0d_type2` | $7.4875\times10^{-4}$, $1.9209\times10^{-4}$, $3.9729\times10^{-5}$ | 2.273 | $9.281\times10^{-6}$, $2.386\times10^{-6}$, $5.750\times10^{-7}$ | 2.053 |

The adiabatic ladder starts at `vnn`/4 because on the first two halvings from the case's `vnn` = 0.25 the
temperature error does not yet fall ($1.6575\times10^{-2}$, $1.9851\times10^{-2}$, $8.5622\times10^{-3}$ K at
0.25, 0.125, 0.0625): those steps are outside the asymptotic range. The spatial non-uniformity stays below
$1.9\times10^{-7}$ K (adiabatic) and $6.8\times10^{-7}$ K (heated). With $K=10^{-7}$ m² the pressure excursion
$\max p_{pg}-p_w$ stays at most 0.075 Pa (adiabatic) and 0.503 Pa (heated), so these cases test the
kinetics/energy splitting under type 2, not the time order of the pressure step.

Reported but not gated: the exact type-2 minus type-1 shift (adiabatic: $0.0696$ K; heated: $-3.914$ K
at the end) is what each type-2 run must reproduce, and (for `heated0d`, when PATO's own TACOT tables
are present) the shift from using PATO's table enthalpy instead of the case's own ($-0.150$ K).

## Pass criteria

- `kinetics_isothermal`: max $\lvert\chi_r-\chi_{\rm exact}\rvert\le10^{-8}$; max $\lvert T-T_0\rvert\le10^{-7}$ K;
  $\rho_s$ recomputed from $\chi$ agrees with the solved $\rho_s$ to $10^{-10}$ relative; the type-1 line
  integral $\dot m_{pg}(s)=s\,\dot m_{pg}(1)$ holds to $10^{-12}$ relative at every cell along the line.
- `adiabatic0d/T1000`, `heated0d` (type 1), on the ladder `vnn`, `vnn`/2, `vnn`/4: max
  $\lvert T-T_{\rm or1}\rvert$ and max $\lvert\chi-\chi_{\rm or1}\rvert$ each fall on every rung and their observed
  order on the finest pair lies in $[1.6, 2.4]$ (expected $p=2$, the band $0.8p$-$1.2p$); no stored FUSS
  value enters. The spatial non-uniformity of each run stays below its own round-off bound,
  $N_{\rm steps}\times3\times\varepsilon\times\max T$ for $T$ and $N_{\rm steps}\times3\times\varepsilon$ for $\chi$
  ($\varepsilon=2^{-52}$; one rounding per RK3 stage per step).
- The type-2 cases reuse the type-1 `verify.py` (it detects `gas-closure` from `input.ini` and switches
  to the type-2 oracle) with the same rule, $p=2$: both errors fall on every rung and their order on the
  finest pair lies in $[1.6, 2.4]$, `heated0d_type2` on `vnn`, `vnn`/2, `vnn`/4 and `adiabatic0d_type2` on
  `vnn`/4, `vnn`/8, `vnn`/16. The non-uniformity bound is $10^{-5}$ K, a floor set from measured draws
  rather than a physics budget (the gas storage $-\epsilon p$ in $E$ makes $T$ follow the sub-pascal
  pressure excursion).
