# Ablating Wall (BC 306/307)

The ablating-wall boundary condition solves a surface energy balance (SEB) with a blowing correction
and, through an equilibrium or finite-rate B′ table, the surface chemistry that turns the wall
temperature and incoming flux into a char/gas mass production. BC 306 is a standalone ablating face
(its own B′ table); BC 307 is the same SEB fed by exchange values from an external coupler (tested
here through FUSS's own `[FUSS-Ablation] coupled-env-file`, a value path used so that FUSS alone can
exercise the coupled form, not a coupler). These cases check the SEB closure itself, its two limiting
and equivalence cases, and the static-mesh wall against PATO.

Setup and run: `test/pyrolysis/ablation/README.md` (`seb_closure`, `zero_blowing`, `ATC21`,
`coupled307`), or the case's own `README.md` (`atc2off_306_static`, `atc2off_306_static_t2`,
`demo_base_pyro`).

---

## Surface energy balance: closure and its fixed point

`ablation/seb_closure`: the ATC 2.1 setup (no-char, blowing parameter $\lambda=0.5$),
with the heating switched off at 1 s instead of 60 s, run to 2 s, so that `OUTPUT/seb.txt` (one row
per 306 face at every solution write) covers three regimes: chemistry on with the $\dot\rho_e u_e C_{H0}$
ramp (0-1 s), PATO's $0<{\rm chemistryOn}<1$ ramp-off window (1-1.1 s), and chemistry off with
re-radiation (1.1-2 s).

**Reference.** The balance closes when the scaled residual $F$ of the Newton solve is at its
tolerance, and when the blowing factor $\beta$ sits at the Picard fixed point PATO uses outside its
own $\beta=1$ exceptions:

$$
\beta = \frac{\zeta}{e^\zeta-1}, \qquad \zeta = \frac{2\lambda(\dot m_{pg}+\dot m_{ch})}{\dot\rho_e u_e C_{H0}},
$$

with $B'_{pg}=\dot m_{pg}/(\dot\rho_e u_e C_{H0}\,\beta)$ and $\dot m_{ch}=B'_{ch}\,\dot\rho_e u_e C_{H0}\,\beta$
consistent to round-off.

**Results.** On every row (all 40): $\lvert F\rvert/\max_k\lvert\text{term}_k\rvert\le10\times$ the
Newton `seb-tol`; on the 18 chemistry-on rows where $\beta\neq1$:
$\lvert\beta-\zeta/(e^\zeta-1)\rvert\le10\times$ `seb-tol`, and $B'_{pg}$, $\dot m_{ch}$ consistent to
$10^{-14}$ relative; every chemistry-off/ramp row has zero convective and gas-advective terms (PATO's
chemistryOn == 0 balance); all three regimes contribute at least one row (18 on, 2 ramp, 20 off).

---

## Zero blowing ($\lambda=0$)

`ablation/zero_blowing`: the same ATC 2.1 setup with $\lambda=0$ on the 306 line, chemistry on
throughout (including the $\dot\rho_e u_e C_{H0}$ ramp over the first 0.1 s). PATO's
`constantLambdaBlowingCorrectionModel` returns $\beta=1$ identically when $\lambda=0$, and FUSS must
match.

**Results.** $\beta=1.0$ bitwise in every `seb.txt` row; the convective term equals
$\dot\rho_e u_e C_{H0}(h_r-h_w)$, recomputed from the same row, to $10^{-15}$ relative; the general
closure identities above hold on every row.

---

## A chemistry-off, zero-optics 306 face is a 303 face

Two cases on `ATC2-chemistryOff`'s mesh and tables (both with the virgin/char emissivity and
absorptivity columns set to 0), face 2 retyped as BC 306 with its env table carrying the same
`hconv = 300`, `Tedge = 2800` row as the BC 303 case: with the chemistry off and zero optics the 306
routine must hand the face to the BC 303 kernel with bit-identical inputs (the env-table
interpolation of two equal rows returns the row value exactly).

| Case | Closure | Compared against | Files identical |
|---|---|---|---|
| `atc2off_306_static` | type 1 | `ATC2-chemistryOff` | 23 of 23 (every `OUTPUT/` file except the wall files) |
| `atc2off_306_static_t2` | type 2 | `darcy/ATC2-chemistryOff-type2` | 24 of 24, `darcy.txt` included; the identical `darcy.txt` and field files show that the 306 face hands the Darcy path the same open-face gas enthalpy and wall pressure as the 303 run |

**Pass criteria.** Every compared file byte-identical (a red control that perturbs the 306 copy's
`hconv` by $10^{-9}$ relative must fail this check).

---

## Ablation Test Case 2.1 against PATO

`ablation/ATC21`: the ablating wall BC 306, TACOT under the type-2 closure, the
equilibrium B′ table and PATO's own `AblationTestCase_2.x` boundary conditions (heating 0.1-60 s,
chemistry off with pure re-radiation to 300 K after 60.1 s), no-char (the char is never consumed, as
PATO's `movingMesh no`).

**Reference.** A local PATO 3.1 run of the same problem with `movingMesh no` and the Ta/p/Xsi
preconditioners set to DIC (the shipped DILU aborts on the symmetric matrices of that setting), at
`maxDeltaT = 0.0005` s (a tenth of the tutorial's own step; at the tutorial's 5 ms step PATO's own
time-discretisation error reaches about 24 K at the wall just after the 60 s switch-off). FUSS samples
PATO's 0.1 s grid: $T$ at the seven probe depths (two bracketing index probes, linear in $x$ then in
$t$); $T_w$ and $\dot m_g$ from `seb.txt` (linear in $t$; a grid time whose two FUSS samples straddle a
breakpoint of the env table — 0.1, 60 or 60.1 s; in these runs the 60 and 60.1 s breakpoints are each
skipped once — interpolating across the 60 s switch-off reads 72 K of sampling error into a wall that
otherwise agrees to 0.1 K just before it); the
virgin-98 %/char-2 % front depths from the field files, by PATO's own rule.

**Results** (averaged/maximised over $t=0.1$-120 s against PATO, and separately against the Amaryllis
organiser reference for the 1.5× check):

| Quantity | max $\lvert\text{FUSS}-\text{PATO}\rvert$ | regression budget (measured + slack) | max $\lvert\text{PATO}-\text{Amaryllis}\rvert$ | $1.5\times$ target met |
|---|---:|---:|---:|---|
| $T$ (7 depths) | 0.880 K | 2.869 K (+2 K) | 25.35 K | yes |
| $T_w$ | 1.391 K | 3.370 K (+2 K) | 98.96 K | yes |
| $\dot m_g$ | $2.30\times10^{-4}$ kg/(m²s) | $2.83\times10^{-3}$ (+10 % of PATO's peak $2.60\times10^{-2}$) | $2.34\times10^{-3}$ | yes |
| virgin-98 % front | $1.15\times10^{-5}$ m | $2.09\times10^{-3}$ m (+10 % of PATO's peak $2.08\times10^{-2}$ m) | $1.13\times10^{-4}$ m | yes |
| char-2 % front | $3.54\times10^{-6}$ m | $2.51\times10^{-4}$ m (+10 % of PATO's peak $2.48\times10^{-3}$ m) | $6.37\times10^{-5}$ m | yes |

Also exactly 0 (no-char): $\dot m_{ch}$. The SEB closure checks of
`seb_closure` above hold on every row.

**Pass criteria.** Each of the five quantities $\le$ its regression budget (the measured maximum,
recorded once, plus 2 K for $T$/$T_w$ or 10 % of PATO's own peak for $\dot m_g$ and the two fronts);
$\dot m_{ch}=0$ exactly; the SEB closure on every row. The $1.5\times$ target
(FUSS within $1.5\times$ PATO's own deviation from Amaryllis) is printed and recorded, not gated.

---

## Coupled ablating wall (BC 307) reproduces BC 306

`ablation/coupled307`: one problem (the base `darcy/ATC2-chemistryOff-type2`, TACOT, type 2, 100
graded cells, shortened to 1 s) solved through two input paths for face 2:

- `run307/`: both face-2 cells typed 307 (the connectivity record, never resolved standalone, plus
  `T_inf` and the B′ file; a coupled face always runs at $\lambda=0$); the exchange values
  ($p_w$, $\dot\rho_e u_e C_H$, $h_r$, $q_{\rm rad}$) come from the test-only
  `[FUSS-Ablation] coupled-env-file = env307.dat`, linear in $t$, clamped;
- `run306/`: face 2 typed 306 with $\lambda=0$ and an env table (`env306.dat`) carrying the same
  values with `chem = 1`, `hconv = 0`, `T_edge = 300`.

Both runs share one `input.ini` (`run307`'s carries one extra line, the `coupled-env-file` key), mesh,
table set, probe list and output cadence (the 306 face consumes no char), and differ only in the
face-2 records
and the env/coupled-env-file values — ATC 2.1's own heating values, $p_w=101325$ Pa,
$\dot\rho_e u_e C_H=0.3$, $h_r=1.5\times10^6$ J/kg, $q_{\rm rad}=0$, held constant to 60 s.

**Results.** Every file of `OUTPUT/` compared (field files, `seb.txt`, wall files,
`bc306-state.dat`): run307 and run306 are byte-identical, 28 of 28 files. A red control (`FUSS_GATE_SELFTEST=1`,
which perturbs `env307.dat`'s $h_r$ from $1.5\times10^6$ to $1.5000001\times10^6$ J/kg, about
$7\times10^{-8}$ relative) is applied by the gate and must fail.

**Multi-block.** `coupled307/run307_4blk`: the same `run307` cut into 4×1 blocks along $x$, the two
307 faces both on block 4; `test/run_mpi_gate.sh` runs it at 4 MPI ranks against 1 and finds every
output byte-identical (27 files). Against the single-block `run307` at one rank: field $T$ relative
$3.2\times10^{-9}$, $p_{pg}$ relative $1.1\times10^{-12}$, $\dot m_{pg}$ relative $1.3\times10^{-9}$,
$\chi$ relative $4.4\times10^{-12}$; the wall's $T_w$ relative $2.8\times10^{-12}$, $\dot m_{pg}$
relative $1.6\times10^{-10}$ (the block cuts give the line-implicit Darcy solve another iterate within
its own tolerance, not a different converged state).

**Pass criteria.** run307 and run306 byte-identical on every output file; the perturbed (red-control)
run307 must differ (FAIL) from run306.

---

## Guard-demo base

`ablation/demo_base_pyro`: a 1-D base case (20 equal cells over 0.05 m, one cell across), face 2 the
ablating wall BC 306 under the type-2 closure with TACOT's tables but the kinetics stripped and a toy
B′ table, chemistry on, run 5 steps. `test/run_guard_demos.sh` copies this base, edits the copy to
violate a specific invariant, and checks that the solver refuses with the expected exit code and
message; the base itself is required only to run its 5 steps to completion unedited (not an accuracy
check against any oracle). It is listed here as the common starting point of those refusal
demonstrations, not as a verification case in its own right.
