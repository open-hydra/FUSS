# Energy and Mass Audit, and Pure Conduction on the Pyrolysis Path

Two checks that do not need an external oracle: the discrete energy and solid-mass identities of the
pyrolysis update, closed to round-off from the solver's own accumulators and cross-checked by an
independent recomputation from the field file (`energy_audit`); and pure conduction on the same
mesh, material tables and boundary setup as [Ablation Test Case 1](pyrolysis-atc1.md), with the
reactions stripped, against the exact step-heated-slab series (`conduction/graded100`).

Setup and run: `test/pyrolysis/energy_audit/README.md`, `test/pyrolysis/conduction/README.md`.

---

## Energy and mass audit

### Problem

A 40-cell TACOT slab (1 mm cells), face 2 ramped 300 → 1644 K in 1 s through a BC 302 time file, face 1
adiabatic (BC 301), 10 s, RK3/`cons`/`vnn = 0.25`, type-1 gas closure. Four probes (`p05mm`, `p15mm`,
`p35mm`, `p75mm`) sit at index positions 40, 39, 37, 33 (i.e. 0.5, 1.5, 3.5, 7.5 mm below the heated
face).

At every solution write the solver appends a row of
`OUTPUT/energy-budget.txt`: `t, E_tot, M_s, Q_out, Q_gas_out, M_gas_out, Q_src`, where `E_tot = sum V
E(T, chi)`, `M_s = sum V rho_s`, `Q_out` is the net outgoing flux (conduction + BC + gas advection)
integrated with the scheme's own Runge-Kutta weights, `Q_gas_out`/`M_gas_out` the energy and mass
carried out through the outflow face, and `Q_src` the integrated volumetric source. Because the
accumulators advance with the same RK weights as the state, the balances below are identities of the
discrete update, not approximations.

### Reference

None external: checks (i)-(iii) are identities of FUSS's own update (the balance must close to
round-off and the field file must reproduce the budget's own state sum); check (iv) is a `vnn`
refinement of FUSS against itself, judged by its observed order (the scheme is first order in $\mathrm dt$
here), not by a stored value.

### Results

| Check | Quantity | Result |
|---|---|---|
| (i) energy balance | $\lvert E(t)-E(0)+Q_{\rm out}-Q_{\rm src}\rvert / \max(\lvert E(0)\rvert,\lvert Q_{\rm out}\rvert)$ | $\le 10^{-10}$ |
| (ii) solid-mass balance | $\lvert M_s(t)-M_s(0)+M_{\rm gas,out}\rvert / M_s(0)$ | $\le 10^{-10}$ |
| (iii) independent state sum | $E_{\rm tot}$, $M_s$ recomputed from `field.tec` (T, $\chi_*$), the input tables and the node-coordinate cell volumes, vs. the solver's own columns | $\le 10^{-10}$ relative |
| (iv) time-discretisation order | the three-`vnn` (0.25/0.125/0.0625) probe-series refinement, interpolated onto $t=2.0,2.1,\dots,10.0$ s: ratio $r=\max\lvert T_{0.25}-T_{0.125}\rvert/\max\lvert T_{0.125}-T_{0.0625}\rvert$ and the observed order $\log_2 r$ | $0.9484$ K / $0.4234$ K, $r=2.240$, order $1.163$ ($\max\lvert T_{0.25}-T_{0.0625}\rvert=1.000$ K, printed) |
| (v) zero-kinetics self-consistency | a `--nokin` rerun: every $\chi_r$ and $\dot m_{pg}$ exactly 0, $\rho_s=\rho_v$ to $10^{-14}$ relative, $Q_{\rm gas,out}=M_{\rm gas,out}=0$ exactly | exact |
| (vi) cross-thread determinism | `field.tec` and `energy-budget.txt`, 1 vs. 8 OpenMP threads | byte-identical |

The 1 s ramp (rather than a shorter one) was chosen because the staircase sampling of a 0.1 s ramp
aliased under mesh/`vnn` refinement and the probe differences *grew* with refinement (3.3, 4.9, 1.9 K);
with the 1 s ramp the coarsest `vnn` resolves it (coarsest $\mathrm dt\approx0.17$ s) and the ratio test
of check (iv) measures an order.

### Pass criteria

(i), (ii), (v): the stated tolerances above, all at $10^{-10}$ or exact. (iii): $10^{-10}$ relative.
(iv): the observed order $\log_2 r$ in $[0.8, 1.2]$ (expected $p=1$, the band $0.8p$-$1.2p$; a sub-check
whose companion directories are not given is skipped, not failed). (vi): byte-identical files.

---

## Pure conduction, graded ATC1 mesh

### Problem

The geometry and graded mesh of ATC1 (100 cells over 0.05 m, graded 0.1 toward the heated face; the check also
runs 50 and 200 cells with the same grading), with
constant-property tables ($c_p = 1000$ J/(kg·K), $k = 0.4$ W/(m·K), $\rho_v = 280$ kg/m³, so
$\alpha = k/(\rho c_p) = 1.4286\times10^{-6}$ m²/s) and every pyrolysis reaction stripped
(`nPyroReac = 0` throughout). Face 2 is stepped to 1644 K at $t=0$ through a BC 302 time file
(`0 1644 / 60 1644`); face 1 is adiabatic (BC 301); $T(0) = 298$ K; 60 s. Ten index-position probes
sit in pairs near the FIAT thermocouple depths of ATC1 (1, 2, 4, 8, 16 mm), plus a back-face probe.

### Reference

The exact solution of the step-heated, back-insulated slab with constant properties
(`test/pyrolysis/tools/slab_series.py`), a cosine series in the depth $d = L - x$ measured from the
insulated face:

$$
T(d, t) = T_w + (T_0-T_w)\sum_{n\ge0}\frac{2(-1)^n}{\lambda_n}\cos\!\Bigl(\frac{\lambda_n(L-d)}{L}\Bigr)
\exp\!\Bigl(-\frac{\lambda_n^2\alpha t}{L^2}\Bigr), \qquad \lambda_n=\frac{(2n+1)\pi}{2},
$$

evaluated with 4000 terms (exact to round-off for $\alpha t/L^2\ge10^{-7}$; checked against the
semi-infinite erf solution to $3\times10^{-12}$ K at $t=0.5$ s, 1 mm). `verify.py` evaluates the series at
the two cell centres that bracket each depth (1, 2, 4, 8, 16 mm) and at the probe sample nearest each
$t\in\{0.5,1,2,4,8,16,30\}$ s, and interpolates the error $T-T_{\rm exact}$ of those two cells linearly to
the depth. The error, not $T$, is interpolated: the exact solution exists where FUSS has values, and an
interpolated $T$ would add the exact profile's own interpolation error, which changes with where the
centres fall on each mesh.

### Results

| Quantity | 50 cells | 100 cells | 200 cells | order (100/200) |
|---|---:|---:|---:|---:|
| (a) max error over every depth and time | 4.3727 K | 1.1331 K | 0.2888 K | 1.972 |
| (b) error at 1 mm, $t=1$ s | 2.0276 K | 0.5211 K | 0.1325 K | 1.975 |
| (c) max error at $\ge2$ mm, every time | 2.0301 K | 0.5251 K | 0.1329 K | 1.983 |

The largest error sits at 1 mm and $t=0.5$ s on every mesh.

### Pass criteria

Each of (a), (b), (c) falls from 50 to 100 to 200 cells and its order on the finest pair lies in
$[1.6, 2.4]$ (expected $p=2$, second order in space, the band $0.8p$-$1.2p$). No stored FUSS value enters.
