# Literature Benchmarks of the Pyrolysis and Ablation Models

Two published decomposition problems with closed-form or digitized-reference solutions, used as
independent checks outside FUSS's own TACOT-based cases: a charring front with Darcy gas flow solved
in closed form, and a decomposing slab solved by three independent finite-element codes.

Setup and run: `test/pyrolysis/lit/baum_atreya2007/README.md`, `test/pyrolysis/lit/zhang2010/README.md`.

---

## Baum & Atreya (2007): charring front with Darcy flow

### Problem

A semi-infinite virgin solid at $T_0=298.16$ K, its surface held at $T_s=1298.16$ K from $t=0$,
decomposing instantly at $T_p=598.16$ K from $\rho_v=700$ to $\rho_c=140$ kg/m³, with conductivities
$k_v=0.209$, $k_c=0.071$ W/(m·K) and specific heats $c_v=1400$, $c_c=1100$ J/(kg·K) (the Di Blasi
property set). The gas flows through the char to the surface at $p_s=101325$ Pa (Darcy, no
storage in the quasi-steady mode; with storage, $\epsilon_c=0.76$, in the full mode). FUSS: type-2 gas
closure, $K_c=10^{-15}$ m² (the virgin is effectively impermeable, $K_v=10^{-20}$ m²), a 0.04 m slab,
BC 302 at the surface, BC 301 (sealed) at the back. Two gas-heat-capacity choices: $c_g=1475$
J/(kg·K) (`--cg 1475`, which makes the decomposition heat $Q=497302.6$ J/kg independent of $T$ and
puts the gas-enthalpy advection through the char into the exact solution; the full-storage cases
`full_N100`/`full_N200` run at this value) and $c_g=0$ (the paper's own case; `qs_N100` is the
quasi-steady, no-storage mode, the paper's $\beta=0$, and runs at $c_g=0$).

### Reference

H.R. Baum, A. Atreya, *Proceedings of the Combustion Institute* **31** (2007) 2633. The front
$X=C^\star\sqrt{\alpha_c t}$, with $C^\star$ from the **corrected** front condition: the paper's printed
eq. 17 has $\lambda_c/\lambda_v$ where $\rho_c c_c/(\rho_v c_v)$ belongs (= Atreya & Baum, 2002, eq. 27);
FUSS computes $C^\star=2\lambda$ by bisection on the corrected Rankine-Hugoniot front condition
(`test/pyrolysis/tools/darcy_exact.py`, problem 3, shared with the [sharp-front exact
benchmark](darcy-exact.md#decomposition-front-with-transpiration-through-the-char)). With the
corrected condition and $c_g=1475$: $\lambda=0.19428991$ ($C^\star=0.388580$); with $c_g=0$ (the
paper's case): $\lambda=0.22215026$ ($C^\star=0.444301$).

### Results

Averages over $30\le t\le300$ s:

| Quantity | N100 | N200 |
|---|---:|---:|
| front error (mass-based) | 15.01 % | 7.36 % |
| $L_2$ temperature error (char + preheat, $x\le3X$) | 27.48 K | 15.09 K |
| gas-mass ledger $\lvert\text{out}+\text{stored}-\text{stored}(0)-\text{lost}\rvert/\text{lost}$ | $-1.10\times10^{-14}$ | $+1.05\times10^{-14}$ |
| $\langle p^2-p_s^2\rangle$ at 1.1 mm, relative to exact | — | $-6.63\,\%$ |

Orders on N100/N200: front error 1.028, $L_2$ temperature error 0.864. The quasi-steady run (`qs_N100`,
$c_g=0$, no storage) gives a gas-mass ledger of $-8.49\times10^{-5}$: the quasi-steady closure carries the
gas released in one step beyond its lagged prediction, so this ledger is not an identity (unlike the
full run's, above); its bound is the run's own one-step share of the solid lost, $1.515\times10^{-3}$
(3432 steps over 31 budget rows). The char pressure at the 1 mm cell is printed, not gated: that cell
sits at 1.0 mm on N100 and 1.1 mm on N200, so the two values are not errors at one point.

### Pass criteria

- Front error and $L_2$ temperature error each fall from N100 to N200 with an observed order in
  $[0.8, 1.2]$ (expected $p=1$: a cell-by-cell front is first order by design, and the temperature error
  follows the front; the band is $0.8p$-$1.2p$). No stored FUSS value enters.
- Every full run's gas-mass ledger $\le10^{-10}$ (round-off: the discrete update's own identity).
- The quasi-steady run's ledger $\le$ its one-step share: (the largest solid loss between two budget
  rows) / (the steps per row, taken uniform) / (the total solid lost), from the run's own budget file and
  step count (one step of carried gas, not an identity).

---

## Zhang (2010): three-code decomposing slab

### Problem

A 12.7 mm slab, $\rho_0=1700\to\rho_f=1255$ kg/m³, $c_p=1100$, $c_g=3000$ J/(kg·K), $k=0.30$ W/(m·K),
single Arrhenius reaction $A=5\times10^{28}$ s⁻¹, $E=3.62\times10^5$ J/mol, $n=4.6$, decomposition
heat $Q=-2\times10^5$ J/kg (Henderson-form model); exposed face
$q_{\rm rad}=25$ kW/m² $-\sigma\epsilon T_s^4$ ($\epsilon=1$, no back-radiation term) $+\,h(T_\infty-T_s)$,
$h=10$ W/(m²K), $T_\infty=300$ K; unexposed face adiabatic; $T_{\rm ini}=300$ K; 2000 s. FUSS: type 1
(the gas leaves through the exposed face), virgin and char sharing $h(T)=Q+c_p(T-T_{\rm ref})$,
$T_{\rm ref}=273$ K (so $h_g-\bar h$ reproduces the thesis's $h_g-h$ term for term); the exposed face is
one BC 305 line whose reference temperature reproduces the three terms of the energy balance for
every $T$. 50 cells (even, so the mid-plane falls between two cell centres).

### Reference

Zhenyu Zhang, *Thermo-Mechanical Behavior of Polymer Composites Exposed to Fire*, PhD dissertation,
Virginia Tech, 29 June 2010, Section 2.4.1 and Table 2.1, Figs. 2.2-2.3: the problem solved by three
independent codes (ABAQUS UMATHT with 20 3-D solid elements, a 1-D UEL element, and the MATLAB FEM
code of Summers). The comparison digitizes the vector paths of Figs. 2.2 (temperature histories) and
2.3 (exposed-face mass-loss rate) from the thesis PDF, calibrated on the axis tick marks (time and
value residuals at most about one device quantum, $\approx0.6$-$0.7$ s and $0.16$-$1.66\times10^{-6}$
in value; the uncertainty model and provenance are recorded beside the digitized data). MATLAB is the
comparison oracle; UEL and UMATHT give the spread between independent codes.

### Results

For $t>20$ s:

| Surface | max $\lvert T_{\rm FUSS}-T_{\rm MATLAB}\rvert$ | mean | the three codes' own spread |
|---|---:|---:|---:|
| exposed | 3.29 K | $-2.81$ K | 4.75 K |
| mid-plane | 2.51 K | $-0.67$ K | 5.26 K |
| unexposed | 2.08 K | $+0.28$ K | 5.10 K |

Peak exposed-face gas mass flux: FUSS $4.369\times10^{-3}$ kg/(m²s) at 904 s, against the three codes'
own peaks of MATLAB $4.598\times10^{-3}$ at 890 s, UEL $4.489\times10^{-3}$ at 910 s, UMATHT
$4.582\times10^{-3}$ at 900 s (MATLAB's peak is a plateau at $4.6\times10^{-3}$ from 880.5 to
909.7 s).

### Pass criteria

- For $t>20$ s, max $\lvert T_{\rm FUSS}-T_{\rm MATLAB}\rvert$ at each plane $\le$ that plane's spread of the
  three codes, max over $t>20$ s of $\lvert T_{\rm UEL}-T_{\rm MATLAB}\rvert$ and
  $\lvert T_{\rm UMATHT}-T_{\rm MATLAB}\rvert$, computed from the digitized curves each time the check runs
  (4.75, 5.26, 5.10 K).
- The FUSS peak gas mass flux lies within $[0.95,1.05]\times$ the lowest/highest of the three codes'
  own peaks, and its time within 30 s of their mean peak time.
