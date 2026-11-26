# Ablating Wall

The boundary conditions 306 (ablating wall) and 307 (coupled ablating wall) close the surface of a decomposing solid with a **surface energy balance**: the wall temperature $T_w$ is the root of the balance between the heat conducted into the solid and the heat the boundary layer, the blowing gas and the radiation bring to the surface. This page gives the balance as the code evaluates it, the blowing correction, the B′ tables, what 306 and 307 take as input, and the solution of the balance. The `bc.txt` lines and the table files are described in [Boundary Conditions](../user/boundary-conditions.md#ablating-wall-306-and-coupled-ablating-wall-307), the keys in [Input File](../user/input.md#ablating-wall-306-307).

The wall **blows its pyrolysis gas and consumes no char**: $B'_c = \dot m_c = 0$, while the wall enthalpy $h_w$ is still read from the B′ table.

---

## Surface Energy Balance

With the chemistry on, the wall temperature is the root of

$$
F(T_w) = q_\text{cond} - \rho_e u_e C_{H0}\,\beta\,(h_r - h_w) - \dot m_{pg}\,\bigl(h_g - h_w\bigr) - \alpha_s\,q_\text{rad} - \varepsilon_s\,\sigma\,T_\infty^4 + \varepsilon_s\,\sigma\,T_w^4 = 0 ,
$$

that is, the heat conducted into the solid is

$$
q_\text{cond} = \rho_e u_e C_{H0}\,\beta\,(h_r - h_w) + \dot m_{pg}\,(h_g - h_w) + \alpha_s\,q_\text{rad} + \varepsilon_s\,\sigma\,\bigl(T_\infty^4 - T_w^4\bigr).
$$

| Term | Expression | Evaluated at |
|------|------------|--------------|
| conduction into the solid | $q_\text{cond} = G\,(T_w - T_c) + C$ | the wall kernel of [Walls](numerics.md#walls-301305), with $k_s$ of the boundary cell ($k(T_c)$ of `properties.dat` without `[FUSS-Pyrolysis]`) |
| convection, blowing-corrected | $\rho_e u_e C_{H0}\,\beta\,(h_r - h_w)$ | $h_w$ from the B′ table at $(p_w, B'_{pg}, T_w)$ |
| pyrolysis-gas advection | $\dot m_{pg}\,(h_g - h_w)$ | $h_g$ from the gas table at $(p, T_w)$ |
| char advection | $\dot m_c\,(h_c - h_w)$ | zero: $\dot m_c = 0$ |
| absorbed radiation | $\alpha_s\,q_\text{rad}$ | $\alpha_s$ at $T_w$ |
| background radiation | $\varepsilon_s\,\sigma\,T_\infty^4$ | $\varepsilon_s$ at $T_w$ |
| emitted radiation | $-\varepsilon_s\,\sigma\,T_w^4$ | |

- $T_c$ is the temperature of the boundary cell. $G$ and $C$ are the half-cell conductance and the cross-diffusion term of the wall faces.
- $\varepsilon_s$ and $\alpha_s$ blend the emissivity and absorptivity columns of `virgin.dat` and `char.dat` at $T_w$ with the virgin share of the cell (see [Effective Properties](pyrolysis.md#effective-properties)).
- $\sigma = 5.6704085581\times10^{-8}$ W/(m²·K⁴) here, the value of PATO. BC 304 and 305 use $5.67\times10^{-8}$.
- The gas enthalpy is taken at $p = p_w$ under the type-2 closure and at `p-gas` under type 1.
- $\rho_e u_e C_{H0}$ is the unblown heat-transfer coefficient; with the chemistry on it is floored at $10^{-6}$ kg/(m²·s).
- The cell receives $A_f\,q_\text{cond}$ at the converged $T_w$; the wall output writes it as $q_w$.

**Blowing mass flux.** $\dot m_{pg}$ (outward positive) comes from the gas closure:

| Closure | $\dot m_{pg}$ |
|---------|---------------|
| type 2 | the outward mass flux of the face from this step's pressure solve |
| type 1 | the line integral of the boundary cell, from the previous evaluation of the stage loop |
| no `[FUSS-Pyrolysis]` | 0 |

In the interior the gas leaves through the same face with $h_g$ at the same $T_w$ ([type 1](pyrolysis.md#type-1-gas-closure); [type 2](numerics.md#type-2-gas-enthalpy-flux)). The term $\dot m_{pg}(h_g - h_w)$ therefore changes the enthalpy of the blown gas from its value at the wall to the wall enthalpy of the table.

### Chemistry off

The env column `chem` switches the chemistry: it is on when $\text{chem} \ge 1 - 10^{-12}$. Otherwise (0, or a ramp $0 < \text{chem} < 1$ between two rows) the balance is

$$
F(T_w) = q_\text{cond} - h_\text{conv}\,(T_\text{edge} - T_w) - \alpha_s\,q_\text{rad} - \varepsilon_s\,\sigma\,T_\infty^4 + \varepsilon_s\,\sigma\,T_w^4 = 0 ,
$$

without the gas advection, the B′ table and the blowing correction. With $\text{chem} = 0$ and surface optics that are zero (no `[FUSS-Pyrolysis]`, or emissivity columns that are zero and absorptivity that is zero or meets $q_\text{rad} = 0$), the face solves the balance of BC 303 with $h_\text{conv}$, $T_\text{edge}$ and $q_\text{rad} = 0$, with the same arithmetic. A 306 face without `[FUSS-Pyrolysis]` must have $\text{chem} = 0$ in every row.

---

## Blowing Correction

The blowing gas thickens the boundary layer and reduces the heat-transfer coefficient by the factor

$$
\beta = \frac{\ln(1 + 2\lambda B')}{2\lambda B'},
\qquad
B' = B'_{pg} + B'_c = \frac{\dot m_{pg}}{\rho_e u_e C_{H0}\,\beta} \quad (B'_c = 0) ,
$$

with $\beta = 1$ when $\lambda = 0$, $B' = 0$, $\dot m_{pg} \le 10^{-6}$ kg/(m²·s) or $\rho_e u_e C_{H0} \le 10^{-4}$ kg/(m²·s).

- $\beta$ appears on both sides. At each Newton iterate of $T_w$ it is brought to its fixed point by Picard iteration, $\beta_{k+1} = \beta(B'(\beta_k))$, until $\lvert\beta_{k+1} - \beta_k\rvert \le$ `seb-tol`. An iteration that does not settle within `seb-max-iter` stops the run. The formula and its exceptions are those of PATO [2]; PATO updates $\beta$ once per call, FUSS iterates it to its fixed point.
- At the fixed point, with $\zeta = 2\lambda\,\dot m_{pg}/(\rho_e u_e C_{H0})$, the correction is $\beta = \zeta/(e^\zeta - 1)$: setting $x = 2\lambda B'$ gives $\zeta = x\beta = \ln(1 + x)$.
- $\lambda$ is the `lambda` of the 306 line ($\lambda = 0$ switches the correction off). A 307 face runs at $\lambda = 0$: the coefficient of the coupled solver already includes the blowing.
- The Picard iteration starts from the last converged $\beta$, which the restart file keeps.

---

## B′ Tables

The B′ table gives the char blowing rate $B'_c$ and the wall enthalpy $h_w$ of the surface chemistry against the wall pressure, the gas blowing rate and the wall temperature. Its rows are `p B'g B'c T h_w`; a four-dimensional table adds $\rho_e u_e C_M$ as a first column.

**Structure.** The file is a sequence of blocks, checked when it is read:

- p-blocks: runs of equal $p$, strictly decreasing from block to block;
- in each p-block, B′g sub-blocks: runs of equal $B'_g$, strictly decreasing;
- in each sub-block at least two rows, $T$ monotonic (an increasing sub-block may repeat a $T$; rows are stored increasing);
- a four-dimensional file stacks such tables, one per $\rho_e u_e C_M$ level, levels increasing.

Each sub-block has its own $T$ rows, so the table need not be rectangular.

**Interpolation.** A query $(p, B'_g, T)$ is interpolated in nested steps:

| Axis | Rule |
|------|------|
| $\rho_e u_e C_M$ (four-dimensional table) | linear in $\ln(\rho_e u_e C_M)$ between the two bracketing levels |
| $p$ | linear in $\ln p$ between the two bracketing p-blocks |
| $B'_g$ | in each p-block, linear between the two bracketing sub-blocks, within that block's own range |
| $T$ | in each sub-block, linear between the row below and the first row at or above $T$ |

A query on a grid value returns the table value. The balance queries $B'_g = \max(B'_{pg}, 0)$, and the fourth axis with $\rho_e u_e C_{H0}\,\beta$.

**Out of range.** A query outside the table is clamped on each axis to the nearest table value, and each face counts, per axis and once per stage, whether the query of its converged balance was clamped (written to `seb.txt`; the queries of the slope evaluations and of the bracketing iterates are not counted, nor those of the discarded first pass of a new type-2 run). With `bprime-out-of-range = abort` a clamp of $p$, $B'_g$ or $\rho_e u_e C_M$ in that query stops the run (exit code 1); a clamp of $T$ is counted only.

The balance uses $h_w$; the $B'_c$ of the table is not applied.

---

## Inputs of 306 and 307

The two faces solve the same balance. They differ in where its boundary-layer values come from.

| Quantity | 306 | 307 |
|----------|-----|-----|
| $p_w$ | env table | exchange with the coupled solver |
| $\rho_e u_e C_{H0}$ | env table | exchange |
| $h_r$ | env table | exchange |
| $q_\text{rad}$ | env table | exchange |
| `chem` | env table | 1 |
| $h_\text{conv}$, $T_\text{edge}$ | env table | 0, $T_\infty$ (unused with the chemistry on) |
| $\lambda$ | `bc.txt` | 0 |
| $T_\infty$ | `bc.txt` | `bc.txt` |
| B′ table | `bc.txt`, or `none` when every env row has `chem = 0` | `bc.txt`, required |

- **306.** The env table holds the rows `t p_w rhoeUeCH0 h_r q_rad chem hconv T_edge`. It is interpolated linearly in $t$ and held at its first and last rows outside its range.
- **307.** The coupled solver writes $(p_w, \rho_e u_e C_H, h_r, q_\text{rad})$ before every FUSS step. A step without them ($p_w = 0$) stops the run; at set-up the face keeps `p-wall`. A standalone run can read them from the test key `coupled-env-file`. A 307 face needs `[FUSS-Pyrolysis]`.
- Both are evaluated **once per step**, at the end time $t^{n+1}$, before the pressure solve of the type-2 closure, which takes the face's $p_w$ as its Dirichlet value. Every Runge–Kutta stage of the step uses the same values.
- Under the type-2 closure the face is open to the gas at $p_w$ whatever `open-bc-types` holds, and $p_w$ must lie inside the gas table: the set-up checks the env rows of a 306 face; a 307 exchange value outside it stops the run at the step. Under type 1 the face must be the exit face of `outflow-direction`.
- Neither face runs with multigrid.

---

## Solution of the Balance

The balance is solved at every Runge–Kutta stage, for each face, by a safeguarded Newton iteration on $T_w$.

1. **Bracket.** $[\max(T_c - 500,\ T_\text{lo} + 1),\ T_\text{hi} - 1]$ K, with $[T_\text{lo}, T_\text{hi}]$ the inversion bracket of the cell's material ([Stored Energy](pyrolysis.md#stored-energy)); without `[FUSS-Pyrolysis]`, $[\max(T_c - 500, 1),\ T_c + 10^4]$ K.
2. **Start.** The last converged $T_w$ of the face (the cell temperature at the first solve), moved into the bracket.
3. **Iterate.** With the chemistry on, $\beta$ is first brought to its fixed point at the current $T_w$. The slope is the central difference $F'(T_w) = [F(T_w + 0.5) - F(T_w - 0.5)]/(1\ \text{K})$ at fixed $\beta$.
4. **Stop.** At $F = 0$, or when $F' > 0$ and both $\lvert F\rvert \le F'\cdot$`seb-tol`$\cdot T_w$ (the Newton step is at most `seb-tol` $T_w$) and $\lvert F\rvert \le$ `seb-tol` $\cdot \max_k \lvert t_k\rvert$, the largest of the seven terms. The current iterate is returned.
5. **Update.** $F$ increases with $T_w$, so the sign of $F$ moves one end of the bracket to $T_w$. The Newton step $T_w - F/F'$ is taken when it falls inside the bracket; otherwise, or when $F' \le 0$, the iterate bisects the bracket.

`seb-max-iter` iterations (default 50) without convergence stop the run (exit code 1) and print the terms of the last iterate.

- The converged $T_w$, $\beta$ and the terms are kept per face: `OUTPUT/<prefix>seb.txt` writes the terms, and the restart file keeps $T_w$ and $\beta$, the starting values of the next solve.
- The ghost cells of the face hold $2\,T_w - T_g$ with the last converged $T_w$ (see [Walls](numerics.md#walls-301305)).
- The balance adds no limit to the time step.

---

## References

1. M. Rotondi, *Modeling of Conjugate Heat Transfer Including Ablation and Pyrolysis for Propulsive Applications*, PhD thesis, Sapienza University of Rome, 2023, ch. 3.
2. J. Lachaud, N. N. Mansour, "Porous-material analysis toolbox based on OpenFOAM and applications," *J. Thermophys. Heat Transfer* 28 (2), 2014.
