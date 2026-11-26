# Pyrolysis

With a `[FUSS-Pyrolysis]` section FUSS solves a porous solid that decomposes into a char and a pyrolysis gas. This page gives the model as the code integrates it: the solid state and its reactions, the stored energy, the effective properties, and the two closures of the gas flow. The keys and files are listed in [Input File](../user/input.md#pyrolysis); the ablating wall that closes the surface is described in [Ablating Wall](ablation.md).

---

## Solid State

### Phases, sub-phases and reactions

A material consists of $N_p$ solid phases (`nSolidPhases`). Phase $i$ has the initial volume fraction $\varepsilon_{i,0}$ (`epsI[i]`) and the intrinsic density $\rho_{i,0}$ (`rhoI[i]`), and decomposes through `nPyroReac[i]` reactions. Reaction $j$ of phase $i$ consumes the sub-phase of mass fraction $F_{ij}$ (`F[i][j]`) of that phase. Below, $r = (i, j)$ runs over all reactions of the material.

The state of reaction $r$ is its **advancement** $\chi_r \in [0, 1]$: 0 virgin, 1 fully decomposed. The constants are fixed per material:

$$
\rho_v = \sum_i \varepsilon_{i,0}\,\rho_{i,0},
\qquad
\rho_c = \rho_v - \sum_r \varepsilon_{i,0}\,\rho_{i,0}\,F_r ,
$$

the densities of the virgin and of the fully charred material. The solid density of a cell and its **virgin fraction** $\tau$ follow from the advancements:

$$
\rho_s = \rho_v - \sum_r \varepsilon_{i,0}\,\rho_{i,0}\,F_r\,\chi_r ,
\qquad
\tau = \sum_r w_r\,(1 - \chi_r),
\quad
w_r = \frac{\varepsilon_{i,0}\,\rho_{i,0}\,F_r}{\rho_v - \rho_c} .
$$

- $\tau = 1$ for the virgin material and $\tau = 0$ for the char; the code clamps $\tau$ to $[0, 1]$.
- $\tau \equiv 1$ for a material without reaction, or whose reactions decompose no mass ($\rho_v = \rho_c$), where the weights $w_r$ are undefined.
- By these definitions $\tau\,\rho_v + (1 - \tau)\,\rho_c = \rho_s$: the solid mass splits into a virgin and a char share.

The gas volume fraction is linear in $\tau$:

$$
\varepsilon_{pg} = \varepsilon_{g,c} + (\varepsilon_{g,v} - \varepsilon_{g,c})\,\tau ,
$$

with `eps_g_v` and `eps_g_c` of `pyrolysis.dat`. The set-up checks $\rho_{i,0} > 0$, $\varepsilon_{i,0}$, $\varepsilon_{g,v}$, $\varepsilon_{g,c} \in (0, 1]$, $\sum_i \varepsilon_{i,0} + \varepsilon_{g,v} \le 1$, $F_{ij} \in [0, 1]$ with $\sum_j F_{ij} \le 1$, and $A$, $E$, $m$, $T_r \ge 0$.

### Kinetics

Each reaction follows an Arrhenius law with a temperature exponent and a threshold:

$$
\frac{\mathrm{d}\chi_r}{\mathrm{d}t} = \kappa_r(T)\,(1 - \chi_r)^{m_r},
\qquad
\kappa_r(T) =
\begin{cases}
A_r\,T^{\,b_r}\,\exp\!\bigl(-E_r/(\Re\,T)\bigr) & T \ge T_r \\
0 & T < T_r
\end{cases}
$$

| Symbol | Key | Quantity | Unit |
|--------|-----|----------|------|
| $A_r$ | `A[i][j]` | pre-exponential factor ($\kappa_r = 0$ when $A_r \le 0$) | $\mathrm{s}^{-1}\,\mathrm{K}^{-b_r}$ |
| $E_r$ | `E[i][j]` | activation energy | J/mol |
| $m_r$ | `m[i][j]` | reaction order | — |
| $b_r$ | `n[i][j]` | temperature exponent | — |
| $T_r$ | `T[i][j]` | threshold temperature | K |

with $\Re = 8.314471469$ J/(mol·K). The key `h[i][j]` (a heat of reaction) is read and not used: see [Stored Energy](#stored-energy).

**Integration.** At frozen $\kappa$ the rate law has a closed-form solution. Over a step $\Delta t$ from $\chi^n$:

$$
(1 - \chi)^{1 - m} = (1 - \chi^n)^{1 - m} - (1 - m)\,\kappa\,\Delta t \quad (m \ne 1),
\qquad
1 - \chi = (1 - \chi^n)\,e^{-\kappa\,\Delta t} \quad (m = 1).
$$

For $m < 1$ a right-hand side $\le 0$ gives $\chi = 1$. The result is kept in $[\chi^n, 1]$: the reactions are irreversible. Each Runge–Kutta stage evaluates this propagator from the start-of-step $\chi^n$ over the full $\Delta t$, with $\kappa$ at the stage temperature, and folds the increments with the weights of the scheme (see [Time Integration](time-integration.md#time-step-with-pyrolysis)).

**Gas production.** One pyrolysis gas is produced, at the rate

$$
\pi = \sum_r \varepsilon_{i,0}\,\rho_{i,0}\,F_r\,\frac{\mathrm{d}\chi_r}{\mathrm{d}t} = -\frac{\mathrm{d}\rho_s}{\mathrm{d}t}
\quad \text{[kg/(m³·s)]}.
$$

There are no gas species: the gas composition is the one the gas table was computed for.

---

## Stored Energy

FUSS integrates the **total stored energy** per unit volume. Without gas storage (the type-1 closure) it is

$$
E(T, \chi) = \tau\,\rho_v\,h_v(T) + (1 - \tau)\,\rho_c\,h_c(T) \quad \text{[J/m³]},
$$

with $h_v$, $h_c$ the enthalpy columns of `virgin.dat` and `char.dat` [J/kg]. The type-2 closure adds the stored gas ([Stored gas and energy](#stored-gas-and-energy)). The energy balance of a cell is

$$
\frac{\mathrm{d}}{\mathrm{d}t}\!\int_{\Omega_i}\!E\,\mathrm{d}V
= \oint_{\partial\Omega_i}\! k_s\,\nabla T\cdot\hat{\mathbf n}\,\mathrm{d}A
\;-\; \oint_{\partial\Omega_i}\! h_g\,\dot{\mathbf m}_g\cdot\hat{\mathbf n}\,\mathrm{d}A
\;+\; \dot q_\text{vol}\,V_i ,
$$

where $\dot{\mathbf m}_g$ is the gas mass flux of the closure and $h_g$ the gas enthalpy. The discrete form is that of [Governing Equations](governing-equations.md#finite-volume-integral-form): the residual $R_i$ collects the conduction and the gas-enthalpy flows leaving the cell, a stage advances $E_i$ by $\Delta t\,(-R_i/V_i + \dot q_{\text{vol},i})$, and $T_i$ is recovered by inverting $E$ at the new $\chi_i$.

**No decomposition-heat source.** The equation has no source term for the heat of pyrolysis. At constant $T$,

$$
\frac{\partial E}{\partial \rho_s}\Big|_T = \frac{\rho_v\,h_v - \rho_c\,h_c}{\rho_v - \rho_c} \equiv \bar h ,
$$

so the decomposition of a mass $\delta m$ of solid releases $\bar h\,\delta m$ from the stored energy while the gas it produces carries $h_g\,\delta m$ out through the gas-enthalpy flux. The net heat absorbed per kilogram is $h_g - \bar h$, set by the tables alone. The three enthalpy tables must therefore share one datum (absolute enthalpies, formation included); `h[i][j]` of `pyrolysis.dat` is ignored. Since $\partial E/\partial t = \rho_s c_{p,s}\,\partial T/\partial t + \bar h\,\partial\rho_s/\partial t$, this is the temperature form of [1], whose right-hand side carries $-\dot E_p = \bar h\,\pi$, written for $E$.

**Inversion.** $T$ follows from $E$ by a bracketed Newton iteration at fixed $\chi$ (and, under type 2, fixed $p$ and $G$).

- The bracket is the intersection of the temperature ranges of the solid and gas tables of the material, one row in from each end. An energy outside the bracket stops the run (exit code 1) and names the cell; so does an initial temperature outside it.
- The temperature range $[T_\text{min}, T_\text{max})$ of `properties.dat` still bounds every cell: an initial temperature, or a temperature after any stage, outside it stops the run as in the plain solver ([Admissible updates](time-integration.md#admissible-updates)). `properties.dat` must therefore cover the bracket.
- The slope is that of the 1 K table segment that holds $T$, taken from the $h$ columns. Under type 1 $E$ is piecewise linear in $T$, so a Newton step from the right segment lands on the root.
- A Newton candidate outside the bracket is replaced by bisection. The iteration ends at a step $\le 10^{-10}$ K or at $\lvert E(T) - E\rvert \le 10^{-13}\lvert E\rvert$; 50 iterations without convergence stop the run.

---

## Effective Properties

### Solid

The conductivity and the specific heat blend the virgin and char tables, both at the cell temperature, with the virgin mass share:

$$
X_s = w_v\,X_v(T) + (1 - w_v)\,X_c(T),
\qquad
w_v = \frac{\tau\,\rho_v}{\max(\rho_s,\, \rho_c)},
\qquad X \in \{c_p,\ k\}.
$$

The surface emissivity and absorptivity of an [ablating wall](ablation.md) blend the same way, at the wall temperature.

- For $\rho_s \ge \rho_c$, $\rho_s\,c_{p,s} = \tau\rho_v c_{p,v} + (1-\tau)\rho_c c_{p,c}$, which is $\partial E/\partial T$ when the $h$ columns integrate the $c_p$ columns. The loader does not compare them; the inversion uses the $h$ columns.
- $k_s$ is the effective conductivity of the porous medium, gas included. It is a scalar: `ki`, `kj` and `kk` must be equal (orthotropic conductivity is refused).
- $\rho_s$, $c_{p,s}$ and $k_s$ are cached per cell after every update. The conduction fluxes use $k_s$; the time step uses $k_s/(\rho_s c_{p,s})$.

**Face conductivity.** Under pyrolysis the conductivity of a face is not the one of [Face Conductivity](numerics.md#face-conductivity):

| Face | Conductivity |
|------|--------------|
| interior | $\tfrac12\bigl(k_{s,L} + k_{s,R}\bigr)$, also between two materials |
| block connection (101) | $\tfrac12(k_{s,P} + k_{s,N})$ for one material; the two half cells in series for two |
| wall | $k_s$ of the cell |

Inside a block a material change therefore takes the arithmetic mean of the two conductivities, not the series value.

### Gas

`gas.dat` tabulates the molar mass $W$, the enthalpy $h_g$ and the viscosity $\mu$ against $T$ at one or more pressure levels $p_1 < \dots < p_{N_p}$. A lookup is linear in $T$ between the 1 K rows and linear in $p$ between the two bracketing levels; a state outside the table, in $T$ or outside $[p_1, p_{N_p}]$, stops the run (exit code 1). With a single level the table holds at that pressure only: any other pressure (`p-gas`, a cell or a wall pressure) stops the run. The gas is ideal: $\rho_g = p\,W/(\Re\,T)$.

### Thermal equilibrium

Solid and gas share one temperature per cell. The gas enthalpy at an interior face is taken at a temperature of the solid (the face mean, or the cell temperature), and at a wall face at the wall temperature. There is no exchange term between a solid and a gas temperature. The heat capacity of the gas enters the stored energy only under type 2 with `quasi-steady = false`.

---

## Type-1 Gas Closure

`gas-closure = type1` (the default) assumes **quasi-steady, one-dimensional** gas flow: the gas produced in a cell leaves at once, along the index direction `outflow-direction`, through the end face of the block. The gas is evaluated at the fixed pressure `p-gas`.

Along each grid line of that direction the cells are numbered $s = 1 \dots N$ from the sealed end. The gas mass flow crossing the downstream face of cell $s$ is the running integral

$$
\dot M_s = \dot M_{s-1} + \pi_s\,V_s,
\qquad \dot M_0 = 0 \quad \text{[kg/s]},
$$

and the cell reports the mass flux $\dot m_{pg} = \dot M_s / A_{s+1/2}$. $\pi_s$ is the production of the stage: the increment of $\chi$ the stage applies, divided by $\Delta t$, so the line carries exactly the mass the solid loses. The gas enthalpy enters the residual of cell $s$ as

$$
R_s \mathrel{+}= \dot M_s\,h_g\bigl(p_\text{gas}, T_{s+1/2}\bigr) - \dot M_{s-1}\,h_g\bigl(p_\text{gas}, T_{s-1/2}\bigr),
$$

with $T_{s+1/2} = \tfrac12(T_s + T_{s+1})$ at an interior face and, at the exit face, the wall temperature its boundary condition computed in the stage (the cell temperature if no wall condition set one: a coupled face, 103). Summed over the line, the flows telescope: the line loses $\dot M_N\,h_g(T_w)$. The discretisation is given in [Spatial Discretisation](numerics.md#type-1-enthalpy-along-the-outflow-lines).

**What type 1 cannot represent.**

- No gas pressure and no gas storage: the gas leaves in the step it is produced, at `p-gas`.
- No flow other than along the index lines of `outflow-direction`: no lateral spreading, no flow around an impermeable insert, no exit through any other face.
- The face opposite the exit is impermeable by construction.
- One block only: the integral restarts at every block, so the set-up refuses block connections (101) and periodic faces (201). The exit face must carry a wall condition (301–307) or be a coupled face (103), and an ablating wall (306, 307) must be that face.
- The central face temperature is free of oscillations while the cell Péclet number $\dot m_{pg}\,c_{p,g}\,\Delta x/k_s$ stays below 2; the code does not check it.

The closure is exact only for a one-dimensional, quasi-steady flow with an impermeable back face.

---

## Type-2 Gas Closure

`gas-closure = type2` solves the gas pressure $p$ (the field `p_pg`) from the gas mass balance with Darcy's law. Its keys are those of `[FUSS-Darcy]` ([Input File](../user/input.md#type-2-gas-closure-darcy)).

### Mass balance and Darcy's law

$$
\frac{\partial(\varepsilon_{pg}\,\rho_g)}{\partial t} + \nabla\cdot\dot{\mathbf m}_g = \pi,
\qquad
\dot{\mathbf m}_g = -\,\frac{\rho_g}{\mu}\,\mathbf K\,\nabla p ,
\qquad
\mathbf u_g = \frac{\dot{\mathbf m}_g}{\varepsilon_{pg}\,\rho_g} .
$$

With $\varepsilon_{pg}\rho_g = \eta\,p$ this is the pressure equation

$$
\frac{\partial(\eta\,p)}{\partial t} - \nabla\cdot\bigl(\boldsymbol\Gamma\,\nabla p\bigr) = \pi,
\qquad
\eta = \frac{\varepsilon_{pg}\,W}{\Re\,T},
\qquad
\boldsymbol\Gamma = \frac{\rho_g}{\mu}\,\mathbf K_\text{eff}.
$$

**Permeability.** $\mathbf K$ is a symmetric $3\times3$ tensor in the axes of the mesh, blended in $\tau$:

$$
\mathbf K = \mathbf K_c + (\mathbf K_v - \mathbf K_c)\,\tau,
\qquad
\mathbf K_\text{eff} = \mathbf K + \frac{\mu}{p}\,\bigl[\mathbf B_c + (\mathbf B_v - \mathbf B_c)\,\tau\bigr]
\quad \text{(Klinkenberg)} .
$$

`K_v`, `K_c` are required (symmetric, positive diagonal) [m²]; the Klinkenberg term is added with `klinkenberg = true`, which requires `B_v`, `B_c` [m²/s]. The face fluxes take $\boldsymbol\Gamma\hat{\mathbf n}$ for $\hat{\mathbf n}\cdot\boldsymbol\Gamma$, so $\mathbf B_v$, $\mathbf B_c$ must be symmetric as well; the set-up does not check them. Without it $\mathbf K_\text{eff} = \mathbf K$. A tensor that follows the grid lines is not available (`axes-follow-grid` is refused).

### Stored gas and energy

The stored energy gains the energy of the gas in the pores:

$$
E = \tau\,\rho_v\,h_v(T) + (1 - \tau)\,\rho_c\,h_c(T) + G\,h_g(p, T) - \varepsilon_{pg}\,p ,
$$

and the gas-enthalpy flux of the energy balance is that of the Darcy mass flux, $h_g\,\dot{\mathbf m}_g = -h_g\,\boldsymbol\Gamma\,\nabla p$.

$G$ [kg/m³] is the stored gas mass. It is a **state of the scheme**, not $\varepsilon_{pg}\rho_g(p, T)$ re-evaluated from the equation of state: the pressure solve sets it ($G^{n+1} = \eta^\ast p^{n+1}$, below) and the end-of-step correction closes it. The gas mass in $E$ is then exactly the mass the mass balance accounted for, and no mismatch between them is ever multiplied by the absolute gas enthalpy.

With `quasi-steady = true` the storage is dropped from both equations: $\partial(\eta p)/\partial t$ from the mass balance and $G h_g - \varepsilon_{pg} p$ from $E$.

### Discrete pressure equation

The pressure is solved **once per step**, before the Runge–Kutta stages, by backward Euler with coefficients lagged at the start-of-step state $(p^n, T^n, \chi^n)$:

$$
\bigl(\eta^\ast\,p_i^{n+1} - G_i^n\bigr)\,\frac{V_i}{\Delta t}
= \sum_f F_f\bigl(p^{n+1}\bigr) + \pi^\ast_i\,V_i .
$$

- $F_f$ is the gas mass flow into the cell through face $f$, with $\boldsymbol\Gamma$ from $(p^n, T^n, \chi^n)$ (see [Spatial Discretisation](numerics.md#type-2-darcy-face-flux)).
- $G^n$ is the gas the previous step left stored.
- $\eta^\ast = 2\eta^n - \eta^{n-1}$ is the storage coefficient extrapolated to the end of the step. Where it falls below $\eta^n/2$ it is replaced by $\eta^n$ for that step. At constant $\varepsilon_{pg}$ and $W$ this happens when the temperature rose by more than half over the previous step ($T^n > 1.5\,T^{n-1}$); the unlimited extrapolation would turn negative once the temperature doubles.
- $\pi^\ast$ is the production over the step from the closed-form kinetics with $\kappa$ at the extrapolated mid-step temperature $T^\ast = T^n + \tfrac12(T^n - T^{n-1})$. The test key `pi-source` adds a constant production.
- The first step of a new run has no $T^{n-1}$ or $\eta^{n-1}$: it is run twice, and the second pass takes them from the first ($T^{n-1} = 2T^n - T^{n+1}$, $\eta^{n-1} = 2\eta^n - \eta^{n+1}$). A restart reads them from the solution file.

After the solve $G^{n+1} = \eta^\ast p^{n+1}$, the cell velocity $\mathbf u_g = -\boldsymbol\Gamma^n\nabla p^{n+1}/(\eta^n p^{n+1})$ is computed for output (central differences, which read the mirrored ghost cells of an open face, not $p_w$), and every open face records its outward mass flux $\dot m_{pg}$, which an [ablating wall](ablation.md) blows. Every Runge–Kutta stage reads $p^{n+1}$ as frozen data.

Under `quasi-steady = true` the storage term is absent. The difference between the gas the solve expels and the gas the stages release is carried to the next step's source (never below `pi-source`), so the gas mass closes with a one-step lag. A sealed quasi-steady domain has no unique pressure and is refused.

### Boundary conditions of the pressure

| Face | Condition |
|------|-----------|
| open | $p = p_w$ at the face (Dirichlet, half a cell from the centre) |
| every other wall, symmetry, null | impermeable, $\dot{\mathbf m}_g\cdot\hat{\mathbf n} = 0$ |
| block connection (101) | $p$ and the gas mass flow continuous |
| chimera (102, 104) | refused |

- A face is open when its code is in `open-bc-types` (default `303 305 306 103`, any of 301–306 and 103). An ablating wall (306, 307) is always open.
- A coupled face (103) is open at `p-wall`: no gas pressure is exchanged with the other solver, and the gas that leaves through it is not passed to it (its mass and enthalpy leave the energy budget of FUSS only).
- $p_w$ is the optional fourth value of a 303 or 305 line, else `p-wall`; an ablating wall takes the $p_w$ of its env table or exchange.
- The set-up refuses a case without an open face while gas is produced.

### Pressure solver

`solver = tridiag` (default) solves a block whose open faces lie in one index direction line by line along it (Thomas algorithm). A block without open faces whose block connection ends the lines of such a block continues those lines. Every other block, and every block with `solver = pseudo`, takes damped Jacobi sweeps, $p \leftarrow p + \omega\,r/D$ with $\omega$ = `pseudo-vnn`.

- Lines are coupled by line Jacobi: the conductance of the faces across the lines is on the diagonal, the neighbour's pressure of the current iterate on the right-hand side. Blocks are coupled the same way through their connection ghosts.
- After each sweep an additive correction removes the error that is uniform across the lines. Line-implicit blocks joined through 101 faces share one correction when, consistently over all their joining faces, their lines continue across the faces or lie side by side; otherwise each block corrects alone. A singular summed system skips the correction. Lines that carry the same data are solved to round-off in one sweep.
- The sweeps stop, after at least one update, when the largest cell mass residual is at most `tol` times the largest flux, source or storage term, or below the round-off floor $64\,\epsilon\,\max_i D_i\lvert p_i\rvert$. `max-iter` sweeps without convergence stop the run (exit code 1). `OUTPUT/<prefix>darcy.txt` records the sweeps and the residual of every step.

### End-of-step mass correction

The pressure solve took the production $\pi^\ast\Delta t$ per unit volume (stored or expelled); the stages decompose the solid by $\sum_r \varepsilon_{i,0}\rho_{i,0}F_r(\chi_r^{n+1} - \chi_r^n)$. After the last stage the difference

$$
\delta m = \sum_r \varepsilon_{i,0}\rho_{i,0}F_r\,(\chi_r^{n+1} - \chi_r^n) - (\pi^\ast - \pi_\text{test})\,\Delta t ,
$$

with $\pi_\text{test}$ the `pi-source` test key, goes into the stored gas, $G^{n+1} \leftarrow G^{n+1} + \delta m$, and $T$ is inverted again at constant $E$. The gas mass, $\int(\rho_s + G)\,\mathrm{d}V$ plus the mass that left, is then conserved to round-off.

When the correction would leave $G^{n+1} \le 0$ (a prediction that crossed a threshold temperature the stages did not reach), the solid supplies the deficit, bringing the stored gas towards a tenth of its start-of-step value $G^n$ as far as its remaining mass allows: the reactions that ran in the step, or can run at $T$, advance further in proportion to the mass each has left (every reaction when those cannot cover it). A solid that cannot keep $G^{n+1}$ positive stops the run.

### Accuracy

Measured on the exact and literature cases of `test/pyrolysis/`:

| Case | Observed order |
|------|----------------|
| steady sealed profile with a source (`darcy/sealed_profile`, 20–160 cells) | 1.94, 1.97, 1.98 in space |
| the same with Klinkenberg | 1.97, 1.99, 1.99 in space |
| transient pressure diffusion (`exact/pressure_transient`) | 2.08 in space; first order in time (backward Euler) |
| steady transpiration through two open faces (`exact/transpiration`) | 2.00 |
| sharp decomposition front, similarity solution (`exact/sharp_front`) | 1.04 |

The pressure is first order in time; the coefficients $\boldsymbol\Gamma$ and $\eta$ are lagged by one step.

---

## Combinations Refused at Set-up

| Closure | Refused |
|---------|---------|
| both | `integration-variables = prim`; `time-accurate = false`; IRS; more than one grid level; chimera faces (102, 104) |
| type 1 | more than one block, block connections (101), periodic faces (201); no `outflow-direction`; an exit face that is neither a wall (301–307) nor a coupled face (103); an ablating wall off the exit face |
| type 2 | `[FUSS-Pyrolysis] p-gas`; no open face while gas is produced; no open face with `quasi-steady = true`; a material without `K_v`, `K_c`; `klinkenberg = true` without `B_v`, `B_c`; `axes-follow-grid = true` |

`[FUSS-Darcy]` without `gas-closure = type2` is refused as well.

---

## References

1. M. Rotondi, *Modeling of Conjugate Heat Transfer Including Ablation and Pyrolysis for Propulsive Applications*, PhD thesis, Sapienza University of Rome, 2023, ch. 3.
2. J. Lachaud, N. N. Mansour, "Porous-material analysis toolbox based on OpenFOAM and applications," *J. Thermophys. Heat Transfer* 28 (2), 2014.
