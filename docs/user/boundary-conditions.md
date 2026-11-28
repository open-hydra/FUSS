# Boundary Conditions

FUSS reads the boundary conditions from `INPUT/bc.txt`, and those of the coarse multigrid levels from `INPUT/bc2.txt`, `INPUT/bc3.txt`, .... The file holds one record for each cell face on the six faces of every block. ATLAS (BCB) writes it from the `[BCB-Block<n>]` sections of `input.ini` (see [ATLAS Input](#atlas-input)); FUSS does not read those sections.

---

## Face Numbering

| Face | Location | Cells on the face |
|:----:|----------|-------------------|
| 1 | $i = 1$ (imin) | $n_j \times n_k$ |
| 2 | $i = n_i$ (imax) | $n_j \times n_k$ |
| 3 | $j = 1$ (jmin) | $n_i \times n_k$ |
| 4 | $j = n_j$ (jmax) | $n_i \times n_k$ |
| 5 | $k = 1$ (kmin) | $n_i \times n_j$ |
| 6 | $k = n_k$ (kmax) | $n_i \times n_j$ |

On a 2-D mesh ($n_k = 1$) faces 5 and 6 must be zero-flux faces (code `0` or `300`). When `face5` and `face6` are left out, ATLAS writes `null` for them on a plane mesh and `axisymmetric` (code 200, which FUSS treats as `0`) on an axisymmetric one.

---

## BC Codes

| Code | ATLAS keyword or keys | Behaviour | Data line |
|:----:|-------|-----------|-----------|
| `0` | `null` | Zero heat flux (adiabatic); the ghost cells mirror the cells inside | none |
| `300` | `symmetry` | Same as `0` | none |
| `101` | `connection` | Block connection | `bs is js ks fs d11 d12 d21 d22` |
| `102` | `chimera` | Overset connection | `n1 n2`, then `n1 + n2` lines `bD iD jD kD vf` |
| `103` | `connection` between blocks of different phases | Coupling with the hydra fluid solver on matching faces | as `101` |
| `104` | — | Coupling with the hydra fluid solver, donors in the other solver | as `102` |
| `201` | `type = periodic`, `blocks`, `faces` | Periodic face: rotational or translational periodicity | as `101` |
| `301` | `q`, `q-time-file` | Prescribed heat flux $q_w$ into the solid | `qw`, or a time table |
| `302` | `T`, `T-time-file` | Prescribed wall temperature $T_w$ | `Tw`, or a time table |
| `303` | `hconv`, `qrad`, `Tref` | Convection plus a prescribed flux $q_\text{rad}$ | `hconv, qrad, Tref`; optionally a fourth value `p_w` [Pa], the gas pressure of the face under the type-2 gas closure ([Input File](input.md#type-2-gas-closure-darcy)), ignored otherwise, > 0 when given |
| `304` | `eps`, `Tref` | Radiation | `eps, Tref` |
| `305` | `hconv`, `eps`, `Tref` | Convection and radiation | `hconv, eps, Tref`; optionally a fourth value `p_w` [Pa], as for `303` |
| `306` | — (written by the case generators, not by ATLAS) | Ablating wall: surface energy balance with B′ tables | `lambda, T_inf, <bprime-file>, <env-file>`: see [Ablating wall](#ablating-wall-306-and-coupled-ablating-wall-307) |
| `307` | — | Coupled ablating wall: 306's balance, exchange values from hydra | `bs is js ks fs d11 d12 d21 d22 T_inf <bprime-file>` (11 values, no `lambda`) |

Units: W/m² for `qw` and `qrad`, K for `Tw` and `Tref`, W/(m²·K) for `hconv`; `eps` is the emissivity.

- **`0` and `300`** add no flux through the face.
- **`101`**: the face of cell $(i,j,k)$ of block `b` touches face `fs` of cell `(is, js, ks)` of block `bs`. `(d11, d12)` and `(d21, d22)` are the index steps on face `fs` that correspond to the first and second tangential index directions of face `f`. The interface temperature follows from the continuity of the heat flux, each side with its own conductivity and metric over half a cell.
- **`102`**: `n1` and `n2` donor cells fill the first and second ghost layer; each donor line gives its block, cell indices and volume fraction. The ghost temperature is the one of the volume-weighted enthalpy of the donors, which must share one material. FUSS stops if the fractions of a layer sum to less than 0.5.
- **`103`, `104`**: the flux is set by hydra; FUSS on its own adds nothing, and its ghost cells mirror the cells inside. Under pyrolysis a `103` face is accepted and the `104` face refused; the pyrolysis gas can leave through a `103` face, and it is not passed to the gas solver ([Pyrolysis](input.md#pyrolysis)).
- **`201`**: a periodic face, as the two side faces of one channel of a cooling jacket (the mid-planes of the ribs between two channels: rotational periodicity about the axis) or the two ends of one period of a translationally periodic part. The record is that of a `101`: face `f` of cell $(i,j,k)$ of block `b` takes as ghost cells the cells across face `fs` of cell `(is, js, ks)` of block `bs`, and the record of that cell points back at it (each of the two faces is the boundary of the other). FUSS treats it as a block connection in everything (ghost temperatures and materials, one flux per face computed by its owner and given to both cells, the exchanges between ranks), except the positions of the cells of the other face, which lie one period away: they are moved by the rigid motion $\mathbf x \mapsto R\,\mathbf x + \mathbf t$ that takes face `fs` onto face `f` (the ghost nodes, hence the ghost metric of `diffusive-metric = mean`, and the ghost centre and its tangential neighbours of the stencil metric). The face is then, for both metrics, the interior face it is in the full mesh. No axis, angle or period is given: FUSS fits the motion of each pair of faces to the centres of the faces of its records and to the points half a cell across them (B. K. P. Horn, J. Opt. Soc. Am. A 4, 1987) and prints it at set-up, one line per direction of each pair (`Periodic faces (201)`: rotation angle and axis, translation, misfit in units of the half cell). FUSS stops at set-up when the misfit exceeds $10^{-3}$ of the smallest half cell (the faces are not images of each other with the index map of the records: a wrong pair of faces, or the cells of the two faces in another order), when the records do not fix the motion (two one-cell faces), when a face is periodic one way only, or when the two directions of a pair do not invert each other. A `201` whose partner is a `101` (or the reverse) stops the run as a connection that does not point back, and so does, when the file is read, a `201` whose `d11 d12` and `d21 d22` are not two unit index steps along the face (BCB before the fix of its periodic orientation rounded the products of the bases of the two turned faces and wrote steps up to 6 on stretched cells). AFUS writes the same records from its periodic BC (type 6 of `bcdef_file`), ATLAS/BCB from a section `type = periodic` (below). On a moving mesh (`[FUSS-MeshMotion]`) a pair keeps the motion fitted at set-up, and a periodic domain can follow only a mesh motion that takes a node of one face and its image to images: every step FUSS measures the misfit of every pair against that motion as at set-up (the centres of the faces and the points half a cell across them) and stops when it has grown past the pair's set-up misfit by more than the set-up bound, $10^{-3}$ of the smallest half cell, naming the law: the problem is no longer periodic. `translation` must move along the axis of a rotational pair (checked at set-up: the copies of the domain would move by $R\,\mathbf v$, not by $\mathbf v$); `prescribed` must give a node and its image images (a displacement along the axis that depends on the position along it, or one that vanishes on the faces, `taper-to-boundary`); `nodeshift` joins a receding node and its image into one surface node. Under `[FUSS-Pyrolysis] gas-closure = type2` every material's permeability tensors must be invariant under the rotation of the pair (checked at set-up); `type1` takes no periodic face.
- Any other code is read without a data line and acts as `0`. If such a record is followed by a data line, that line is read as the next record and the run stops; the message names the first record with a code outside this table.

---

## Wall Model (301–305)

The heat flux into the solid through a wall face is

$$
q_\text{in} = G\,(T_w - T_P) + C,
$$

where $T_P$ is the temperature of the cell, $T_w$ the wall temperature, $G$ the conductance of the half cell between the cell centre and the face ($G = 2k/\Delta x$ on a Cartesian mesh, $\Delta x$ the cell width normal to the wall, $k$ at $T_P$), and $C$ the cross-diffusion term from the temperature differences between the neighbouring cells along the face (zero on an orthogonal mesh). The wall condition prescribes one of $T_w$, $q_\text{in}$:

| Code | Prescribed |
|:----:|------------|
| `301` | $q_\text{in} = q_w$ |
| `302` | $T_w$ |
| `303` | $q_\text{in} = h_\text{conv}\,(T_\text{ref} - T_w) + q_\text{rad}$ |
| `304` | $q_\text{in} = \varepsilon\,\sigma\,(T_\text{ref}^4 - T_w^4)$ |
| `305` | $q_\text{in} = h_\text{conv}\,(T_\text{ref} - T_w) + \varepsilon\,\sigma\,(T_\text{ref}^4 - T_w^4)$ |

with $\sigma = 5.67\times10^{-8}$ W/(m²·K⁴). For `302` the relation gives $q_\text{in}$; for the others FUSS solves the prescribed $q_\text{in}(T_w) = G\,(T_w - T_P) + C$ for $T_w$ by Newton's method, exact in one step without radiation.

- $q_\text{in}$ and $q_w$ are positive into the solid on every face, whatever its orientation: `q = 1.0e5` heats the solid.
- `303` with `qrad = 0` is pure convection.
- The ghost cells hold $2\,T_w - T$ of the mirrored cells, exact for a linear field.
- `OUTPUT/wall.*` holds $T_w$ (`Tw`) and $q_\text{in}$ (`qw`) of every wall face.
- `null` and `q = 0` both give zero flux; with `q = 0` the face is a wall: it appears in `wall.*` and its ghost cells are $2\,T_w - T$.

The [skewed plate](../vv/2D-steady-state-skewed-plate.md) test checks the model on a sheared mesh with walls 301, 302, 303 and 305.

---

## Time Tables (301, 302)

When the data line of a `301` or `302` record is not a number, FUSS reads it as the name of a time table for $q_w$ or $T_w$:

```
periodic
0.000000e+00  2.000000e+03
1.000000e+00  2.078459e+03
2.000000e+00  2.156434e+03
...
```

- Two columns, time (s) and value. An optional first line containing `periodic` repeats the table with period $t_\text{last} - t_\text{first}$.
- Blank lines are skipped, also at the end of the file. Error messages give the line number in the file.
- Linear interpolation; before the first time the first value, after the last time the last value. Steady runs (`time-accurate = false`) run at time −1 and get the first value.
- Each Runge–Kutta stage reads the table at its own time (`t^n + c_k Δt`: RK3 `0, 1, ½`), so time-dependent boundary data keep the order of the scheme; the wall output reports the value at the output time.
- A steady run (`time-accurate = false`, and hydra not time-accurate) stays at time −1: each table gives its value there for the whole run (its first value; a `periodic` table its value at −1 within the period), and set-up prints a `[WARNING]` line saying so.
- The file is opened relative to the directory FUSS runs in; the name is the whole data line (up to 256 characters, a trailing comma dropped), so it may point into a subdirectory, e.g. `tables/Tw.dat`. A missing, empty or non-increasing table stops the run with a message.

---

## `bc.txt` Format

Each record is a header line followed by the data line(s) of its code (none for `0` and `300`):

```
<b> <i> <j> <k> <f> <code>
<data line>
```

`b` is the block, `(i, j, k)` the cell (1-based) and `f` its face (1–6). Values are free format, separated by blanks or commas.

- The number of records must equal the number of boundary faces, $\sum_\text{blocks} 2\,(n_i n_j + n_j n_k + n_i n_k)$, or FUSS stops with `Boundary conditions number (...) is different than the one of the initial conditions (...)`.
- Every block, face and cell index of a record is checked on every grid level: `b`, `i`, `j`, `k`, `f`; the partner `bs`, `is`, `js`, `ks`, `fs` of a connection (`101`) or of a periodic face (`201`, the record of a `101`); every donor `bD`, `iD`, `jD`, `kD` of a chimera record (`102`). A value out of range stops every rank with exit code 1, for example `bc.txt line 22 (record 11, type 101): partner block bs = 9 is outside the admissible range 1 to 2`. The partners of `103`, `307` and `104` belong to the other solver and are not checked.
- Records are grouped by block (1, 2, ...) and, within a block, by face (1 to 6), the order ATLAS writes; the wall output relies on it.

Records from the test cases, as ATLAS writes them:

`test/transient/oscillating_Rod`: a fixed and a time-table wall temperature.

```
       1       1       1       1       1     302
    0.200000E+04,
       1     100       1       1       2     302
 Tw_time.dat,
```

`test/steady-state/radiative_Rod`: a zero-flux face and a radiating wall.

```
       1       1       1       1       3       0
       1     100       1       1       2     304
    0.980000E+00,    0.300000E+03,
```

`test/steady-state/convective_Plate`: convection.

```
       1      60       1       1       2     303
    0.750000E+03,    0.000000E+00,    0.273150E+03,
```

`test/transient/multimat_Plate`: face 2 of block 1 connected to face 1 of block 2.

```
       1      15       1       1       2     101
       2       1       1       1       1       1       0       0       1
```

`test/numerics-features/chimera`: two donors for each ghost layer.

```
       1     100       2       1       2     102
       2       2
       2       1       1       1    0.4285714290E+00
       2       1       2       1    0.5714285710E+00
       2       2       1       1    0.4285714290E+00
       2       2       2       1    0.5714285710E+00
```

### Ablating wall (306) and coupled ablating wall (307)

ATLAS does not write these two codes yet; the case generators under `test/pyrolysis/` write `bc.txt`
directly. Both faces solve the surface energy balance of the material response (thesis eq. 3.49, PATO's
`BprimeBoundaryConditions`) for the wall temperature, with the char and pyrolysis-gas blowing of an equilibrium or
finite-rate B′ table. They need `[FUSS-Pyrolysis]` when the chemistry is on, and they set the gas pressure of the face
under the type-2 closure (open faces). Keys and outputs: [Input File](input.md#ablating-wall-306-307).

**306 line 2**, list-directed: `lambda  T_inf  <bprime-file>  <env-file>`.

| Field | Description |
|-------|-------------|
| `lambda` | blowing-correction exponent (`β = ζ/(e^ζ − 1)`, `ζ = 2λ(ṁ_pg + ṁ_ch)/(ρ_e u_e C_H0)`); `0` = no correction, `β = 1` |
| `T_inf` | background temperature of the re-radiation, K |
| `<bprime-file>` | bare file name under `INPUT/`, or `none` when every env row has `chem = 0` |
| `<env-file>` | bare file name under `INPUT/`: the time table of the face |

File names are bare (no `/`, `,` or blanks; at most 255 characters): an undelimited list-directed value ends at the
first `/`.

**Env table** (`//` lines are comments): one row per time, `t p_w rhoeUeCH0 h_r q_rad chem hconv T_edge`, linear in t,
clamped outside its range, evaluated once per step at the step's end time.

| Column | Units | Meaning |
|--------|-------|---------|
| `p_w` | Pa | wall pressure: the B′ lookup, and the gas pressure of the face under type 2 |
| `rhoeUeCH0` | kg/(m²·s) | unblown heat-transfer coefficient (Stanton number times ρ_e u_e) |
| `h_r` | J/kg | recovery enthalpy |
| `q_rad` | W/m² | incident radiative flux (absorbed with the material's absorptivity) |
| `chem` | — | 1 = chemistry on (the B′ balance). Below 1 (0, or the switch-off ramp between two rows) the chemistry-off balance, PATO's `chemistryOn == 0`: convection `hconv (T_edge − T_w)`, no gas advection, radiation kept; while `0 < chem < 1` the B′_pg of the last chemistry-on solve and the blowing correction stay frozen. With `chem = 0` and a zero emissivity the face is BC 303 `hconv, T_edge` bit for bit |
| `hconv`, `T_edge` | W/(m²·K), K | the chemistry-off convection |

**B′ table** (PATO's `BprimeTable` layout): rows `p B'g B'c T h_w` (`//` comments), p-blocks of decreasing `p`, B′g
sub-blocks of decreasing B′g, `T` monotonic within a sub-block (repeated `T` allowed). The header line `// kind:
equilibrium | finite-rate | forced-equilibrium` records the provenance. A 4-D finite-rate table adds `rhoeueCM` as the
first column (`// axes: rhoeueCM p Bprime-g T`). Queries outside the table are clamped and counted in `seb.txt`, or
stop the run with `[FUSS-Ablation] bprime-out-of-range = abort` (on `p`, B′g or `rhoeueCM`; a query outside the `T` range
is always clamped). The equilibrium TACOT table of the test cases and its
provenance: `test/pyrolysis/ablation/tables/`.

**307 line 2**: the 103 connectivity record of the MOSE partner face, then `T_inf <bprime-file>`. hydra writes the
exchange values `p_w, ρ_e u_e C_H, h_r, q_rad` before every FUSS step (`bc%env_now`); the face runs with `chem = 1` and
`lambda = 0` (the CFD's coefficient is already blown). A standalone run reads those values from the
TEST ONLY `[FUSS-Ablation] coupled-env-file`, a table with the rows `t p_w rhoeUeCH h_r q_rad`, or one such table per
307 face, each opened by a line `// face b i j k f` (see [Input File](input.md#ablating-wall-306-307)).

The wall output's `qw` of a 306/307 face (and `bc%qw_last`, which hydra writes as column 7 of a 307 zone) is the conduction flux into the solid, positive when it heats the solid, on every face.

---

## ATLAS Input

ATLAS BCB builds `bc.txt` from one `[BCB-Block<n>]` section per block, which names the condition of each face:

```ini
[BCB-Block1]
face1 = Twall-hot
face2 = connection
face3 = adiabatic
face4 = convection
face5 = null
face6 = null
```

A face takes one of the keywords `null` (code 0), `symmetry` (300), `connection` (101, or 103 between blocks of different phases), `chimera` (102, with `BC-chimera = true` in `[ATLAS-Parameters]`), or the name of a section with `type = wall`. The keys present in the section select the code:

| Keys in the section | Code |
|---------------------|:----:|
| `q` (without `T`, `qrad`) | `301` |
| `q-time-file` (without `T`, `qrad`) | `301` |
| `T` (without `q`, `qrad`) | `302` |
| `T-time-file` (without `q`, `qrad`) | `302` |
| `hconv`, `qrad`, `Tref` | `303` |
| `eps`, `Tref` (without `hconv`) | `304` |
| `hconv`, `eps`, `Tref` (without `qrad`) | `305` |

ATLAS checks the rows in this order and takes the first that matches.

!!! warning "A section that matches no row is a zero-flux face"
    For example `hconv` and `Tref` without `qrad`. ATLAS writes code `0` for its faces and prints a
    warning. Pure convection needs `qrad = 0`.

Sections of the test cases:

```ini
[adiabatic]
type = wall
q    = 0.0

[Twall-hot]
type = wall
T    = 1273.15

[Tw_variable]
type        = wall
T-time-file = Tw_time.dat

[convection]
type  = wall
hconv = 750.0      ; W/(m2 K)
Tref  = 273.15     ; K
qrad  = 0.0        ; W/m2, needed for pure convection

[radiation]
type = wall
eps  = 0.98
Tref = 300.0
```

A section can be shared by several faces and blocks. ATLAS also splits a face into patches with different conditions: see the ATLAS documentation.

### Periodic faces (201)

A section `type = periodic` pairs two faces, given by their numbers in the mesh (`blocks`, the numbers of the zones of `mesh.tec`, and `faces`); both faces name the same section. ATLAS pairs the cells of the two faces in the order of their indices on the face, as AFUS does, and writes the `201` records of both faces:

```ini
[BCB-Block2]
phase = solid
face3 = per_1_6        ; side face of the sector ...
...
[BCB-Block7]
phase = solid
face4 = per_1_6        ; ... and the side face of the next sector, its image by the rotation of one sector
...
[per_1_6]
type   = periodic
blocks = 2 7
faces  = 3 4
```

The two faces must be images of each other by a rigid motion with that order of the cells (the side faces of a sector of a mesh built by rotation, or the two ends of one period): FUSS fits the motion and stops when they are not (see `201` above). A face periodic with another face of its own block gives the block twice (`blocks = 5 5`). In a mesh with the blocks of more than one phase (hydra) the records carry the numbers of the blocks in the BC file of the phase, as the connections.
