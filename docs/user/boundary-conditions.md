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
| `301` | `q`, `q-time-file` | Prescribed heat flux $q_w$ into the solid | `qw`, or a time table |
| `302` | `T`, `T-time-file` | Prescribed wall temperature $T_w$ | `Tw`, or a time table |
| `303` | `hconv`, `qrad`, `Tref` | Convection plus a prescribed flux $q_\text{rad}$ | `hconv, qrad, Tref` |
| `304` | `eps`, `Tref` | Radiation | `eps, Tref` |
| `305` | see the [warning](#atlas-input) | Convection and radiation | `hconv, eps, Tref` |

Units: W/m² for `qw` and `qrad`, K for `Tw` and `Tref`, W/(m²·K) for `hconv`; `eps` is the emissivity.

- **`0` and `300`** add no flux through the face.
- **`101`**: the face of cell $(i,j,k)$ of block `b` touches face `fs` of cell `(is, js, ks)` of block `bs`. `(d11, d12)` and `(d21, d22)` are the index steps on face `fs` that correspond to the first and second tangential index directions of face `f`. The interface temperature follows from the continuity of the heat flux, each side with its own conductivity and metric over half a cell.
- **`102`**: `n1` and `n2` donor cells fill the first and second ghost layer; each donor line gives its block, cell indices and volume fraction. The ghost temperature is the one of the volume-weighted enthalpy of the donors, which must share one material. FUSS stops if the fractions of a layer sum to less than 0.5.
- **`103`, `104`**: the flux is set by hydra; FUSS on its own adds nothing, and its ghost cells mirror the cells inside.
- Any other code is read without a data line and acts as `0`: FUSS does not reject it.

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

- Two columns, time (s) and value; an optional first line containing `periodic` repeats the table with period $t_\text{last} - t_\text{first}$.
- Linear interpolation; before the first time the first value, after the last time the last value. Steady runs (`time-accurate = false`) run at time −1 and get the first value.
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
- Every block number, face number and cell index that a record carries is checked when the file is read, on every grid level: `b` (1 to the number of blocks), `i`, `j`, `k` (1 to the cells of block `b` in that direction) and `f` (1 to 6); of a connection (`101`) `bs`, `is`, `js`, `ks` (the cells of block `bs`) and `fs` (1 to 6); of a chimera record (`102`) `bD`, `iD`, `jD`, `kD` of every donor. A value outside its range stops every rank with exit code 1 and writes no solution, with the line, the record and the range, for example `bc.txt line 22 (record 11, type 101): partner block bs = 9 is outside the admissible range 1 to 2`. The partner of a `103` and the donors of a `104` are cells of the other solver and are not tested.
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
| `eps`, `Tref` | `304` |

ATLAS checks the rows in this order and takes the first that matches. A section that matches none, for example `hconv` and `Tref` without `qrad`, gets code `0`, a zero-flux face, without a message: pure convection needs `qrad = 0`.

!!! warning "Convection with radiation (305)"
    ATLAS 7d3c72d checks `eps`, `Tref` before `hconv`, `eps`, `Tref`: a section with all three is written as `304` and `hconv` is lost. Until ATLAS is fixed, write the `305` records in `bc.txt` by hand (data line `hconv, eps, Tref`).

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
