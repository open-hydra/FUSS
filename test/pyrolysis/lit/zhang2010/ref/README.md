# Zhang (2010), Figs 2.2 and 2.3: digitized reference of the zhang2010 gate

Extracted on 2026-10-13 from the vector paths of the thesis PDF (the extraction scripts and the page SVGs are not in
the repository; they are archived with the benchmark record). Every value below comes from the drawn curves.

## Source

- **Document.** Zhenyu Zhang, *Thermo-Mechanical Behavior of Polymer Composites Exposed to Fire*, PhD dissertation, Virginia Tech, 29 June 2010.
  - File: `bench/papers/Zhang_2010_VTthesis.pdf`, 160 pages, produced from Word by PScript5 and Acrobat Distiller 7.0.
- **Figures.**
  - **Fig 2.2** (temperature histories) is on **PDF page 59**, printed page 49.
  - **Fig 2.3** (mass loss per unit area at the exposed surface) is on **PDF page 60**, printed page 50.
  - The mapping was verified from the form feeds of the text dump and the printed page numbers. The earlier renders `img/laneA/zhang_ver-059/060.png` are these two pages.
- **Case.** Section 2.4.1 and Table 2.1 (printed pp. 35-36; `txt/Zhang_2010_VTthesis.txt`, lines 2019-2075):
  - geometry: slab 0.0127 m;
  - exposed face: q_rad = 25 kW/m², ε = 1, σ = 5.67e-11 kW/m²K⁴, h = 0.01 kW/m²K, T∞ = 300 K;
  - unexposed face adiabatic; T_ini = 300 K;
  - Table 2.1: ρ0 = 1700 kg/m³, ρf = 1255 kg/m³, cp = 1100 J/kgK, cg = 3000 J/kgK, k = 0.30 W/mK, Q = -2e5 J/kg, A = 5.0e28 1/s, E = 3.62e5 J/mol, n = 4.6;
  - codes compared: UMATHT with 20 3-D solid elements, the 1-D UEL element, and the MATLAB FEM code of Summers [134].

## Files

| file | content |
|---|---|
| `zhang_fig2_2_T.csv` | `time_s` (0 to 2000, step 10) plus 9 columns: `{MATLAB,UEL,UMATHT}_{exposed,middle,unexposed}_K`. Empty cell = no drawn curve at that time. |
| `zhang_fig2_3_mloss.csv` | `time_s` (0 to 2000, step 10) plus `MATLAB_kg_s_m2`, `UEL_kg_s_m2`, `UMATHT_kg_s_m2`. |
| `zhang_fig2_2_T_unc.csv`, `zhang_fig2_3_mloss_unc.csv` | Same layout. Holds the **bound on the uncertainty** of each cell, in the same units (model below). |
| `zhang_fig2_2_vertices.csv`, `zhang_fig2_3_vertices.csv` | Raw vertices of the drawn paths, one row per vertex: `curve, stroke, kind, time_s, value, x_pt, y_pt`. `kind` is `stroke_start`, `stroke_end` (dash ends), `L_end` or `C_end` (segment ends). These are the unresampled data. |
| `zhang_fig2_2_overlay.png`, `zhang_fig2_3_overlay.png` | Overlay check. Left: the page rendered at 600 dpi, with the **CSV values** mapped back through the calibration (o = MATLAB, □ = UEL, ◇ = UMATHT). Right: four zoom windows rendered at 2400 dpi (axes stretched), which also show the raw vertices. |

The CSVs have a plain header row and no `#` comment lines, so they load directly. Provenance is in this file.

All curve values are tagged **[S]**: scaled from the vector figure, with the tick marks as the scale reference. The `_unc` values are tagged **[C]**: computed here with the formula below.

## Method

1. **Export.**
   - `pdftocairo -svg` on pages 59 and 60.
   - `pdftotext -bbox` for the tick labels.
   - The SVG paths were parsed together with their transform matrices into page points (y down).
   - Every coordinate lies on a **0.06 pt grid (1200 dpi PostScript device)**, to within 0.0025 pt (cairo's 1/256 rounding). This is the resolution limit of the source.
2. **Path identification.**
   - **Color** gives the surface (black = exposed, red = middle, blue = unexposed).
   - **Dash pattern** gives the code. Contiguous sub-paths were merged into strokes, and each stroke's on-length and gap were measured:

     | | solid | dashed | dotted |
     |---|---|---|---|
     | curves, Fig 2.2 | 1 continuous stroke | on 1.32-1.37 pt / gap 1.40-1.45 pt, duty 0.49 | on 0.32-0.35 / gap 1.30-1.34, duty 0.20-0.21 |
     | legend, Fig 2.2 | solid | duty 0.45 | duty 0.18 |
     | curves, Fig 2.3 | | duty 0.46 | duty 0.18 |
     | legend, Fig 2.3 | | duty 0.46 | duty 0.20 |

   - The legend headers (MATLAB, UEL, UMATHT) and the surface names were matched to the sample lines from the text layer, so: solid = MATLAB, dashed = UEL, dotted = UMATHT.
   - Paths inside the white legend box were excluded. Zero curve points fall inside it.
   - Fig 2.3 holds all curves in two path elements, drawn in the order MATLAB, UEL, UMATHT. They were split where the drawing jumps back to t = 0, and each group's class was confirmed from its dash duty.
3. **Calibration.**
   - **Tick marks** come from the axis path: major ticks are the long ones.
   - **Values** come from the nearest label in the text layer. The label-centre-to-tick offsets are tiny next to the tick spacing (40.7 pt and 36.8 pt in Fig 2.2; 49.2 pt and 36.3 pt in Fig 2.3):
     - Fig 2.2: x +0.40 to +0.44 pt, y +0.16 pt.
     - Fig 2.3: x -0.02 to +0.05 pt, y +0.12 to +0.18 pt.
   - **Linear 2-point map on the end ticks**: Fig 2.2 uses 0 and 2000 s, 300 and 800 K; Fig 2.3 uses 0 and 2000 s, 0.000 and 0.005.
4. **Curve evaluation.**
   - The drawn path is evaluated exactly: cubic Bezier segments at 0.005 pt steps, straight segments as drawn.
   - For dashed and dotted curves, the gaps are bridged linearly between drawn dash ends.
   - Every curve starts exactly on the tick intersection of the axis origin: (0 s, 300 K), or (0 s, 0.000) in Fig 2.3.
   - Resampling onto the 10 s grid uses linear interpolation along the evaluated path.

### Calibration residuals (2-point map, other labelled ticks)

| | Fig 2.2 | Fig 2.3 |
|---|---|---|
| time, labelled ticks 500 / 1000 / 1500 s | +0.743 / +0.719 / -0.024 s | +0.586 / +0.615 / -0.030 s |
| value, labelled ticks | 400 K +0.066, 500 K +0.132, 600 K -0.132, 700 K -0.055 K | 0.001: +2e-8; 0.002: +1.66e-6; 0.003: -1.66e-6; 0.004: -1.64e-6 |
| unlabelled minor ticks (max) | 1.11 s, 0.16 K | 1.21 s, 2.5e-6 |
| least-squares fit on all majors vs 2-point fit (max) | 0.44 s, 0.05 K | 0.36 s, 8.6e-7 |

All residuals are at most about 1 device quantum, because the ticks are rounded to the same 0.06 pt grid as the curves.

### Pixel size in data units

| | device quantum 0.06 pt (resolution of the data) | 600 dpi pixel (0.12 pt) | 2400 dpi pixel (zooms) |
|---|---|---|---|
| Fig 2.2 time | 0.736 s | 1.473 s | 0.368 s |
| Fig 2.2 T | 0.163 K | 0.326 K | 0.081 K |
| Fig 2.3 time | 0.610 s | 1.220 s | 0.305 s |
| Fig 2.3 mass loss | 1.65e-6 kg/(s m²) | 3.31e-6 | 8.3e-7 |

### Uncertainty model (the `_unc.csv` files: bounds, not 1-sigma)

`u = u_v + |dv/dt| * u_t + gap + spread/2`

- `u_t = max(1 quantum, max time residual)` = 0.743 s (Fig 2.2) or 0.615 s (Fig 2.3).
- `u_v = max(1 quantum, max value residual)` = 0.163 K or 1.66e-6.
- `gap` applies to dashed and dotted curves only. It is the error of bridging a gap linearly, measured by cutting the same dash pattern (at 12 phases) out of the solid MATLAB curve of that surface and comparing with the full curve.
- `spread` is the height of a vertical run of the drawn path at that time.
- At t = 0 the start point is pinned to the axis origin, so `u = u_v` there.

Results on the 10 s grid:

| | Fig 2.2 MATLAB | Fig 2.2 UEL / UMATHT | Fig 2.3 MATLAB | Fig 2.3 UEL / UMATHT |
|---|---|---|---|---|
| median u | 0.21-0.24 K | 0.35-0.41 K | 4.1e-6 | 4.0e-6 |
| max u, t > 100 s | 0.61-0.73 K | 0.96-1.16 K | 1.4e-5 (at the steep rise near 160 s) | 2.4e-5 / 2.3e-5 |
| max u, all t | 4.4 K (exposed at 10 s, slope about 5 K/s) | 3.4 K (exposed at 10 s) | 1.4e-5 | 2.5e-5 (onset kink, t = 100 s) |
| gap-bridging part (max) | — | 0.49-0.65 K | — | 2.2e-5 / 2.1e-5 |

### Independent raster check (`work/raster_check.py`, 600 dpi, pixels only)

- **Ticks.** The ticks found in the raster sit within 0.02 px of the vector ticks, which confirms px = pt·600/72 − crop.
- **Raster-only calibration.** The residuals are x 0.74 s (Fig 2.2) and 0.61 s (Fig 2.3), and y ≤ 0.13 K or 1.65e-6, i.e. at most half a pixel. It differs from the vector calibration by at most 0.012 s and 0.005 K (Fig 2.2), and 0.014 s and 5e-8 (Fig 2.3).
- **Vector points on ink.**
  - Every vector point lies on ink of its own color: median 0.4-0.5 px, p99 0.7 px.
  - Exception: red points at t < 12 s lie under the blue line, which Origin drew later; the maximum there is 7.8 px.
- **No missed or extra curves.** Every colored ink pixel inside the plot is within 1.9 px of a drawn vector stroke of that color. 0 % are farther than 2 px, in both figures and all three colors.
- **Black band check.** In Fig 2.2, the top black run in each column, traced from the raster alone, matches the vector MATLAB-exposed curve with median 0.003 px, p95 0.51 px (0.17 K) and max 1.5 px (0.49 K).

## Overlap of the three codes (from the vector data, 10 s grid)

Fig 2.2, difference in K (median, with the range in brackets). One line width (0.3 pt) is 0.81 K, or 2.5 px at 600 dpi.

| surface | MATLAB − UEL | MATLAB − UMATHT | UEL − UMATHT |
|---|---|---|---|
| exposed | +3.97 [0, +4.75] | +3.21 [0, +4.30] | −0.75 [−1.59, +1.13] |
| middle | +2.61 [0, +5.26] | +0.71 [−0.21, +3.24] | −1.28 [−4.34, +1.13] |
| unexposed | +1.57 [−0.12, +5.10] | −0.14 [−1.11, +1.97] | −0.96 [−5.77, +0.52] |

- Three pairs stay within one line width (0.81 K, 2.5 px at 600 dpi) on 60 % of the grid: MATLAB/UMATHT middle, MATLAB/UMATHT unexposed, and UEL/UMATHT exposed. Their largest separations are 9.9, 6.0 and 4.9 px.
- They are nevertheless separated exactly here, by line style in the vector data.
- MATLAB exposed is consistently 3.0-4.75 K (1.1-1.75 pt) above UEL and UMATHT after 200 s. This separation is visible in the figure.

Fig 2.3, difference in kg/(s m²):

| pair | median | range |
|---|---|---|
| MATLAB − UEL | +1.7e-5 | [−1.2e-4, +2.4e-4] |
| MATLAB − UMATHT | +9.4e-6 | [−3.2e-4, +1.6e-4] |
| UEL − UMATHT | −8.3e-6 | [−5.2e-4, +2.0e-4] |

The largest differences are at the onset: |d| > 1e-4 over 130-460 s, with the maximum at t = 180-200 s. The next largest are around 660-910 s and in the decay after about 1370 s, where UMATHT is the lowest curve at 70 of the 71 grid times. All three curves are separable.

## Caveats (read before comparing)

1. **The MATLAB curves in Fig 2.2 are Origin's drawn line, not samples.**
   - The path mixes cubic Bezier and straight segments, with irregular vertex spacing.
   - Some spans are a single cubic with no vertex inside: exposed 0-56.7 s and 56.7-112.7 s; middle 0-33.9, 33.9-68.5, 68.5-136.2 and 136.2-272.4 s; unexposed 176.0-264.3 s.
   - The CSV evaluates the drawn cubic exactly. A straight chord between the same vertices would differ by up to 67.8 K (exposed, t = 20 s), 3.7 K (middle, t = 200 s) and 0.85 K (unexposed).
   - The drawn cubic follows the densely drawn UEL/UMATHT exposed curves within 1.4-3.4 K over 10-50 s, the same sign and size as their late-time offset. It is therefore a faithful curve, but the MATLAB time step cannot be recovered.
   - Three cubics of MATLAB unexposed have control points beyond their end points, so the drawn line steps back in time by ≤ 0.10 s. Handled by a cumulative maximum.
   - MATLAB middle dips to 299.92 K during 0-12 s. This is a spline overshoot below the 0.16 K quantum.
2. **Fig 2.3 MATLAB data have 4-decimal precision.**
   - All 42 vertices are integer multiples of 1e-4 kg/(s m²), within ±2.4e-6. The staircase at the peak is a property of the source data, not of the extraction.
   - Peak: **0.0046, a plateau from t = 880.5 to 909.7 s**. The curve ramps up from 0.0045 at 870.1 s and back down to 0.0045 at 920.1 s.
   - **There is no vertex between 300.6 s (0.0021) and 650.0 s (0.0035).** That span is one straight segment (a cubic with its control points at 1/3 and 2/3 of the chord), so values there are a straight-line connection, not data.
   - Vertices, as t (s) : k with value = k·1e-4:

     `0.0:0, 100.0:0, 150.0:5, 170.1:9, 180.5:11, 190.9:12, 201.2:13, 210.4:14, 251.2:18, 300.6:21, 650.0:35, 700.0:38, 750.0:40, 800.0:43, 850.0:45, 870.1:45, 880.5:46, 909.7:46, 920.1:45, 959.7:45, 970.1:44, 1000.6:43, 1050.6:42, 1100.6:40, 1150.6:39, 1200.6:37, 1250.6:36, 1300.0:34, 1350.0:33, 1400.0:32, 1450.0:30, 1500.0:29, 1550.0:27, 1600.0:26, 1649.4:24, 1699.4:23, 1749.4:21, 1799.4:19, 1850.6:18, 1900.6:16, 1950.0:15, 2000.0:13`

   - The times carry the ±0.6 s quantum.
   - The UEL and UMATHT vertices are **not** on this lattice: their median offset is 0.2 of a step.
3. **Peaks of the other codes.** UEL reaches 0.0044893, flat within one quantum over 888-918 s. UMATHT reaches 0.0045835 at 893-896 s.
4. **Exposed curves at t = 0.** UEL rises 300 → 311.7 K and UMATHT 300 → 310.9 K inside the first device column (< 0.74 s). The CSV gives 300 K at t = 0, the start point, which is also T_ini. Values at t = 10-50 s on the exposed surface carry the slope term: up to 4.4 K at 10 s.
5. **Dashed and dotted curves are bridged across gaps.** Inside a dash gap the value is linear between the drawn ends. The error is included in `_unc`, and the vertices CSV holds only what is actually inked.
6. **Grid cells after a curve's last drawn point are empty.**
   - Fig 2.2, row 2000 s: only MATLAB exposed and MATLAB unexposed have values.
   - Fig 2.2, row 1990 s: UMATHT middle and UMATHT unexposed are empty.
   - Fig 2.3: UEL is empty at 1990 and 2000 s.
   - The last drawn points are listed below.

## Spot values

All values are [S] from `zhang_fig2_2_T.csv` and `zhang_fig2_3_mloss.csv`. At these times u ≈ 0.2-0.4 K for MATLAB and ≤ 1.2 K for UEL/UMATHT (exact bounds in `_unc`).

Fig 2.2 temperature (K):

| t (s) | MATLAB exp / mid / unexp | UEL exp / mid / unexp | UMATHT exp / mid / unexp |
|---|---|---|---|
| 500 | 680.11 / 544.11 / 496.68 | 675.98 / 541.42 / 494.65 | 676.79 / 541.96 / 494.94 |
| 1000 | 717.21 / 616.72 / 593.58 | 713.14 / 615.03 / 592.70 | 714.02 / 616.21 / 593.68 |
| 2000 | 749.33 / — / 668.16 | — | — |

Last drawn points:

| curve | t (s) | T (K) |
|---|---|---|
| MATLAB middle | 1997.8 | 692.44 |
| UEL exposed | 1999.2 | 744.75 |
| UEL middle | 1999.2 | 687.55 |
| UEL unexposed | 1999.2 | 662.78 |
| UMATHT exposed | 1999.2 | 746.23 |
| UMATHT middle | 1986.7 | 690.65 |
| UMATHT unexposed | 1987.5 | 667.18 |

Fig 2.3 mass loss (kg/(s m²)):

- MATLAB peak: 0.0046, on the plateau from 880.5 to 909.7 s.
- MATLAB: 0 up to t = 100 s; 0.0013 at 2000 s.
- Last drawn points: UEL 0.0014332 at 1985.4 s; UMATHT 0.0012034 at 1998.8 s.
