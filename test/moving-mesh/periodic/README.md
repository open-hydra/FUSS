# Periodic faces (201) with the moving mesh and type-2 pyrolysis (plan 14 S9d, Part B)

One sector of upstream's `PeriodicSector` (`test/numerics-features/periodic_sector/make_sector.py`: a helical annular wall
about the x axis, 6 sectors in the ring, 6 x 6 x 4 cells per sector) with ours' features switched on. `generate.py`
imports `make_sector.py` (it edits nothing there) and writes:

| deck | variant | what moves or flows across the pair |
|---|---|---|
| `translate/` | `translate` | `law = translation` along x, the axis of the pair |
| `prescribed/` | `prescribed` | `law = prescribed`, a displacement along x that depends on x only |
| `receding/`, `receding1/` | `receding` | `law = nodeshift`: the inner wall is a 306 receding at a uniform rate (`sdot-override`); two blocks, and one block periodic with itself |
| `../../pyrolysis/darcy/periodic_static/`, `periodic_static1/` | `darcy` | `gas-closure = type2` on a static mesh, `K` transversely isotropic about x; two blocks, and one block periodic with itself |
| `../../pyrolysis/ablation/periodic_receding/` | `ablation` | type 2 and `nodeshift`: the inner wall a 306 with the TACOT B' table, `env_hi.dat` on the first third of every sector and `env_lo.dat` on the rest (w5_2blk's env_hi with rhoeUeCH0 × 0.1 and × 0.05: the rate differs 2× across the seam, the recession stays ~11 % of the line budget), `smoothing-width = 5`; 20 steps |

Every motion here takes a node of one periodic face and its image to images: that is what a periodic domain can
follow (`src/lib/numerics/space/Mod_Periodic_Checks.f90`). The ring of 6 sectors joined by block connections (101) is
the oracle: `run_physics_gates.sh` (gates `periodic_motion`, `periodic_receding`, `periodic_darcy`) writes it with
`generate.py <variant> ring <dir>`, runs it beside the committed sector and compares with `compare.py` (T, and the nodes
of a moving mesh, over the six sectors of the ring; type 2 also `p_pg`); the sector with symmetry faces (300) instead of
the pair must differ from the ring. The guard demos (`run_guard_demos.sh`) run the committed decks (edited for
`across`, `wobble`, `aniso`) and write one base with `generate.py` (`openseam`).

Regenerate from the repository root: `python3 test/moving-mesh/periodic/generate.py --decks`. Do not edit the decks: the
physics gate regenerates each sector and stops when the committed one differs.
