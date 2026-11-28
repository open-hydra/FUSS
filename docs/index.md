---
title: FUSS
hide:
  - navigation
  - toc
---

<div class="fuss-hero" markdown>

<div class="fuss-hero__text" markdown>

<p class="fuss-hero__eyebrow">HYDRA CFD SUITE</p>

<h1 class="fuss-hero__title">FUSS</h1>

<p class="fuss-hero__lead">Open-source finite-volume solver for steady and time-accurate heat conduction in solids,
on multi-block structured grids with temperature-dependent materials.</p>

<p class="fuss-hero__cmd"><span>$</span> ./install.sh build --compilers=gnu --use-openmp</p>

[Get started](getting-started/index.md){ .md-button .md-button--primary }
[Verification](vv/index.md){ .md-button .fuss-hero__ghost }

<p class="fuss-hero__name"><strong>F</strong>ourier <strong>U</strong>n-<strong>S</strong>teady <strong>S</strong>olver</p>

</div>

<!-- The hero figure is docs/images/home-hero.png, drawn by extra.css as a background: without the file the
     hero shows only its own background. -->
<div class="fuss-hero__shot" aria-hidden="true"><div class="fuss-hero__figure"></div></div>

</div>

<p class="fuss-kicker">The solver</p>

## Heat conduction in solids { .fuss-heading }

<p class="fuss-sub" markdown>FUSS solves the Fourier heat equation on the solid parts of a case, alone or coupled
to the gas inside Hydra. Start with the [Quick Start](getting-started/quick-start.md): the 2-D plate case runs in
seconds and checks itself.</p>

<div class="fuss-tools" markdown>

[:material-thermometer: <b>Solid heat conduction</b> <i>Steady or time-accurate, with a volumetric source; k, ρ and c<sub>p</sub> tabulated against temperature per material.</i>](theory/governing-equations.md){ .fuss-tool }
[:material-grid: <b>Cell-centred finite volume</b> <i>Multi-block structured grids; heat-flux continuity at block connections and between materials.</i>](theory/numerics.md){ .fuss-tool }
[:material-wall: <b>Boundary conditions</b> <i>Flux or temperature (constant or in time), convection, radiation, symmetry, periodic faces, chimera blocks and coupling to the gas.</i>](user/boundary-conditions.md){ .fuss-tool }
[:material-fire: <b>Pyrolysis and ablation</b> <i>Pyrolysing porous materials with a quasi-steady or Darcy gas closure; an ablating wall against B′ tables.</i>](user/input.md#pyrolysis){ .fuss-tool }
[:material-timer-sand: <b>Accelerated time marching</b> <i>Euler, RK2 and SSP RK3; local time step, residual smoothing and FAS multigrid for steady runs.</i>](theory/time-integration.md){ .fuss-tool }
[:material-chip: <b>Parallel runs</b> <i>OpenMP inside the blocks, MPI over whole blocks; one process and N ranks write identical files.</i>](user/parallel.md){ .fuss-tool }
[:material-file-chart: <b>Flexible I/O</b> <i>Tecplot ASCII and binary, VTK, wall temperature and heat-flux files, point probes and restart.</i>](user/using.md){ .fuss-tool }

</div>

<p class="fuss-kicker"></p>

## Hydra CFD Suite { .fuss-heading }

<p class="fuss-sub">Hydra couples solvers of gases, real fluids, solids, particles and radiation in one
simulation. FUSS is its solver of heat conduction in solids.</p>

<div class="hydra" markdown>

<div class="hydra__case" markdown>

<p class="hydra__label">Your case</p>

<p class="hydra__file" markdown>:material-grid: structured mesh</p>

<p class="hydra__file" markdown>:material-database: [databases](https://open-hydra.github.io/ATLAS/databases/)</p>

</div>

<div class="hydra__atlas" markdown>

<p class="hydra__name">ATLAS</p>

<p class="hydra__tools" markdown>[GPB](https://open-hydra.github.io/ATLAS/user-guide/gpb/) [BCB](https://open-hydra.github.io/ATLAS/user-guide/bcb/) [ICB](https://open-hydra.github.io/ATLAS/user-guide/icb/) [STB](https://open-hydra.github.io/ATLAS/user-guide/stb/) [MDB](https://open-hydra.github.io/ATLAS/user-guide/mdb/) [KAnT](https://open-hydra.github.io/ATLAS/user-guide/kant/)</p>

</div>

<div class="hydra__link">
<span class="hydra__link-line"></span>
<span class="hydra__link-bottom">input files for each run</span>
</div>

<div class="hydra__suite" markdown>

<p class="hydra__label" markdown>[Hydra](https://github.com/open-hydra) · the coupled simulation</p>

<div class="hydra__solvers" markdown>

[<b>MOSE</b> compressible flow <i>ideal-gas</i>](https://github.com/open-hydra/MOSE){ .hydra__solver .hydra__solver--gas }
[<b>GROOT</b> radiative transfer <i>no ATLAS phase</i>](https://github.com/open-hydra/GROOT){ .hydra__solver .hydra__solver--radiation }
[<b>ARES</b> real-fluid flow <i>real-fluid</i>](https://github.com/open-hydra/ARES){ .hydra__solver .hydra__solver--fluid }
[<b>FUSS</b> heat conduction in solids <i>solid</i>](https://github.com/open-hydra/FUSS){ .hydra__solver .hydra__solver--solid .is-here }
[<b>IGLOO</b> Lagrangian particles <i>dispersed</i>](https://github.com/open-hydra/IGLOO){ .hydra__solver .hydra__solver--particle }
[<b>ICE</b> Eulerian particles <i>dispersed</i>](https://github.com/open-hydra/ICE){ .hydra__solver .hydra__solver--particle }

</div>

</div>

<div class="hydra__libs" markdown>

<p class="hydra__label">Built on</p>

[<b>ORION</b> Tecplot, VTK and PLOT3D I/O](https://MarcoGrossi92.github.io/ORION/){ .hydra__lib }
[<b>FLINT</b> thermochemistry](https://MarcoGrossi92.github.io/FLINT/){ .hydra__lib }
[<b>OSLO</b> ODE integrators](https://github.com/MarcoGrossi92/OSLO){ .hydra__lib }
[<b>FiNeR</b> INI file parser](https://github.com/szaghi/FiNeR){ .hydra__lib }

</div>

</div>

<div class="fuss-links" markdown>

[Quick Start](getting-started/quick-start.md)
[User Guide](user/index.md)
[Input Parameters](user/registry.md)
[Theory](theory/index.md)
[Verification](vv/index.md)
[Developer Guide](development/index.md)

</div>
