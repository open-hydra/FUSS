"""plot_convergence.py <runs dir> [<out prefix>]

Residual histories of the runs of run_configs.sh, in two figures on linear x axes:
<prefix>_iterations.png against the finest-level iterations and <prefix>_time.png against the
wall time (the Elapsed line of each run's FUSS timing summary, spread uniformly over its
iterations). One column per case; the top row shows the full range, the bottom row zooms on the
range of the four FAS cycles (the anisotropic plate, where the cycles end with the single-grid
runs, keeps the full range). Colour = cycle, dashed = IRS; a run that did not end with status 0
is marked with a cross where it stopped. The default prefix is <runs dir>/convergence. A summary
(iterations, wall time, last residual, status) is printed for every run found.
"""
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402

CASES = [('plate', 'Two-material plate\n200 × 10, 2 blocks (101), 2 levels'),
         ('liner', 'Four-layer liner wall\n64 × 128, 2 blocks (101), jumps inside, 4 levels'),
         ('chimera_mg', 'Overlapping chimera\n2 blocks (102), 3 levels'),
         ('anisotropic', 'Anisotropic plate\n128 × 16, aspect ratio 100, 3 levels')]
# colour = cycle (three slots of a palette validated for colour-vision deficiency), dash = IRS
CYCLE = {'none': ('#2a78d6', 'none'), 'v': ('#eb6834', 'V'), 'w': ('#1baf7a', 'W')}
CFGS = [('none', False), ('none', True), ('v', False), ('v', True), ('w', False), ('w', True)]
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'


def label(cyc, irs):
    return ('IRS' if irs else 'none') if cyc == 'none' else CYCLE[cyc][1] + (' + IRS' if irs else '')


def run_name(case, cyc, irs):
    return f'{case}_' + (('irs' if irs else 'none') if cyc == 'none' else cyc + ('_irs' if irs else ''))


def load(run):
    """iterations, residuals, times, ok, elapsed, threshold of one run; None without a history."""
    hist = os.path.join(run, 'OUTPUT', 'residual-history.dat')
    if not os.path.isfile(hist):
        return None
    it, r = [], []
    for line in open(hist):
        p = line.split()
        if len(p) == 3 and p[0].isdigit():
            it.append(int(p[0])); r.append(float(p[2]))
    if not it:
        return None
    log = open(os.path.join(run, 'run.log'), errors='replace').read()
    m = re.search(r'Elapsed\s+([0-9.E+-]+)', log)
    status = open(os.path.join(run, 'run.time')).read().split() if os.path.isfile(os.path.join(run, 'run.time')) \
        else ['exit', '?', '0']
    elapsed = float(m.group(1)) if m else float(status[2])
    ini = open(os.path.join(run, 'input.ini')).read()
    thr = float(re.search(r'^res-threshold\s*=\s*(\S+)', ini, re.M).group(1))
    it, r = np.array(it), np.array(r)
    return it, r, elapsed * it / it[-1], status[1] == '0', elapsed, thr


def main(runs, prefix):
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK2,
                         'xtick.color': INK2, 'ytick.color': INK2, 'axes.titlesize': 9.5,
                         'axes.titlecolor': INK, 'figure.facecolor': SURF, 'axes.facecolor': SURF,
                         'savefig.facecolor': SURF})
    data, summary = {}, []
    for case, _ in CASES:
        for cyc, irs in CFGS:
            d = load(os.path.join(runs, run_name(case, cyc, irs)))
            if d is None:
                continue
            data[(case, cyc, irs)] = d
            it, r, _, ok, el, _ = d
            summary.append(f'{run_name(case, cyc, irs):22s} iterations {it[-1]:7d}  wall {el:8.1f} s  '
                           f'residual {r[-1]:.2e}  {"ok" if ok else "ERROR"}')
    if not data:
        sys.exit(f'plot_convergence.py: no runs of run_configs.sh under {runs}')
    thr = min(d[5] for d in data.values())

    for what, xlabel, unit in (('iterations', 'finest-level iterations', ''), ('time', 'wall time [s]', ' s')):
        fig, axes = plt.subplots(2, len(CASES), figsize=(15, 7.6), sharey=True)
        for c, (case, title) in enumerate(CASES):
            ends = {}
            for cyc, irs in CFGS:
                if (case, cyc, irs) not in data:
                    continue
                it, r, t, ok, _, _ = data[(case, cyc, irs)]
                x = it if what == 'iterations' else t
                ends[(cyc, irs)] = x[-1]
                for row in (0, 1):
                    axes[row, c].plot(x, r, color=CYCLE[cyc][0], lw=1.6, ls=(0, (4, 2)) if irs else '-')
                    if not ok:
                        axes[row, c].plot(x[-1], r[-1], marker='x', ms=9, mew=2, color=CYCLE[cyc][0], ls='none')
            axes[0, c].set_title(title, loc='left')
            if not ends:
                continue
            full = max(ends.values())
            cycles = [v for (cyc, _), v in ends.items() if cyc != 'none']
            zoom = 1.12 * max(cycles) if cycles and case != 'anisotropic' else 1.12 * full
            axes[0, c].set_xlim(0, 1.03 * full)
            axes[1, c].set_xlim(0, zoom)
            txt = '\n'.join(f'{label(cyc, irs):9s} {v:>9,.0f}{unit}' if what == 'iterations' or v >= 10
                            else f'{label(cyc, irs):9s} {v:>9.1f}{unit}' for (cyc, irs), v in ends.items())
            axes[0, c].text(0.98, 0.97, txt, transform=axes[0, c].transAxes, ha='right', va='top',
                            family='monospace', fontsize=7.5, color=INK2)
        for ax in axes.flat:
            ax.set_yscale('log')
            ax.grid(True, which='major', color=GRID, lw=0.6)
            ax.axhline(thr, color=INK2, lw=0.8, ls=':')
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))
            if ax.get_xlim()[1] >= 10000:
                fmt = lambda v, _: '0' if v == 0 else f'{v / 1000:g}k'  # noqa: E731
            else:
                fmt = lambda v, _: f'{v:,.0f}' if abs(v) >= 10 or v == 0 else f'{v:g}'  # noqa: E731
            ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt))
            for s in ('top', 'right'):
                ax.spines[s].set_visible(False)
        for c, (case, _) in enumerate(CASES):
            axes[0, c].set_xlabel(xlabel)
            axes[1, c].set_xlabel(xlabel + (', zoom on the FAS cycles' if case != 'anisotropic' else ''))
        for row in (0, 1):
            axes[row, 0].set_ylabel('residual  ||ΔT|| per iteration [K]')
        h = [Line2D([], [], color=CYCLE[cyc][0], lw=1.6, ls=(0, (4, 2)) if irs else '-') for cyc, irs in CFGS]
        l = [label(cyc, irs) for cyc, irs in CFGS]
        h.append(Line2D([], [], color=INK2, lw=0.8, ls=':')); l.append(f'threshold {thr:g}')
        fig.legend(h, l, loc='upper center', ncol=7, frameon=False, bbox_to_anchor=(0.5, 1.0),
                   handlelength=3, labelcolor=INK)
        note = (f'RK3, cons, threshold {thr:g}. Without IRS vnn 1.0, with IRS vnn 2.0 and irs-beta 0.5; none and '
                'IRS on one grid, cycles on the case levels with a 1-iteration descent. Numbers top right: where '
                'each run reaches the threshold.')
        if what == 'time':
            note += ('\nWall time: the Elapsed line of each run, one core per run, several runs at a time on a '
                     'shared node, spread uniformly over the iterations of the run.')
        fig.text(0.5, 0.004, note, ha='center', color=INK2, fontsize=8)
        fig.tight_layout(rect=(0, 0.035, 1, 0.95))
        out = f'{prefix}_{what}.png'
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print(out)
    print('\n'.join(summary))


if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else os.path.join(sys.argv[1], 'convergence'))
