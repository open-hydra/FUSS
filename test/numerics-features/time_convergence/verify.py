import os
import re
import matplotlib as mpl
import matplotlib.pyplot as plt

# Match aesthetic of 1D-transient-oscillating-rod-probes.svg
mpl.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        9,
    "axes.titlesize":   10,
    "axes.titleweight": "bold",
    "axes.labelsize":   11,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "legend.fontsize":  9,
})

SOLUTION_DIR = os.path.join(os.path.dirname(__file__), "SOLUTION")

_mpi_re = re.compile(r"MPI-(\d+)x(\d+)$")
_omp_re = re.compile(r"openMP-(\d+)$")
_time_re = re.compile(r"Time of operation was\s+([\d.]+(?:[Ee][+-]?\d+)?)")

# --- read all files in SOLUTION/ ---------------------------------------------
mpi_points = []   # (cores, time, label)
omp_points = []   # (cores, time)

for filename in sorted(os.listdir(SOLUTION_DIR)):
    filepath = os.path.join(SOLUTION_DIR, filename)
    if not os.path.isfile(filepath):
        continue

    # extract timing from the last matching line in the file
    time_val = None
    with open(filepath) as f:
        for line in f:
            m = _time_re.search(line)
            if m:
                time_val = float(m.group(1))
    if time_val is None:
        print(f"WARNING: no timing line found in {filename}, skipping.")
        continue

    m = _mpi_re.match(filename)
    if m:
        ranks, threads = int(m.group(1)), int(m.group(2))
        mpi_points.append((ranks * threads, time_val, f"{ranks}x{threads}", ranks))
        continue

    m = _omp_re.match(filename)
    if m:
        omp_points.append((int(m.group(1)), time_val))
        continue

    print(f"WARNING: could not parse case type from filename '{filename}', skipping.")

mpi_points.sort(key=lambda p: (p[3], p[0]))
omp_points.sort(key=lambda p: p[0])

# --- plot --------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(568.8 / 72, 280.8 / 72))

# openMP: blue squares + dashed blue line
if omp_points:
    omp_x, omp_t = zip(*omp_points)
    ax.plot(omp_x, omp_t, color="tab:blue", linestyle="--", linewidth=1,
            zorder=1, label="_nolegend_")
    ax.scatter(omp_x, omp_t, color="tab:blue", marker="s", s=60, zorder=2,
               label="openMP")

# MPI: red markers (distinct per case) + one dashed red line per rank group
if mpi_points:
    # group by rank count
    rank_groups = {}
    for point in mpi_points:
        rank_groups.setdefault(point[3], []).append(point)

    group_markers = ["o", "^", "D", "v", "P", "*"]
    for i, ranks in enumerate(sorted(rank_groups)):
        group  = rank_groups[ranks]
        marker = group_markers[i % len(group_markers)]
        gx     = [p[0] for p in group]
        gt     = [p[1] for p in group]
        ax.plot(gx, gt, color="tab:red", linestyle="--", linewidth=1, zorder=1,
                label="_nolegend_")
        for j, (cores, t, lbl, _) in enumerate(group):
            ax.scatter(cores, t, color="tab:red", marker=marker, s=60, zorder=2,
                       label=f"{ranks} x n" if j == 0 else "_nolegend_")

all_cores = sorted(set(
    ([p[0] for p in mpi_points] if mpi_points else []) +
    ([p[0] for p in omp_points] if omp_points else [])
))
ax.set_xticks(all_cores)
ax.set_xticklabels([str(c) for c in all_cores])

ax.set_xlabel("Number of cores")
ax.set_ylabel("Time (min)")
ax.set_title("Execution Time Comparison")
ax.legend()
ax.grid(True, linestyle="--", alpha=0.5)

plt.tight_layout()
output_path = os.path.join(os.path.dirname(__file__), "time_convergence.svg")
plt.savefig(output_path)
print(f"Figure saved to {output_path}")
plt.show()
