import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Match aesthetic of 1D-transient-oscillating-rod-probes.svg
# Figure: 568.8 x 280.8 pt  |  tick labels: 9pt  |  axis labels: 11pt  |  title: 10pt bold
mpl.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         9,          # tick labels
    "axes.titlesize":    10,
    "axes.titleweight":  "bold",
    "axes.labelsize":    11,
    "xtick.labelsize":   9,
    "ytick.labelsize":   9,
    "legend.fontsize":   9,
})

SOLUTION_DIR = os.path.join(os.path.dirname(__file__), "SOLUTION")

CASES = {
    "Nominal":                          "nominal.dat",
    "2 level MG (25000 + conv.)":       "mg1.dat",
    "2 level MG (50000 + conv.)":       "mg2.dat",
    "2 level MG (75000 + conv.)":       "mg3.dat",
    "IRS":                              "irs.dat",
    "2 level MG (conv. + conv.) + IRS": "mg_irs.dat"
}

# 568.8 x 280.8 pt converted to inches (1 pt = 1/72 in)
fig, ax = plt.subplots(figsize=(568.8 / 72, 280.8 / 72))

for label, filename in CASES.items():
    filepath = os.path.join(SOLUTION_DIR, filename)
    if not os.path.isfile(filepath):
        print(f"WARNING: {filepath} not found, skipping.")
        continue
    data = np.loadtxt(filepath)
    iterations = data[:, 0]
    residuals  = data[:, 2]
    ax.semilogy(iterations, residuals, label=label)

ax.set_xlabel("Iteration")
ax.set_ylabel("Residual")
ax.set_title("Residual History Comparison")
ax.legend()
ax.grid(True, which="both", linestyle="--", alpha=0.5)

plt.tight_layout()
output_path = os.path.join(os.path.dirname(__file__), "residual_convergence.svg")
plt.savefig(output_path)
print(f"Figure saved to {output_path}")
plt.show()
