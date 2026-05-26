import numpy as np

f = 0.0125          # Hz
T = 1.0 / f         # period = 80 s
dt = 1.0            # time step [s]
t = np.arange(0.0, T, dt)
values = 2000.0 + 1000.0 * np.sin(2.0 * np.pi * f * t)

with open("Tw_time.dat", "w") as fout:
    fout.write("periodic\n")
    for ti, vi in zip(t, values):
        fout.write(f"{ti:.6e}  {vi:.6e}\n")
