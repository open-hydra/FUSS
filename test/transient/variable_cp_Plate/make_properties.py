#!/usr/bin/env python3
"""INPUT/properties.dat of the variable_cp_Plate test: a steel-like solid with a specific heat linear in T,
cp(T) = A + B T, constant density and conductivity. The Energy column is the volumetric enthalpy
h(T) = int_0^T rho cp dT' = rho (A T + B T^2 / 2), exact at every integer temperature (ATLAS GPB writes
constant-property tables only, so this file is written here).

  python3 make_properties.py
"""
import os

A, B = 300.0, 0.4           # cp = A + B T  [J/(kg K)]: 420 at 300 K, 828 at 1320 K
RHO, K = 7850.0, 52.0       # [kg/m3], [W/(m K)]
T_MAX = 2000                # table from 1 K to T_MAX K, 1 K apart (FUSS indexes it by integer kelvin)

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, 'INPUT', 'properties.dat'), 'w') as f:
    f.write('TITLE = "Mass Thermodynamic Properties"\n')
    f.write('VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\n')
    f.write('ZONE T="mat"\n')
    f.write(f'I={T_MAX}, F=POINT\n')
    for T in range(1, T_MAX + 1):
        f.write(f'{float(T):.1f} {A + B*T:.6f} {RHO:.6f} {K:.6f} {RHO*(A*T + 0.5*B*T*T):.6f}\n')
