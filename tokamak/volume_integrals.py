#!/usr/bin/env python
#
# Calculate and print volume integrated quantities from Hermes-2 simulation

import argparse

parser = argparse.ArgumentParser(
    description="Calculate volume integrated quantities from Hermes-2 simulation"
)
parser.add_argument("gridfile", type=str, help="The grid file used")
parser.add_argument("datapath", type=str, help="The path to the data (BOUT.dmp files)")

args = parser.parse_args()
gridfile = args.gridfile
path = args.datapath

from boutdata import collect
from boututils.datafile import DataFile
import numpy as np

t_arr = collect("t_array", path=path)
tind = len(t_arr) - 2

Nnorm = collect("Nnorm", path=path)
Tnorm = collect("Tnorm", path=path)
wci = collect("Omega_ci", path=path)
rho_s = collect("rho_s0", path=path)

power_norm = Nnorm * Tnorm * 1.602e-19 * rho_s ** 3 * wci  # Power [W]
pdens_norm = Nnorm * Tnorm * 1.602e-19 * wci  # Power density [W/m^3]

J = collect("J", path=path)[2:-2, :]  # Note: Remove X guard cells
dx = collect("dx", path=path)[2:-2, :]
dy = collect("dy", path=path)[2:-2, :]
dz = 2 * np.pi

# Volume element
dV = np.abs(J * dx * dy * dz)

# Create a mask which is 1 in the core
# Excludes X guard cells
def core_mask(grid):
    ixseps = min([grid["ixseps1"], grid["ixseps2"]])  # Innermost separatrix
    j11 = grid["jyseps1_1"]
    j12 = grid["jyseps1_2"]
    j21 = grid["jyseps2_1"]
    j22 = grid["jyseps2_2"]

    mask = np.zeros((grid["nx"] - 4, grid["ny"]))
    # Inner core
    mask[: ixseps - 2, (j11 + 1) : (j21 + 1)] = 1.0
    # Outer core
    mask[: ixseps - 2, (j12 + 1) : (j22 + 1)] = 1.0
    return mask


with DataFile(gridfile) as grid:
    ixseps = grid["ixseps1"]
    mask = core_mask(grid)

Spe = collect("Spe", path=path, tind=tind).squeeze()[2:-2, :]
Spi = collect("Spi", path=path, tind=tind).squeeze()[2:-2, :]
Sn = collect("Sn", path=path, tind=tind).squeeze()[2:-2, :]

try:
    Rzrad = collect("Rzrad", path=path, tind=tind).squeeze()[
        2:-2, :
    ]  # Impurity radiation
except:
    Rzrad = np.zeros(Sn.shape)
Rp = collect("Rp", path=path, tind=tind).squeeze()[
    2:-2, :
]  # Hydrogenic radiation + ionization cost
Qi = collect("Qi", path=path, tind=tind).squeeze()[2:-2, :]  # Transfer ions -> neutrals
Wi = collect("Wi", path=path, tind=tind).squeeze()[
    2:-2, :
]  # Transfer electrons -> ions

Spe_tot = np.sum(Spe * dV)
Spi_tot = np.sum(Spi * dV)
Sn_tot = np.sum(Sn * dV)

# Watts
input_power = (3.0 / 2) * (Spe_tot + Spi_tot) * power_norm
input_pot = Sn_tot * Nnorm * rho_s ** 3 * 13.6 * 1.602e-19 * wci

### Core input

print("\n==== CORE ====")
print(f"Simulation time: {t_arr[-1]/wci * 1e3} ms")

print("\n==== CORE ====")
print(
    f"Input power : {input_power*1e-6:.2f} MW (heat) + {input_pot*1e-6:.2f} MW (potential energy)"
)
print(f"Peak input power density : {np.max(Spe + Spi)*pdens_norm * 1e-6:.2f} MW/m^3")
print(f"Electron heating: {(3./2)*(Spe_tot) * power_norm * 1e-6:.2f} MW")
print(f"Peak electron input power density : {np.max(Spe)*pdens_norm * 1e-6:.2f} MW/m^3")

print(
    f"Core electron->ion transfer: {np.sum(Wi * dV * mask) * power_norm * 1e-6:.2f} MW"
)
print(
    f"Max core electron->ion transfer: {np.max(Wi * mask) * pdens_norm * 1e-6:.2f} MW/m^3"
)
print(
    f"Min core electron->ion transfer: {np.min(Wi * mask) * pdens_norm * 1e-6:.2f} MW/m^3"
)

zrad_core = np.sum(Rzrad * dV * mask) * power_norm
print(f"Core radiated impurity power: {zrad_core * 1e-6:.2f} MW")
print(f"Peak core impurity power: {np.max(Rzrad * mask)*pdens_norm * 1e-6:.2f} MW/m^3")

hyd_core = np.sum(Rp * dV * mask) * power_norm
print(f"Core radiated hydrogen power: {hyd_core * 1e-6:.2f} MW")
print(f"Peak core hydrogen power: {np.max(Rp * mask)*pdens_norm * 1e-6:.2f} MW/m^3")

print("\n==== SOL ====")

sol_mask = 1 - mask

print(
    f"SOL radiated impurity power: {np.sum(Rzrad * dV * sol_mask) * power_norm * 1e-6:.2f} MW"
)
print(
    f"SOL peak impurity power: {np.max(Rzrad * sol_mask) * pdens_norm * 1e-6:.2f} MW/m^3"
)
print(
    f"SOL radiated hydrogen power: {np.sum(Rp * dV * sol_mask) * power_norm * 1e-6:.2f} MW"
)
print(
    f"SOL peak hydrogen power: {np.max(Rp * sol_mask) * pdens_norm * 1e-6:.2f} MW/m^3"
)

print(
    f"SOL transfer ions->neut: {np.sum(Qi * dV * sol_mask) * power_norm * 1e-6:.2f} MW"
)
print(
    f"SOL max ion->neutral transfer power: {np.max(Qi * sol_mask) * pdens_norm * 1e-6:.2f} MW/m^3"
)
print(
    f"SOL min ion->neutral transfer power: {np.min(Qi * sol_mask) * pdens_norm * 1e-6:.2f} MW/m^3"
)
