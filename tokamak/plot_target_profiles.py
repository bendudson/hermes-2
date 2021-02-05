#!/usr/bin/env python
#
# Make plots of the density, temperature and heat flux on the inner and outer targets
# Saves to outer_target_profiles and outer_target_profiles PDF and PNG files

import argparse

parser = argparse.ArgumentParser(description="Make plots of target profiles")
parser.add_argument("gridfile", type=str, help="The grid file used")
parser.add_argument("datapath", type=str, help="The path to the data (BOUT.dmp files)")

args = parser.parse_args()
gridfile = args.gridfile
path = args.datapath

# sheath heat transmission
gamma_e = 4.0
gamma_i = 2.5

from boutdata import collect

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from boututils.datafile import DataFile
import numpy as np


with DataFile(gridfile) as grid:
    Rxy = grid["Rxy"]
    Zxy = grid["Zxy"]
    Bpxy = grid["Bpxy"]
    Bxy = grid["Bxy"]
    ixsep = grid["ixseps1"]

J = collect("J", path=path)
g_22 = collect("g_22", path=path)
dx = collect("dx", path=path)
dy = collect("dy", path=path)
dz = 2 * np.pi

dV = abs(J * dx * dy * dz)  # Volume element
dS = dV / (dy * np.sqrt(g_22))  # Area element

# Calculate distance along target

# Outer target
dR = Rxy[3:-1, -1] - Rxy[2:-2, -1]
dZ = Zxy[3:-1, -1] - Zxy[2:-2, -1]

dl = np.sqrt(dR ** 2 + dZ ** 2)
distance = np.cumsum(dl)

# Inner target
dR = Rxy[3:-1, 0] - Rxy[2:-2, 0]
dZ = Zxy[3:-1, 0] - Zxy[2:-2, 0]

dl_inner = np.sqrt(dR ** 2 + dZ ** 2)
distance_inner = np.cumsum(dl_inner)

# Set to zero at separatrix, which is half-way between
# grid centers
distance -= 0.5 * (distance[ixsep] + distance[ixsep - 1])
distance_inner -= 0.5 * (distance_inner[ixsep] + distance_inner[ixsep - 1])

t_arr = collect("t_array", path=path)
tind = len(t_arr) - 2

Nnorm = collect("Nnorm", path=path)
Tnorm = collect("Tnorm", path=path)
wci = collect("Omega_ci", path=path)
cs0 = collect("Cs0", path=path)
rho_s = collect("rho_s0", path=path)

ni = collect("Ne", path=path, tind=tind, yguards=True).squeeze() * Nnorm  # m^-3
nn = collect("Nn", path=path, tind=tind, yguards=True).squeeze() * Nnorm  # m^-3
te = collect("Telim", path=path, tind=tind, yguards=True).squeeze() * Tnorm  # eV
ti = collect("Ti", path=path, tind=tind, yguards=True).squeeze() * Tnorm  # eV
pn = collect("Pn", path=path, tind=tind, yguards=True).squeeze()
tn = pn / (nn / Nnorm) * Tnorm  # eV

# Outer target
ni_targ = 0.5 * (ni[:, -3] + ni[:, -2])[2:-2]
nn_targ = 0.5 * (nn[:, -3] + nn[:, -2])[2:-2]
te_targ = 0.5 * (te[:, -3] + te[:, -2])[2:-2]
ti_targ = 0.5 * (ti[:, -3] + ti[:, -2])[2:-2]
tn_targ = 0.5 * (tn[:, -3] + tn[:, -2])[2:-2]

cs = np.sqrt((ti_targ + te_targ) / Tnorm) * cs0  # m/s

q_heat_e = gamma_e * 1.602e-19 * te_targ * ni_targ * cs  # W/m^2
q_heat_i = gamma_i * 1.602e-19 * ti_targ * ni_targ * cs  # W/m^2
q_pot = ni_targ * cs * 1.602e-19 * (13.6 + 2.2)  # Potential energy

# Inner target

ni_targ_inner = 0.5 * (ni[:, 1] + ni[:, 2])[2:-2]
nn_targ_inner = 0.5 * (nn[:, 1] + nn[:, 2])[2:-2]
te_targ_inner = 0.5 * (te[:, 1] + te[:, 2])[2:-2]
ti_targ_inner = 0.5 * (ti[:, 1] + ti[:, 2])[2:-2]
tn_targ_inner = 0.5 * (tn[:, 1] + tn[:, 2])[2:-2]

cs_inner = np.sqrt((ti_targ_inner + te_targ_inner) / Tnorm) * cs0  # m/s

q_heat_e_inner = gamma_e * 1.602e-19 * te_targ_inner * ni_targ_inner * cs_inner  # W/m^2
q_heat_i_inner = gamma_i * 1.602e-19 * ti_targ_inner * ni_targ_inner * cs_inner  # W/m^2
q_pot_inner = ni_targ_inner * cs_inner * 1.602e-19 * (13.6 + 2.2)  # Potential energy

###########

# Convert parallel power flux to perpendicular
qperp = (q_heat_e + q_heat_i + q_pot) * np.abs(Bpxy[2:-2, -1] / Bxy[2:-2, -1])

# Integrate to get total power
power_total = np.trapz(qperp * dl * 2 * np.pi * Rxy[2:-2, -1])
power_total2 = np.sum(
    (q_heat_e + q_heat_i + q_pot) * dS[2:-2, -1] * rho_s ** 2
)  # Cross-check

print("Total power outer target power: {} MW".format(power_total * 1e-6))
print("Total power outer target power (alt.): {} MW".format(power_total2 * 1e-6))


qperp = (q_heat_e_inner + q_heat_i_inner + q_pot_inner) * np.abs(
    Bpxy[2:-2, 0] / Bxy[2:-2, 0]
)
power_total = np.trapz(qperp * dl_inner * 2 * np.pi * Rxy[2:-2, 0])
print("Total power inner target power: {} MW".format(power_total * 1e-6))

############
# Outer target

fig, ax = plt.subplots(1, 3, figsize=(10, 4))

ax[0].plot(distance, ni_targ, "-k", label="Ion")
ax[0].plot(distance, nn_targ, "--b", label="Atom")
ax[0].set_yscale("log")
ax[0].set_xlabel("Perpendicular distance [m]")
ax[0].axvline(0.0, linestyle="--", color="k")
ax[0].legend()
ax[0].set_ylabel(r"Density [m$^{-3}$]")

ax[1].plot(distance, ti_targ, "-k", label="Ion")
ax[1].plot(distance, te_targ, "-r", label="Electron")
ax[1].plot(distance, tn_targ, "--b", label="Atom")
ax[1].set_ylabel(r"Temperature [eV]")
ax[1].set_xlabel("Perpendicular distance [m]")
ax[1].axvline(0.0, linestyle="--", color="k")
ax[1].set_ylim(bottom=0.0)
ax[1].legend()

ax[2].plot(distance, q_heat_i * 1e-6, "-k", label="Ion")
ax[2].plot(distance, q_heat_e * 1e-6, "-r", label="Electron")
ax[2].plot(distance, q_pot * 1e-6, "-g", label="Potential")
ax[2].set_ylabel(r"Parallel heat flux [MW/m$^2$]")
ax[2].set_xlabel("Perpendicular distance [m]")
ax[2].set_ylim(bottom=0.0)
ax[2].axvline(0.0, linestyle="--", color="k")
ax[2].legend()

fig.tight_layout()

plt.savefig("outer_target_profiles.pdf")
plt.savefig("outer_target_profiles.png")
plt.close("all")

############
# Inner target

fig, ax = plt.subplots(1, 3, figsize=(10, 4))

ax[0].plot(distance_inner, ni_targ_inner, "-k", label="Ion")
ax[0].plot(distance_inner, nn_targ_inner, "--b", label="Atom")
ax[0].set_yscale("log")
ax[0].set_xlabel("Perpendicular distance [m]")
ax[0].axvline(0.0, linestyle="--", color="k")
ax[0].legend()
ax[0].set_ylabel(r"Density [m$^{-3}$]")

ax[1].plot(distance_inner, ti_targ_inner, "-k", label="Ion")
ax[1].plot(distance_inner, te_targ_inner, "-r", label="Electron")
ax[1].plot(distance_inner, tn_targ_inner, "--b", label="Atom")
ax[1].set_ylabel(r"Temperature [eV]")
ax[1].set_xlabel("Perpendicular distance [m]")
ax[1].axvline(0.0, linestyle="--", color="k")
ax[1].set_ylim(bottom=0.0)
ax[1].legend()

ax[2].plot(distance_inner, q_heat_i_inner * 1e-6, "-k", label="Ion")
ax[2].plot(distance_inner, q_heat_e_inner * 1e-6, "-r", label="Electron")
ax[2].plot(distance_inner, q_pot_inner * 1e-6, "-g", label="Potential")
ax[2].set_ylabel(r"Parallel heat flux [MW/m$^2$]")
ax[2].set_xlabel("Perpendicular distance [m]")
ax[2].set_ylim(bottom=0.0)
ax[2].axvline(0.0, linestyle="--", color="k")
ax[2].legend()

xticks = [-0.05, 0.0, 0.05]
for i in range(3):
    ax[i].set_xticks(xticks)


fig.tight_layout()

plt.savefig("inner_target_profiles.pdf")
plt.savefig("inner_target_profiles.png")
plt.close("all")
