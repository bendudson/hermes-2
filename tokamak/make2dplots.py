#!/usr/bin/env python
#
# Makes 2D plots:
#  - ne_nn_te_2d.png   Electron density, neutral density and electron temperature
#  - radiation2d.png   Impurity and hydrogenic radiation
#  - transfer2d.png    Ion <-> Electron and Ion <-> Neutral energy exchange
#

import argparse

parser = argparse.ArgumentParser(description="Make 2D plots of Hermes-2 data")
parser.add_argument("gridfile", type=str, help="The grid file used")
parser.add_argument("datapath", type=str, help="The path to the data (BOUT.dmp files)")

args = parser.parse_args()
gridfile = args.gridfile
path = args.datapath

from boutdata import collect
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter

from boututils.datafile import DataFile
from boutdata.griddata import gridcontourf
import numpy as np

t_arr = collect("t_array", path=path)
tind = len(t_arr) - 2

Nnorm = collect("Nnorm", path=path)
Tnorm = collect("Tnorm", path=path)
wci = collect("Omega_ci", path=path)

grid = DataFile(gridfile)

ni = collect("Ne", path=path, tind=tind).squeeze() * Nnorm  # m^-3
nn = collect("Nn", path=path, tind=tind).squeeze() * Nnorm  # m^-3
te = collect("Telim", path=path, tind=tind).squeeze() * Tnorm  # eV

#######################

fig, ax = plt.subplots(1, 3, figsize=(10, 8))

c = gridcontourf(grid, ni, ax=ax[0], show=False)
fig.colorbar(c, ax=ax[0])
ax[0].set_title(r"Ion dens [m$^{3}$]")
ax[0].axes.get_xaxis().set_ticks([1, 2, 3, 4])

c = gridcontourf(grid, nn, ax=ax[1], show=False, ylabel=None)
fig.colorbar(c, ax=ax[1])
ax[1].set_title(r"Atom dens [m$^{3}$]")
ax[1].axes.get_yaxis().set_ticks([])
ax[1].axes.get_xaxis().set_ticks([1, 2, 3, 4])

c = gridcontourf(grid, te, ax=ax[2], log=True, nlevel=20, show=False, ylabel=None)
cbar = fig.colorbar(c, ax=ax[2])
cbar.set_ticks([0.1, 0.5, 1.0, 5.0, 10.0, 50, 100, 200, 500])
ax[2].set_title("Elec. temp [eV]")
ax[2].axes.get_yaxis().set_ticks([])
ax[2].axes.get_xaxis().set_ticks([1, 2, 3, 4])

fig.tight_layout(pad=3.0)

plt.savefig("ne_nn_te_2d.png")
plt.show()
plt.close("all")

#######################

rp = collect("Rp", path=path, tind=tind).squeeze() * Nnorm * Tnorm * 1.602e-19 * wci

try:
    rzrad = (
        collect("Rzrad", path=path, tind=tind).squeeze()
        * Nnorm
        * Tnorm
        * 1.602e-19
        * wci
    )  # W/m^3
except:
    rzrad = np.zeros(rp.shape)

# Same scale for all plots
maxd = max(np.max(rzrad), np.max(rp))
mind = min(np.min(rzrad), np.min(rp))

fig, ax = plt.subplots(2, 1, figsize=(10, 7), gridspec_kw={"height_ratios": [10, 1]})
from mpl_toolkits.axes_grid1 import make_axes_locatable

ne_ax = ax[0]
cbar_ax = ax[1]
divider = make_axes_locatable(ne_ax)
h_ax = divider.append_axes("right", size="100%", pad=0.01)

c = gridcontourf(grid, rzrad, ax=ne_ax, show=False, mind=mind, maxd=maxd)
ne_ax.set_title(r"Impurity radiation [Wm$^{-3}$]")
ne_ax.axes.get_xaxis().set_ticks([1, 2, 3, 4])
ne_ax.set_ylim(bottom=4.0)

ne_ax.xaxis.grid(True)
ne_ax.yaxis.grid(True)

c = gridcontourf(grid, rp, ax=h_ax, show=False, ylabel=None, mind=mind, maxd=maxd)
h_ax.set_title(r"H radiation [Wm$^{-3}$]")
h_ax.axes.get_yaxis().set_ticklabels([])
h_ax.axes.get_xaxis().set_ticks([1, 2, 3, 4])
h_ax.set_ylim(bottom=4.0)

h_ax.xaxis.grid(True)
h_ax.yaxis.grid(True)

cbar = fig.colorbar(c, cax=cbar_ax, orientation="horizontal")

fig.tight_layout()

plt.savefig("radiation2d.png")

plt.show()
plt.close("all")

###############################
# Power transfer

qi = collect("Qi", path=path, tind=tind).squeeze() * Nnorm * Tnorm * 1.602e-19 * wci
wi = (
    collect("Wi", path=path, tind=tind).squeeze() * Nnorm * Tnorm * 1.602e-19 * wci
)  # elec -> ion

cmap = plt.cm.get_cmap("bwr")

maxd = max([np.max(np.abs(qi)), np.max(np.abs(wi))]) * 1e-6
mind = -maxd

fig, ax = plt.subplots(2, 1, figsize=(10, 7), gridspec_kw={"height_ratios": [10, 1]})
from mpl_toolkits.axes_grid1 import make_axes_locatable

qi_ax = ax[0]
cbar_ax = ax[1]
divider = make_axes_locatable(qi_ax)
wi_ax = divider.append_axes("right", size="100%", pad=0.01)

c = gridcontourf(grid, qi * 1e-6, ax=qi_ax, show=False, mind=mind, maxd=maxd, cmap=cmap)
qi_ax.set_title(r"Ion$\rightarrow$Neutral transfer [MWm$^{-3}$]")
qi_ax.axes.get_xaxis().set_ticks([1, 2, 3, 4])
qi_ax.set_ylim(bottom=4.0)

qi_ax.xaxis.grid(True)
qi_ax.yaxis.grid(True)

c = gridcontourf(
    grid, -wi * 1e-6, ax=wi_ax, show=False, ylabel=None, mind=mind, maxd=maxd, cmap=cmap
)
wi_ax.set_title(r"Ion$\rightarrow$Electron transfer [MWm$^{-3}$]")
wi_ax.axes.get_yaxis().set_ticklabels([])
wi_ax.axes.get_xaxis().set_ticks([1, 2, 3, 4])
wi_ax.set_ylim(bottom=4.0)

wi_ax.xaxis.grid(True)
wi_ax.yaxis.grid(True)

cbar = fig.colorbar(c, cax=cbar_ax, orientation="horizontal")

fig.tight_layout()

plt.savefig("transfer2d.png")

plt.show()
plt.close("all")
