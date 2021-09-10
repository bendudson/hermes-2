#!/usr/bin/env python3
#
# Change profiles in input grid
# - Keep pressure gradient the same, so Grad-Shafranov is unchanged
# - Set separatrix density and temperature

import argparse

parser = argparse.ArgumentParser(description="Set profiles in a BOUT++ grid")
parser.add_argument("gridfile", type=str, help="The file to modify")
parser.add_argument("Tsep", type=float, help="Separatrix temperature [eV]")
parser.add_argument("Nsep", type=float, help="Separatrix density [m^-3]")

args = parser.parse_args()
filename = args.gridfile
Tsep = args.Tsep
Nesep = args.Nsep

from boututils.datafile import DataFile
import numpy as np

qe = 1.602e-19

with DataFile(filename) as f:
    P = 0.0
    if "pressure" in f.keys():
        P = f["pressure"]  # In Pascals
    # Check that P's not zero. This is because some grid files
    # may have pressure = 0.

    if np.amax(P) < 1e-3:
        print("WARNING: input grid has no pressure, or pressure is zero")
        Ne = f["Ni0"] * 1e20  # Note: Grid file in 1e20 m^-3
        Te = f["Te0"]
        Ti = f["Ti0"]
        P = qe * Ne * (Te + Ti)

    ixseps = f["ixseps1"]
    # Y index around the midplane. Doesn't matter as long as it's in the core region
    yind = int(0.5 * (f["jyseps1_2"] + f["jyseps2_2"] + 1))

# Subtract pressure, to get required value at separatrix
# Clip to impose minimum of zero
P = np.clip(P - (P[ixseps, yind] - qe * Nesep * 2.0 * Tsep), 1e-3, None)

shape = np.sqrt(P) / np.sqrt(P[ixseps, yind])  # 1 at the separatrix

Te = shape * Tsep
Ti = Te
Ne = shape * (Nesep * 1e-20)

with DataFile(filename, write=True, create=False) as f:
    f["Ni0"] = Ne
    f["Ne0"] = Ne
    f["Te0"] = Te
    f["Ti0"] = Ti
