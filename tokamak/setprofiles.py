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

with DataFile(filename) as f:
    Ne = f["Ni0"]  # Note: In 1e20 m^-3
    Te = f["Te0"]
    Ti = f["Ti0"]
    ixseps = f["ixseps1"]
    # Y index around the midplane. Doesn't matter as long as it's in the core region
    yind = int(0.5 * (f["jyseps1_2"] + f["jyseps2_2"] + 1))

# Total pressure
P = Ne * (Te + Ti)

# Subtract pressure, to get required value at separatrix
# Clip to impose minimum of zero
P = np.clip(P - (P[ixseps, yind] - (Nesep * 1e-20) * 2.0 * Tsep), 1e-3, None)

shape = np.sqrt(P) / np.sqrt(P[ixseps, yind])  # 1 at the separatrix

Te = shape * Tsep
Ti = Te
Ne = shape * (Nesep * 1e-20)

with DataFile(filename, write=True, create=False) as f:
    f["Ni0"] = Ne
    f["Ne0"] = Ne
    f["Te0"] = Te
    f["Ti0"] = Ti
