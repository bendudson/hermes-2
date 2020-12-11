#!/usr/bin/env python
#
# Run this with a directory to print the input sources
#
#     $ ./input_source.py   data/

import numpy as np

from boutdata import collect


def total_sources(path, tind=-1):
    """Return total sources integrated over the domain in SI units:
    - Volume [m^3]
    - Particles [s^-1]
    - Power into electrons and ions [Watts]

    Inputs
    ------
    path    The path to the data files
    tind    Time index
    """

    # Normalisation
    qe = 1.602e-19  # Electron charge [C]
    Nnorm = collect("Nnorm", path=path)  # Reference density [m^-3]
    Tnorm = collect("Tnorm", path=path)  # Reference temperature [eV]
    rho_s0 = collect("rho_s0", path=path)  # Reference length [m]
    wci = collect("Omega_ci", path=path)  # Reference frequency [s^-1]

    # Metrics
    J = collect("J", path=path)
    dx = collect("dx", path=path)
    dy = collect("dy", path=path)
    dz = collect("dz", path=path)

    dV = J * dx * dy * dz  # Normalised cell volume

    dV_si = dV * rho_s0 ** 3  # Cell volume [m^3] as a function of [x,y]

    result = {}

    # Read and convert units for each source
    # Note the 3/2 factor multiplying the pressure source to convert to internal energy
    for source_name, norm in [
        ("NeSource", Nnorm * wci),
        ("PeSource", (3.0 / 2) * qe * Nnorm * Tnorm * wci),
        ("PiSource", (3.0 / 2) * qe * Nnorm * Tnorm * wci),
    ]:
        # Read source, select single time, exclude X guard cells, convert to SI units
        source = norm * collect(source_name, path=path, tind=tind)[-1, 2:-2, :, :]

        # Get number of Z points
        nz = source.shape[-1]

        # Multiply by cell volume [m^3]
        for z in range(nz):
            source[:, :, z] *= dV_si[2:-2, :]

        # Sum over all cells
        result[source_name] = np.sum(source)

    result["volume"] = np.sum(dV_si[2:-2, :]) * nz  # Exclude guard cells
    return result


if __name__ == "__main__":

    import sys

    if len(sys.argv) != 2:
        # Print usage information
        print("Usage: {0} path\n e.g. {0} data".format(sys.argv[0]))
        sys.exit(1)

    path = sys.argv[1]
    sources = total_sources(path)

    print("Volume of the simulation: {} m^3".format(sources["volume"]))
    print("Particle source: {} /s".format(sources["NeSource"]))
    print("Electron power: {} W".format(sources["PeSource"]))
    print("Ion power: {} W".format(sources["PiSource"]))
