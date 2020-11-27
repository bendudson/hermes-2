TCV tokamak example
===================

Input grid
----------

- Lower single null plasma
- Horizontal divertor leg
- Separatrix density around 9e19 m^-3

Simulations
-----------

First run the transport case. This includes neutrals, cross-field
diffusion, but no currents or drifts.

    $ mpirun -np 16 ./hermes-2 -d 1-no-currents

This will run for approx. 40 minutes, for a simulation time of 1e5/wci,
sufficient to reach steady state.

