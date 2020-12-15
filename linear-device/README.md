Linear device simulations
=========================

Simulations without magnetic curvature. The diamagnetic current
term is usually turned off, because it has zero divergence.

Annulus
-------

Sets up a simulation where the Z coordinate is the azimuthal angle,
and X is a radial coordinate. This coordinate system has a pole
at the axis, and so X covers a range of radius with a hole in the middle.

The parameters are set up to be approximately LAPD-like.

To run the simulation:

    $ mpirun -np 16 ./hermes-2 -d annulus


Annulus with neutrals
---------------------

This is the same geometry as the Annulus case, but has a fixed neutral
background. The density and pressure of the neutrals are controlled by the `Nn:function`
and `Pn:function` settings. The mean free path of the neutrals is assumed to be much
larger than the plasma, so the effect of the plasma on the neutrals is neglected.

This simulation has no external particle source, only an external power
input to the electrons. The particle source is from ionisation of the neutral
background. Recombination and charge exchange are also included.

To run the simulation:

    $ mpirun -np 16 ./hermes-2 -d annulus-with-neutrals

