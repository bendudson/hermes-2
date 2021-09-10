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

