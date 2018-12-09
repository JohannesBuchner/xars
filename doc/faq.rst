===========================
Frequently asked questions
===========================

What does the NH parameter mean?
---------------------------------

This is not a property of the geometry, but the line-of-sight (LOS) column density.

Basically, the sky as seen from the corona is segmented by column density (and, corsely, by viewing angle).
Output spectra are accumulated in these column density bins.

This makes spectral fitting less degenerate, 
because LOS NH predominantly shapes the spectrum.

It also allows very similar viewing angles to have very different column densities.

How is the photon path solved in XARS? Is it a Monte Carlo integration?
------------------------------------------------------------------------

No. The optical depth a photon travels is a Monte carlo draw, however the end point
is always determined analytically.

In simple geometries (such as spheres, cones), the necessary line integral 
can be computed and programmed analytically.

For geometries with many spheres, there is optimized, parallelised C code
to determine which spheres intersect, and how they are ordered. 
This is implemented in https://github.com/JohannesBuchner/LightRayRider

For general geometries based on grids, there is optimized, parallelised C code
to traverse the grid.
This is implemented in https://github.com/JohannesBuchner/LightRayRider


I have a question
---------------------

Please open a `Github issue <http://github.com/JohannesBuchner/xars/issues>`_ 
or, if it should not be public, send me an email.



