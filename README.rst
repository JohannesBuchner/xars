XARS X-ray Monte-carlo simulator
------------------------------------

XARS simulates X-rays passing through matter in user-defined geometries.

Photo-electric absorption, compton scattering and fluorescent line processes are
modelled.

How to cite XARS correctly
---------------------------

Please reference Buchner et al (in prep).

Outline of the code
----------------------

binning: Defines the energy binning used. 

xsects/__init__.py: computes Compton scattering cross-sections, loads cross-sections for lines from xsects.dat

xsects/xscats.dat contains for each outcoming line its energy, yield, and a table of cross-sections leading to its emission over the energy binning.

If you change the binning, xsects_convert.py can help you rebin xsects_orig.dat into xscats.dat.
If you need different abundances, the Fortran code in xsects/generate/ can help you create a new xsects.dat file.
You will need to insert the yields and line energies in the file header.

To see a example usage, look at disk.py. 
Most of it specifies that the user can run this code e.g. with "python disk.py --nevents 600000 --output=output/disk",
making the output files being stored in output/disk* and using 600000 events in each energy bin for estimating the energy response.

The geometry is loaded "from geometries.disk import DiskGeometry". Several geometries are already defined in geometries.

The binmapfunction defines how output spectra are binned across the sky.

montecarlo.py: montecarlo.run runs the monte-carlo code and returns green function response matrices (These can be stored efficiently into HDF5 files with montecatlo.store).
montecarlo.run goes through each energy bin and creates a package of many photons. The package is pumped through the geometry.
Two things can happen in a pumping step: Either it is still scattering around at some location in the matter of the geometry,
or it escapes to infinity (at which point montecarlo.py records its output energy and direction).
Optionally you can also plot and print what the photons are doing. If you need something more specialised, 
you should write your own version of montecarlo.run.

photons.py: All interactions are modelled here. This is actually fast because we are dealing with large numpy arrays of photon packages at once.
PhotonBunch.__init__ sets up the photon package at the origin in random, isotropic directions.
PhotonBunch.pump asks the geometry to compute how far the photon may travel, computes the interaction taking place.
If compton scattered (marked "stuck"), a new direction is computed in the lower half of PhotonBunch.pump.
PhotonBunch.pump returns photons that escaped to infinite, and should be run in a loop until no photons are left.

Geometries
---------------

It is easy to define geometries in XARS. geometries/disk.py shows an example.
A class needs to be defined with a function compute_next_point.
It receives the photon location (xi, yi, zi) and its direction (dist, beta, alpha).
Dist is how far, in units of column density NH [1e22cm^-2], it should travel.
Given this information, compute_next_point must compute where it ends up, 
in linear (xf,yf,zf) and spherical coordinates (rad, phi, theta)
It also returns whether the photon has left the geometry to infinitey (inside).
All operations work on (large) arrays.

Another example is the sphere geometry in geometries/spheretorus.py. It is 
good practice to visualise the geometry as well. This is 

torus2.py shows how the visualisation is stored in this more elaborate example.

Parallelisation
-------------------

runtorus.sh shows how an array of simulations is run, exploring a grid of 
geometry configurations.

Xspec table models
-------------------

At the bottom of runtorus.sh, commands are shown how to transform rdata output
arrays into fits model tables that xspec can read.
These scripts (in xspecexport, e.g. createtorustable.py) assume a input 
photon spectrum (e.g. a powerlaw) and store the output spectrum into a fits file.



