XARS X-ray Monte-carlo simulator
------------------------------------

.. image:: ../logo3-mid.png
  :align: right

XARS simulates X-rays propagating through matter in user-defined geometries.

This code tutorial explains how to use and modify XARS.

To find existing models, go back to `Models <README.rst>`_.

Tutorial: Part I: Irradiating a programmed geometry
---------------------------------------------------

In part I we look at a irradiating a user-specified geometry. Part II will look at 
a geometry made up of spherical clumps/blobs. 

The bash script runsphere.sh simulates a spherical obscurer with various
column densities. To run a single one, run, for example::

	$ python torus2.py --log10nh=24.2 --opening-angle=0 --nevents=1000000 --output=myoutput

You can study torus2.py, and in particular its geometry definition geometries/spheretorus.py
to understand how the code works. See below for detailed description.


Tutorial: Part II: Irradiating a geometry made up of spheres
---------------------------------------------------------------

In part II, we assume that your geometry can be expressed as many spheres.

The example-blobs/generate_blobs.py demonstrates how to generate a hdf5 input 
file which describes the blobs, with their x/y/z positions and column densities.
In this simple case, there is only one; in general, just enlarge the arrays.

To irradiate such models, you need to install the LightRayRider library
which performs fast photon propagation via optimized C functions.
Download from https://github.com/JohannesBuchner/LightRayRider, for example to
$HOME/Downloads/LightRayRider/. Then compile with::

	$ make -C $HOME/Downloads/LightRayRider/

This will create a ray.so object.

To irradiate the output files, e.g. torusblob23.0.hdf5, run::

	$ PYTHONPATH=$HOME/Downloads/LightRayRider/ python torusC.py --geometry=torusblob23.0.hdf5 --nevents=1000000

To use parallelise over 10 CPUs, run with::

	$ OMP_NUM_THREADS=10 PYTHONPATH=$HOME/Downloads/LightRayRider/ python torusC.py --geometry=torusblob23.0.hdf5 --nevents=1000000

Tutorial: Part III: Irradiating a simulation grid
-------------------------------------------------------------

In part III we assume that you have created a hydrodynamic simulation on a 
3d uniform grid, and want to irradiate this with X-rays.

The example-grid/generate_warpeddisk.py demonstrates how to generate a hdf5 input 
file which describes the grid and its density, as well as the irradiation 
location.

To irradiate such models, you need to install the LightRayRider library
which performs fast photon propagation via optimized C functions.
Download from https://github.com/JohannesBuchner/LightRayRider, for example to
$HOME/Downloads/LightRayRider/. Then compile with::

	$ make -C $HOME/Downloads/LightRayRider/

This will create a ray.so object.

To irradiate the output files, e.g. warpeddisk_1.hdf5, run::

	$ PYTHONPATH=$HOME/Downloads/LightRayRider/ python torusG.py --geometry=warpeddisk_1.hdf5 --nevents=1000000

To use parallelise over 10 CPUs, run with::

	$ OMP_NUM_THREADS=10 PYTHONPATH=$HOME/Downloads/LightRayRider/ python torusC.py --geometry=warpeddisk_1.hdf5 --nevents=1000000



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

* Irradiating different geometries is embarrassingly parallel. 
* For irradiating the same geometry, XARS can take advantage of multiple CPUs (see OMP_NUM_THREADS).
* To parallelise over multiple machines, make sure the output files are named differently. You can combine the rdata output files with the rdataaddmultiple.py script.

Xspec table models
-------------------

At the bottom of runtorus.sh, commands are shown how to transform rdata output
arrays into fits model tables that xspec can read.
These scripts (in the xspecexport folder, e.g. createtorustable.py) assume a input 
photon spectrum (e.g. a powerlaw) and store the output spectrum into a fits file.
Adjust to additional parameters and input spectra as needed.

Questions and Problems
--------------------------------------------

For any questions or problems with the software, please open an issue.
This helps other people google the same question.




