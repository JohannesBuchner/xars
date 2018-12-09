==================
Warped disk model
==================

If you have any issues or questions, please check the FAQ or open a github issue!

To go to the corresponding infrared model go back to ´models.rst´_.

Visualisation: https://zenodo.org/record/1193541/files/warpgeometry.pdf

Components
--------------

Download: https://doi.org/10.5281/zenodo.823728

``atable{warpeddisk.fits}``:

	Warped disk transmitted and reflected component with fluorescent lines
	
	Incident radiation parameters:
		
	- PhoIndex: Photon Index
	- Ecut: Energy cutoff [keV]
	- norm: Photon Flux normalisation at 1keV
	
	Viewing angle parameters:
	
	- NHLOS: LOS column density
	- Theta_inc: Viewing angle, relative to the inner (flat) disk portion.
	
	Geometry parameters:
	
	- NHdisk: Column density through the disk (in particular at the center). In log
	- diskfrac: Warp extent, i.e. how far the disk reaches. 0.1 means a relatively flat disk, 1 is a stronger warp.
	
``atable{warpeddisk-omni.fits}``:

	Warm mirror emission. This is the angle-averaged (omni-directional) spectrum, 
	containing mostly the incident powerlaw from unobscured sightlines.
	
	Given space-filling ionised gas (e.g. in the narrow-line region), 
	Thomson scattering can mirror this emission into otherwise obscured LOS.
	
	The parameters are the same as for the main component, and should always
	be linked. A fraction (with const) should be multiplied onto this component,
	with a maximum of 0.1.

Model setup
-------------

1) **Heavily Obscured configuration**::
	Here you view Compton-thick through the disk, edge-on to the inner disk.
	
	``atable{warpeddisk.fits} + atable{warpeddisk-omni.fits}*const``
	
	NHLOS=1000, PhoIndex=2, Ecut=400, diskfrac=0.5, Theta_inc=90
	
	Freeze Ecut=400, Theta_inc=90. Link the omni component parameters to the main model, const should be free between 1e-5 and 0.1.

2) **Mildly Obscured configuration**::
	The warped disk can not produce Compton-thin obscured AGN by itself.
	Clouds in the host ISM (in particular the inner 100s of pc), can provide 
	this obscuration.
	Here you view the disk unobscured (e.g. perpendicular), with some polar absorber.
	
	``tbabs * atable{warpeddisk.fits} + atable{warpeddisk-omni.fits} * const``
	
	NHTBABS=1, NHLOS=0, PhoIndex=2, Ecut=400, diskfrac=0.5 Theta_inc=0
	
	Freeze Ecut=400, NHLOS=0. Link the omni component parameters to the main model, const should be free between 1e-5 and 0.1.

Also try 1) with a different viewing angle (e.g. polar, Theta_inc=0).


Fitting
-------------


AGN obscurer models, including MYTORUS and BNTORUS have highly degenerate parameter spaces.
These are not easy to explore with simple fitting methods, or MCMC.
I generally recommend a global search algorithms. These include Multinest (through BXA).

If you are stuck in a situation without such algorithms, here is some strategies to escape local minima.


1) Freeze some parameters that are less influential. For the warpeddisk model, freeze Ecut=400, NHdisk=25 and diskfrac=0.5. 
2) Limit the data. Local minima are difficult to escape because they are surrounded by steep walls. Use fewer spectra, and start with a high-energy data (e.g. 20-50keV). Freeze the omni constant to 1e-10. Fit and gradually add more data (15keV, 8keV, 5keV, etc). Then allow the omni constant to vary.
3) Use the error command to explore more. This is most helpful on NH and geometry parameters.
4) Plot the model -- if the model component is not present, it is because this viewing angle and this NHLOS does not exist in this geometry. You need to change the geometry or viewing angle until you get back the model.

Additionally, some parameter regions are simply discontinuous.

1) Try freezing NH to Compton-thin (e.g. 30) and refit using the above steps, then thaw and fit.
2) Try freezing NH to Compton-thick (e.g. 300) and refit using the above steps, then thaw and fit.
3) Try freezing NH to 0 and use a LOS absorber (see Model setup).



Failure states
---------------

- If you get a very low photon index (<1.5), you are probably in a bad local minimum. Start from scratch. Maybe freeze PhoIndex=2 and see how far you get.

- Plot the model. If the AGN obscurer component is not present, it is because this viewing angle and this NHLOS does not exist in this geometry (zero photons). You need to change the geometry or viewing angle until you get back the model.










