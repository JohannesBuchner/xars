===========================
CAT3D+WIND clumpy model
===========================

Contact: Johannes Buchner <johannes.buchner.acad@gmx.com>

Geometry information: http://www.sungrazer.org/cat3d.html

If you have any issues or questions, please check the `FAQ <faq.rst>`_ or open a `Github issue <http://github.com/JohannesBuchner/xars/issues>`_!

To go to the corresponding infrared model go back to `Models <README.rst>`_.

(This geometry is the one of `Leftley et al. 2018<http://adsabs.harvard.edu/abs/2018ApJ...862...17L>`_, matching e.g., ESO323-G77)

Components
--------------

Download: https://zenodo.org/record/2211263

``atable{CAT3D-WIND.fits}``:

	Clumpy model transmitted and reflected component with fluorescent lines
	
	Incident radiation parameters:
		
	- PhoIndex: Photon Index
	- Ecut: Energy cutoff [keV]
	- norm: Photon Flux normalisation at 1keV
	
	Viewing angle parameters:
	
	- NHLOS: Total LOS column density
	- Theta_inc: Viewing angle, relative to the inner (flat) disk portion.
	
	Geometry parameters:
	
	- NH_cloud: Column density of each cloud. In log
	
Model setup
-------------

``atable{CAT3D-WIND.fits} + zpow*const``

Initially, freeze Ecut=400, Theta_inc=90. 

Link the powerlaw parameters to the clumpy model parameters, const should be free between 1e-5 and 0.1.

Fitting
-------------

AGN obscurer models, including MYTORUS and BNTORUS have highly degenerate parameter spaces.
These are not easy to explore with simple fitting methods, or MCMC.
I generally recommend a global search algorithms. These include Multinest (through BXA).

If you are stuck in a situation without such algorithms, here is some strategies to escape local minima.


1) Freeze some parameters that are less influential. For this model, freeze Ecut=400, Theta_inc=90, NH_cloud=24. 
2) Limit the data. Local minima are difficult to escape because they are surrounded by steep walls. Use fewer spectra, and start with a high-energy data (e.g. 20-50keV). Freeze the powerlaw constant to 1e-10. Fit and gradually add more data (15keV, 8keV, 5keV, etc). Then allow the powerlaw constant to vary.
3) Use the error command to explore more. This is most helpful on NH and geometry parameters.
4) Plot the model -- if the model component is not present, it is because this viewing angle and this NHLOS does not exist in this geometry. You need to change the geometry or viewing angle until you get back the model.

Additionally, some parameter regions are simply discontinuous.

1) Try freezing NH to Compton-thin (e.g. 30) and refit using the above steps, then thaw and fit.
2) Try freezing NH to Compton-thick (e.g. 300) and refit using the above steps, then thaw and fit.


Failure states
---------------

- If you get a very low photon index (<1.5), you are probably in a bad local minimum. Start from scratch. Maybe freeze PhoIndex=2 and see how far you get.

- Plot the model. If the AGN obscurer component is not present, it is because this viewing angle and this NH does not exist in this geometry (zero photons). You need to change the geometry or viewing angle until you get back the model.










