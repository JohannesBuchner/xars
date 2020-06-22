=================================
UXCLUMPY - Unified Clumpy model
=================================

.. image:: uxclumpy.png
  :target: https://vimeo.com/218031864
  :align: right

This model was presented in `Buchner et al (2019) <https://ui.adsabs.harvard.edu/abs/2019A%26A...629A..16B/abstract>`_. If you have any issues or questions, please check the `FAQ <faq.rst>`_ or open a `Github issue <http://github.com/JohannesBuchner/xars/issues>`_!

To go to the corresponding infrared model go back to `Models <README.rst>`_.

Visualisation
---------------

- 360° VR video: https://vimeo.com/253036759
- normal video: https://vimeo.com/218031864

Components
--------------

Download: https://doi.org/10.5281/zenodo.602282

``atable{uxclumpy.fits}``:

	Clumpy torus transmitted and reflected component with fluorescent lines
	
	Incident radiation parameters:
		
	- PhoIndex: Photon Index
	- Ecut: Energy cutoff [keV]
	- norm: Photon Flux normalisation at 1keV
	
	Viewing angle parameters:
	
	- NHLOS: Total LOS column density
	- Theta_inc: Viewing angle, relative to the inner (flat) disk portion.
	
	Geometry parameters:
	
	- TORsigma: vertical extent of the cloud population. sigma is the width of a gaussian distribution (see `CLUMPY model <https://www.clumpy.org/pages/model-description.html>`_). The number of clouds remains constant, so low sigmas yield slightly higher covering factors (2%-5%).
	- CTKcover: covering factor of inner Compton-thick ring of clouds. If low, many small clouds form a thin ring. If high, few large clouds are used. The column densities of these clouds are logNH=25+-0.5.
	
``atable{uxclumpy-omni.fits}``:

	Warm mirror emission. This is the angle-averaged (omni-directional) spectrum, 
	containing mostly the incident powerlaw from unobscured sightlines.
	
	Given space-filling ionised gas (e.g. in the narrow-line region), 
	Thomson scattering can mirror this emission into otherwise obscured LOS.
	
	The parameters are the same as for the main component, and should always
	be linked. A fraction (with const) should be multiplied onto this component,
	with a maximum of 0.1.

Model setup
-------------

``atable{uxclumpy.fits} + atable{uxclumpy-omni.fits}*const``

Initially, freeze Ecut=400, Theta_inc=90. 

Link the omni component parameters to the main model parameters, const should be free between 1e-5 and 0.1.


Fitting
-------------


AGN obscurer models, including MYTORUS and BNTORUS have highly degenerate parameter spaces.
These are not easy to explore with simple fitting methods, or MCMC.
I generally recommend a global search algorithms. These include Multinest (through BXA).

If you are stuck in a situation without such algorithms, here is some strategies to escape local minima.


1) Freeze some parameters that are less influential. For the uxclumpy model, freeze Ecut=400, and maybe CTKcov=0.4.
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










