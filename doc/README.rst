====================================
XARS X-ray Monte-carlo simulator
====================================

XARS simulates X-rays passing through matter in user-defined geometries.

Photo-electric absorption, compton scattering and fluorescent line processes are
modelled.


Usage
--------------------------
see `Code Documentation <xars.rst>`_

How to cite XARS correctly
---------------------------

Please reference Buchner et al (in prep).


Models
==================

In Buchner et al. (in prep) we irradiated the following geometries,
and you can download xspec table models for them. 

They also have infrared models associated with them.

See the paper for description of the parameters and model assumptions.


UXCLUMPY
--------------------

.. image:: uxclumpy.png
  :target: https://vimeo.com/218031864
  :align: right

The Unified X-ray Clumpy model (UXCLUMPY) features:

* Unification: Can produce unobscured, obscured and Compton-thick AGN in the right proportions.
* Eclipse events: The frequency and sizes of eclipses are reproduced by clumps on Keplerian orbits.
* X-ray spectra can fit observations well
* Compatible with CLUMPY infrared models

Here you can access:

* Geometry movies: https://vimeo.com/218031864 and 360 VR https://vimeo.com/253036759
* X-ray table model available at: `UXCLUMPY page <uxclumpy.rst>`_
* Infrared model available at: http://clumpy.org 
* More information: Buchner et al. (submitted) (or send me an email)

Warped disk
--------------------

.. image:: warpgeometry.png
  :target: warpeddisk.rst
  :align: right

A simple warped disk geometry

* Geometry images: `Warped Disk page <warpeddisk.rst>`_
* X-ray table model available at: `Warped Disk page <warpeddisk.rst>`_
* Infrared model: see `Jud et al (2017) <http://cdsads.u-strasbg.fr/abs/2017MNRAS.465..248J>`_
* More information: Buchner et al. (in prep) (or send me an email)


Radiative fountain (Wada 2012)
-------------------------------

.. image:: wadageometry.png
  :target: wada.rst
  :align: right

* Geometry images: `Radiative fountain page <wada.rst>`_
* X-ray table model available at: `Radiative fountain page <wada.rst>`_
* Infrared model: `Radiative fountain page <wada.rst>`_
* More information: Buchner et al. (in prep) (or send me an email)

CAT3D+WIND
---------------------------

* Geometry images: ...
* Infrared model: http://www.sungrazer.org/cat3d.html
* X-ray model: To be uploaded.



Response of a single spherical blob
-------------------------------------

Reflection from a single sphere with

* Isotropic density
* Exponential density profile
* Gaussian density profile

X-ray spectral model available on request (open a `Github issue <http://github.com/JohannesBuchner/xars/issues>`_).

* Geometry images: (imagine a sphere)
* X-ray table model available on request (open a `Github issue <http://github.com/JohannesBuchner/xars/issues>`_).
* Infrared model: https://en.wikipedia.org/wiki/Planck%27s_law



