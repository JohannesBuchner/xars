XARS X-ray Monte-carlo simulator
------------------------------------

.. image:: logo3-large.png

XARS simulates X-rays propagating through matter in user-defined geometries.

Features:

* Photo-electric absorption
* Compton scattering 
* Fluorescent line emission (self-consistent with absorption above)
* Arbitrary user-defined geometries (included examples: toroid, sphere, disk)
* Arbitrary geometries made from many spherical blobs/clumps (when combined with LightRayRider)
* Arbitrary grid geometries from hydrodynamic simulations (when combined with LightRayRider)

XARS stands for X-ray Absorption, Re-emission and Scattering.

The code-base is small (few hundred lines) and written in pure Python. New contributions are welcome.

How to cite XARS correctly
---------------------------

Please reference Buchner et al (in prep). 

Models
--------------------------------------

In Buchner et al. (in prep) we irradiated the following geometries,
and you can download xspec table models.

* Unified X-ray Clumpy model UXCLUMPY: https://zenodo.org/record/582674
* Warped Disk obscurer: https://zenodo.org/record/823729
* Wada hydrodynamic simulation
* Response of a single blob

See `Models <doc/README.rst>`_

See the paper for description of the parameters and model assumptions.

Usage
---------------------------------------------------

See `Code Tutorial <doc/xars.rst>`_

Questions and Problems
--------------------------------------------

For any questions or problems with the software, please open an issue.
This helps other people google the same question.

License
-------------------

AGPLv3. Contact me if you need a different license.



