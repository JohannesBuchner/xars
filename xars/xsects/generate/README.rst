===========================
Cross-section computation
===========================

The xsects program computes cross-sections as a function of energy for
photo-electric absorption and fluorescent emission (including line energies and yields).


Compiling and Running
=======================

Compile with::

	$ make 

Run with::

	$ ./xsects 
	SYNOPSIS: xsects <Feabundance> <Zabundance>

	Abundances are relative. Use 1 for local ISM abundances.

To create a xsects file with standard abundances, use::

	$ ./xsects 1.0 1.0

This produces a xsects.dat file, which should be placed in the xsects/ folder.

Caveats
===========

The total photo-electric cross-section has a bad behaviour at E>70keV. This is
corrected in XARS when loading the xsects.dat file.

See Brightman et al. 2011 for the definition of abundances, cross-sections, 
energies and yields.

Updated and alternative cross-sections are welcome.


