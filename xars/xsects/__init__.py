"""
Cross sections
---------------

Loading and computation of neutral, solar-abundance cross-sections.
"""

import os

import numpy
import scipy
from numpy import arccos as acos
from numpy import (cos, exp, log, log10, logical_and, logical_or, pi, round,
                   sin, tan)

from xars import binning

electmass = 511.  # electron rest mass in keV/c^2

xscatt = numpy.zeros(binning.nbins)

energy_lo, energy_hi = binning.bin2energy(numpy.arange(binning.nbins))
energy = (energy_hi + energy_lo) / 2.
deltae = energy_hi - energy_lo

xthom = 6.7e-4          # Thomson cross section
x = energy / electmass  # energy in units of electron rest mass

# compute scattering cross section
# Thomson regime
xscatt_thomson = xthom * (1. - (2. * x) + (26. * (x**2.) / 5.))
t1 = 1. + 2. * x
t2 = 1. + x
term1 = t2 / (x**3.)
term2 = (2. * x * t2) / t1
term3 = (1. / (2. * x)) * log(t1)
term4 = (1. + 3. * x) / (t1**2)
xscatt_compton = 0.75 * xthom * (term1 * (term2 - log(t1)) + term3 - term4)
xscatt = numpy.where(x < 0.05, xscatt_thomson, xscatt_compton)

# convert to units of 1e-22 cm^2 (from 1e-21)
xscatt *= 10
xscatt_thomson *= 10
xscatt_compton *= 10

# When applied the cross section is 120% larger
xscatt_thomson *= 1.2
xscatt_compton *= 1.2
xscatt *= 1.2

# photoelectric and line cross-sections
# find xsects.dat file next to this file
xsectsdata = numpy.loadtxt(os.path.join(os.path.dirname(__file__), 'xsects.dat'))
xlines_energies = xsectsdata[0,2:]
xlines_yields = xsectsdata[1,2:]
xsects = xsectsdata[2:,:]
# convert to units of 1e-22 cm^2 (from 1e-21)
xsects[:,1:] *= 10
e1 = xsects[:,0]
xphot = xsects[:,1]
e70 = e1 > 70.
if e70.any():
    lines_max = numpy.max(xsects[:,2:], axis=1)
    xphot[e70] = lines_max[e70] * xphot[e70][0] / lines_max[e70][0]
assert (xphot >= 0).all()
assert (xscatt >= 0).all()

xlines = xsects[:,2:] * xlines_yields
xlines_relative = xlines / xphot[:,None]
# compute probability to come out as fluorescent line
xlines_cumulative = numpy.cumsum(xlines_relative, axis=1)
assert (xlines >= 0).all()
assert (xlines_relative >= 0).all()
assert (xlines_cumulative >= 0).all()

xboth = xphot + xscatt
absorption_ratio = xphot / xboth


def test():
    assert e1.shape == energy.shape, (e1.shape, energy.shape)
    for i in range(len(energy)):
        assert (numpy.isclose(e1[i], energy[i])), (e1[i], energy[i], energy_lo[i], energy_hi[i])

