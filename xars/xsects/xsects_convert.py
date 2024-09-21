"""
Cross sections
---------------

Convert cross-section input file to right energy grid
"""

import numpy
import scipy
from numpy import arccos as acos
from numpy import (cos, exp, log, log10, logical_and, logical_or, pi, round,
                   sin, tan)

from xars import binning
from xars.binning.bn import bin2energy_lo, nbins


def interpolate(etarget, e, y):
    ytarget = 10**numpy.interp(x=log10(etarget), xp=log10(e), fp=log10(y))
    ytarget[~numpy.isfinite(ytarget)] = 0
    return ytarget


energy_lo, energy_hi = binning.bin2energy(numpy.arange(binning.nbins))
energy = (energy_hi + energy_lo) / 2.
deltae = energy_hi - energy_lo


def bin2energy_hi(i):
    return bin2energy_lo(i + 1)


i = numpy.arange(nbins)
emid = (bin2energy_lo(i) + bin2energy_hi(i)) / 2.

# photoelectric and line cross-sections
xsectsdata = numpy.loadtxt('xsects_orig.dat')
xlines_energies = xsectsdata[0, 2:]
xlines_yields = xsectsdata[1, 2:]
xsects = xsectsdata[2:,:]
e1 = xsects[:,0]
assert len(e1) == len(emid)
xphot = xsects[:,1]
e70 = e1 > 70.
lines_max = numpy.max(xsects[:,2:], axis=1)
xphot[e70] = lines_max[e70] * xphot[e70][0] / lines_max[e70][0]
assert (xphot >= 0).all()
xsects[:,1] = xphot

# now rebin
xsects_orig = xsects
xsects = [energy]
for i in range(1, xsects_orig.shape[1]):
    xsects.append(interpolate(energy, emid, xsects_orig[:,i]))
xsects = numpy.transpose(xsects)

# write out
with open('xsects_orig.dat') as fin:
    lines = fin.readlines()
nheader = max([i for i, l in enumerate(lines) if l.startswith('#')])
with open('xsects.dat', 'w') as f:
    f.write(''.join(lines[:nheader + 1]))
    numpy.savetxt(f, xsects)
