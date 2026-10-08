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


def interpolate(etarget, e, y, nslope=50):
    x = log10(e)
    xtarget = log10(etarget)
    logy = log10(y)
    logytarget = numpy.interp(x=xtarget, xp=x, fp=logy)
    # continue with a power law beyond the end of the source table
    # (the extrapolated photoelectric cross-sections are orders of magnitude
    # below the analytically computed scattering cross-section there)
    beyond = xtarget > x[-1]
    if beyond.any() and numpy.isfinite(logy[-1]):
        finite = numpy.isfinite(logy)
        if finite.sum() >= 2:
            xf, yf = x[finite][-nslope:], logy[finite][-nslope:]
            slope = numpy.polyfit(xf, yf, 1)[0]
            logytarget[beyond] = logy[finite][-1] + slope * (xtarget[beyond] - x[-1])
    ytarget = 10**logytarget
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
xsects = xsectsdata[4:,:]
e1 = xsects[:,0]
assert len(e1) == len(emid), (len(e1), len(emid))
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

# write out, keeping the comments and the line energy/yield/width/asymmetry
# header rows
with open('xsects_orig.dat') as fin:
    lines = fin.readlines()
data_rows = [i for i, l in enumerate(lines) if not l.startswith('#')]
nheader_rows = 4  # line energies, yields, fwhm, asymmetries
table_start = data_rows[nheader_rows]
with open('xsects.dat', 'w') as f:
    f.write(''.join(lines[:table_start]))
    numpy.savetxt(f, xsects)
