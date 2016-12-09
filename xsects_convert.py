"""
Cross sections
---------------

Convert cross-section input file to right energy grid
"""

import numpy
import scipy
from numpy import pi, arccos as acos, tan, round, exp, log, log10, sin, cos, logical_and, logical_or
import binning

def interpolate(etarget, e, y):
	ytarget = 10**numpy.interp(x=log10(etarget), xp=log10(e), fp=log10(y))
	ytarget[~numpy.isfinite(ytarget)] = 0
	return ytarget

energy_lo, energy_hi = binning.bin2energy(range(binning.nbins))
energy = (energy_hi + energy_lo)/2.
deltae = energy_hi - energy_lo

nbins = 1250 # number of energy bins

def bin2energy_lo(i):
	i = numpy.asarray(i).reshape((-1,))
	dn = 0.01 * numpy.ones_like(i)
	n = i * dn
	m1 = n > 8.34
	dn[m1] = 0.022
	n[m1] = 8.34 + (i[m1] - 834.0) * dn[m1]
	m2 = n > 9.308
	dn[m2] = 0.05
	n[m2] = 9.308+(i[m2]-878.0)*dn[m2]
	m3 = n > 10.258
	dn[m3] = 0.1
	n[m3] = 10.258 + (i[m3]-897.0)*dn[m3]
	m4 = n > 11.158
	dn[m4] = 0.35
	n[m4] = 11.158+(i[m4] - 906.0)*dn[m4]
	
	return 23.4 * 10**(n/70.0)-23.3

def bin2energy_hi(i):
	return bin2energy_lo(i+1)

i = numpy.arange(1250)
emid = (bin2energy_lo(i) + bin2energy_hi(i)) / 2.

# photoelectric and line cross-sections
xsectsdata = numpy.loadtxt('xsects_orig.dat')
xlines_energies = xsectsdata[0,2:]
xlines_yields = xsectsdata[1,2:]
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
lines = open('xsects_orig.dat').readlines()
nheader = max([i for i, l in enumerate(lines) if l.startswith('#')])
with open('xsects.dat', 'w') as f:
	f.write(''.join(lines[:nheader+1]))
	numpy.savetxt(f, xsects)

