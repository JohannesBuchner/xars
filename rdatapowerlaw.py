from __future__ import print_function, division
"""
converts a rdata fits file into spectra using a input powerlaw
"""

import numpy
import scipy
from numpy import pi, arccos, tan, round, exp, log, log10, sin, cos, logical_and, logical_or, arctan
from binning import nbins, energy2bin, bin2energy
import matplotlib.pyplot as plt

energy_lo, energy_hi = bin2energy(numpy.arange(nbins))
energy = (energy_hi + energy_lo)/2.
deltae = energy_hi - energy_lo

import argparse
import sys, os

parser = argparse.ArgumentParser(
	description="""Spectrum generator using powerlaw input onto rdata matrix""",
	epilog="""(C) Johannes Buchner, 2013 """)

parser.add_argument('--PhoIndex', type=float, help='Photon Index / Gamma', required=True)

#parser.add_argument('--Ecut', type=float, help='Cutoff energy', default=numpy.inf)
parser.add_argument('--rdata', type=str, help='energy/energy/viewing angle matrix in fits format', required=True)

parser.add_argument('--prefix', type=str, default=None, help='output file prefix')

args = parser.parse_args()

prefix = args.rdata + "_" if args.prefix is None else args.prefix
PhoIndex = args.PhoIndex

import pyfits
f = pyfits.open(args.rdata)
rdata = f[0].data
nphot = f[0].header['NPHOT']

nmu = rdata.shape[2]

"""
Integral of x^-G from a to b:
a/(a^y*y-a^y)-b/(b^y*y-b^y)
"""
#a = self.xdata - self.xdelta/2.
#b = self.xdata + self.xdelta/2.
#G = self.PhoIndex.val

weights = (
	energy_lo / (energy_lo**PhoIndex * PhoIndex - energy_lo**PhoIndex) - 
	energy_hi / (energy_hi**PhoIndex * PhoIndex - energy_hi**PhoIndex))
#A, B = numpy.meshgrid(deltae, deltae)
weights = energy**-PhoIndex

plt.figure(figsize=(7,5))

results = []

for mu in range(nmu):
	print('viewing pos', mu)
	# we assume a quadratic matrix, i.e. x * x
	o = []
	for j in range(rdata.shape[1]):
		oi = (rdata[:,j,mu] / nphot * energy**-PhoIndex * deltae).sum() * 10
		o.append(oi / deltae[j])
	
	#matrix = numpy.multiply(weights, rdata[:,:,mu] / nphot)
	#print weights, rdata[:,:,mu], matrix
	
	#y = matrix.sum(axis=1)
	y = numpy.array(o)
	results.append(y)
	plt.plot(energy, y, label='viewing bin %d' % mu)
plt.xlabel('energy')
plt.ylabel('photon flux [cts]')
plt.xlim(0.1, 100)
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.ylim(1e4/1e9, 1e11/1e9)
plt.legend(loc='best', prop=dict(size=6))
plt.savefig(prefix + "spec.pdf")
plt.savefig(prefix + "spec.png")
plt.close()
	
numpy.savetxt(prefix + "spec.txt", numpy.vstack([energy_lo, energy_hi] + results).transpose())






