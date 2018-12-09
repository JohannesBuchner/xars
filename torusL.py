from __future__ import print_function, division
"""
Monte-Carlo simulator for X-ray obscurer geometries

Literature:

   * Brightman & Nandra (2011)
   * Leahy & Creighton (1993)
"""

import numpy
import scipy
from numpy import pi, arccos as acos, tan, round, log, log10, sin, cos, logical_and, logical_or, arctan as atan
import progressbar
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

rng = scipy.random

from geometries.layeredconetorus import LayeredConeTorusGeometry
import montecarlo

#rng.seed(0)

import argparse
import sys, os

parser = argparse.ArgumentParser(
	description="""Monte-Carlo simulator for X-ray obscurer geometries""",
	epilog="""(C) Johannes Buchner, 2013-2016. Based on work by Murray Brightman & Kirpal Nandra (see 2011 publication)""")

parser.add_argument('--geometry', type=str, required=True, help='Geometry file')
parser.add_argument('--nevents', type=int, default=1000000, help='number of input photons per energy bin')
parser.add_argument('--verbose', default=False, help='Be more talkative, show debug statistics on interactions', action='store_true')
args = parser.parse_args()


Theta_tors, NHs = numpy.loadtxt(args.geometry).transpose()
nmu = len(Theta_tors)     # number of viewing angle bins
prefix = args.geometry + 'layered_'

geometry = LayeredConeTorusGeometry(Theta_tors = Theta_tors, NHs = NHs, verbose=args.verbose)
geometry.viz()
plt.savefig(prefix + "geometry.pdf")
plt.savefig(prefix + "geometry.png")
plt.close()

def binmapfunction(beta, alpha):
	beta = numpy.where(beta > pi/2, pi - beta, beta)
	slot = numpy.zeros(beta.shape, dtype=int)
	for i, theta in enumerate(Theta_tors[:-1]):
		slot[beta > theta] = i+1
	return slot

rdata, nphot = montecarlo.run(prefix, nphot = args.nevents, nmu = nmu, geometry=geometry, 
	binmapfunction = binmapfunction,
	plot_paths=False, plot_interactions=False, verbose=args.verbose)

rdata_transmit, rdata_reflect = rdata
header = dict()
montecarlo.store(prefix + 'transmit', nphot, rdata_transmit, nmu, extra_fits_header = header, plot=False)
montecarlo.store(prefix + 'reflect', nphot, rdata_reflect, nmu, extra_fits_header = header, plot=False)
rdata_transmit += rdata_reflect
del rdata_reflect
montecarlo.store(prefix, nphot, rdata_transmit, nmu, extra_fits_header = header, plot=True)




