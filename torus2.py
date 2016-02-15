"""
Monte-Carlo simulator for X-ray obscurer geometries

Literature:

   * Brightman & Nandra (2011)
   * Leahy & Creighton (1993)
"""

import numpy
import scipy
from numpy import pi, arccos as acos, tan, round, log, log10, sin, cos, logical_and, logical_or, arctan as atan
from binning import nbins, energy2bin, bin2energy
import progressbar
import matplotlib.pyplot as plt

energy_lo, energy_hi = bin2energy(range(nbins))
energy = (energy_hi + energy_lo)/2.
deltae = energy_hi - energy_lo

rng = scipy.random

from conetorus import ConeTorusGeometry
import montecarlo

rng.seed(0)

import argparse
import sys, os

parser = argparse.ArgumentParser(
	description="""Monte-Carlo simulator for X-ray obscurer geometries""",
	epilog="""(C) Johannes Buchner, 2013. Based on work by Murray Brightman & Kirpal Nandra (see 2011 publication)""")

parser.add_argument('--log10nh', type=float, help='column density (10^X cm^-1)', default=None)
parser.add_argument('--nh', type=float, help='column density (X/10^22 cm^-1)', default=None)
parser.add_argument('--covering-fraction', type=float, help='cosine of opening angle', default=None)
parser.add_argument('--opening-angle', type=float, help='opening angle in degrees', default=None)
parser.add_argument('--nevents', type=int, default=1000000, help='number of input photons per energy bin')
parser.add_argument('--plot-paths', default=False, help='plot the paths taken?', action='store_true')
parser.add_argument('--plot-interactions', default=False, help='plot the points at each interaction?', action='store_true')
parser.add_argument('--verbose', default=False, help='Be more talkative, show debug statistics on interactions', action='store_true')
parser.add_argument('--output', type=str, default=None, help='Prefix for output files. Chosen based on log10nh and opening angle if not specifified')
args = parser.parse_args()

nmu = 10     # number of viewing angle bins

if args.log10nh is not None:
	nh_in = args.log10nh
	nh = 10**(nh_in-22)   # current NH value
else:
	if args.nh is None:
		sys.stderr.write("ERROR: NH not given.\n\n")
		parser.print_help()
		sys.exit(-1)
	nh = args.nh
	nh_in = log10(nh)+22
print '  NH                : 10^%.1f (%.3fe22)' % (nh_in, nh)
if args.covering_fraction is not None:
	cone = acos(args.covering_fraction) # current torus opening angle
	cone_in = args.covering_fraction
else:
	if args.opening_angle is None:
		sys.stderr.write("ERROR: opening angle not given.\n\n")
		parser.print_help()
		sys.exit(-1)
	cone = args.opening_angle / 180. * pi # current torus opening angle
	cone_in = cos(cone)
print '  opening angle: %.1f degrees (covering fraction %.1f)' % (cone * 180. / pi, cone_in)
assert cone * 180. / pi >= 0, 'Opening angle should be above 0 degrees'
assert cone * 180. / pi <= 90, 'Opening angle should be below 90 degrees'
assert nh_in >= 19, 'NH should be above 20'
assert nh_in <= 27, 'NH should be below 27'

prefix = '%.1f_%.1f_' % (nh_in, cone_in) if args.output is None else args.output
  # total number of photons to send in

geometry = ConeTorusGeometry(Theta_tor = cone, NH = nh, verbose=args.verbose)
geometry.viz()
plt.savefig(prefix + "geometry.pdf")
plt.savefig(prefix + "geometry.png")
plt.close()

#binmapfunction = lambda beta: numpy.round(0.5 + nmu * numpy.abs(cos(beta))) - 1
def mapper(beta):
	#beta[beta > pi/2] = pi - beta[beta > pi/2]
	slot = numpy.floor(nmu * beta / pi)
	#print beta * 180 / pi, slot
	return slot

binmapfunction = lambda beta: numpy.floor((nmu - 2) * beta / pi)
binmapfunction = mapper

rdata, nphot = montecarlo.run(prefix, nphot = args.nevents, nmu = nmu, geometry=geometry, 
	binmapfunction = binmapfunction,
	plot_paths=args.plot_paths, plot_interactions=args.plot_interactions, verbose=args.verbose)

header = dict(NH='%f' % nh, OPENING='%f' % cone)
montecarlo.store(prefix, nphot, rdata, nmu, extra_fits_header = header)




