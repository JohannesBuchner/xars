"""
Monte-Carlo simulator for X-ray obscurer geometries

Literature:

   * Brightman & Nandra (2011)
   * Leahy & Creighton (1993)
"""

import numpy
from numpy import cos
import matplotlib as mpl
mpl.use('Agg')

rng = numpy.random

from xars.geometries.disk import DiskGeometry
from xars import montecarlo

import argparse

parser = argparse.ArgumentParser(
	description="""Monte-Carlo simulator for X-ray obscurer geometries""",
	epilog="""(C) Johannes Buchner, 2013. Based on work by Murray Brightman & Kirpal Nandra (see 2011 publication)""")

parser.add_argument('--nevents', type=int, default=1000000, help='number of input photons per energy bin')
parser.add_argument('--plot-paths', default=False, help='plot the paths taken?', action='store_true')
parser.add_argument('--plot-interactions', default=False, help='plot the points at each interaction?', action='store_true')
parser.add_argument('--verbose', default=False, help='Be more talkative, show debug statistics on interactions', action='store_true')
parser.add_argument('--output', type=str, default='disk_', help='Prefix for output files. ')
parser.add_argument('--plot-every', type=int, default=1, help='Plot only for every Nth energy bins')
args = parser.parse_args()

nmu = 10     # number of viewing angle bins
geometry = DiskGeometry(verbose=args.verbose)
prefix = args.output

def binmapfunction(beta, alpha): 
	return (numpy.round(0.5 + nmu * numpy.abs(cos(beta))) - 1).astype(int)

rdata, nphot = montecarlo.run(prefix, nphot = args.nevents, nmu = nmu, geometry=geometry, 
	binmapfunction = binmapfunction,
	plot_paths=args.plot_paths, plot_interactions=args.plot_interactions, verbose=args.verbose,
	plot_every=args.plot_every)

rdata_transmit, rdata_reflect = rdata
rdata_both = rdata_transmit + rdata_reflect
montecarlo.store(prefix + 'transmit', nphot, rdata_transmit, nmu, plot=True)
montecarlo.store(prefix + 'reflect', nphot, rdata_reflect, nmu, plot=True)
montecarlo.store(prefix, nphot, rdata_both, nmu, plot=True)




