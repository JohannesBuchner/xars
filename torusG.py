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
from binning import nbins, energy2bin, bin2energy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

rng = scipy.random

from geometries.hydrotorus import HydroTorusGeometry
import montecarlo
from photons import PhotonBunch

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

nmu = 3     # number of viewing angle bins
nh_bins = numpy.array([9.99999978e-03, 1.41000003e-02, 1.99999996e-02,
	2.82000005e-02,   3.97999994e-02,   5.62000014e-02,
	7.94000030e-02,   1.12000003e-01,   1.58000007e-01,
	2.24000007e-01,   3.16000015e-01,   4.46999997e-01,
	6.30999982e-01,   8.90999973e-01,   1.25999999e+00,
	1.77999997e+00,   2.50999999e+00,   3.54999995e+00,
	5.01000023e+00,   7.07999992e+00,   1.00000000e+01,
	1.41000004e+01,   2.00000000e+01,   2.82000008e+01,
	3.97999992e+01,   5.62000008e+01,   7.94000015e+01,
	1.12000000e+02,   1.58000000e+02,   2.24000000e+02,
	3.16000000e+02,   4.47000000e+02,   6.31000000e+02,
	8.91000000e+02,   1.26000000e+03,   1.78000000e+03,
	2.51000000e+03,   3.55000000e+03,   5.01000000e+03,
	7.08000000e+03,   1.00000000e+04])

n_nh_bins = len(nh_bins)
binmapfunction = lambda beta, alpha: (numpy.round(0.5 + nmu * numpy.abs(cos(beta))) - 1).astype(int)

prefix = args.geometry + '_out'

geometry = HydroTorusGeometry(args.geometry, verbose=args.verbose)
geometry.viz()
plt.savefig(prefix + "geometry.pdf")
plt.savefig(prefix + "geometry.png")
plt.close()

def compute_normalisation(prefix, binmapfunction, verbose=False, nphot=1000000):
	# use 40000 rays in random directions
	import astropy.io.fits as pyfits
	if verbose:
		print 'computing normalisation ...'
	photons = PhotonBunch(i=100, nphot=nphot, verbose=verbose, geometry=geometry)
	# vertical bin, i.e. which viewing angle can see this photon?
	mbin = numpy.asarray(binmapfunction(beta=photons.beta, alpha=photons.alpha)).astype(numpy.uint)
	# highest bin exceeded due to rounding
	mbin[mbin == nmu] = nmu - 1

	# bin in NH
	if verbose:
		print '   computing LOS NH ...'
	nh = geometry.compute_los_nh(photons.beta, photons.alpha)
	if verbose:
		print '   computing LOS NH ... done'
	nh[nh<1e-2] = 1e-2
	kbin = ((log10(nh) + 2) * n_nh_bins / (4 + 2)).astype(int)
	kbin[kbin == n_nh_bins] = n_nh_bins - 1

	mkbin = kbin * nmu + mbin

	# compute fraction in the given NH/mu bins
	counts, xedges = numpy.histogram(mkbin, bins=range(nmu*n_nh_bins+1))
	normalisation = counts * 1. / len(mkbin)
	print normalisation, normalisation.shape
	hdu = pyfits.PrimaryHDU(normalisation)
	import datetime, time
	now = datetime.datetime.fromtimestamp(time.time())
	nowstr = now.isoformat()
	nowstr = nowstr[:nowstr.rfind('.')]
	hdu.header['CREATOR'] = """Johannes Buchner <johannes.buchner.acad@gmx.com>"""
	hdu.header['DATE'] = nowstr
	hdu.header['METHOD'] = 'Monte-Carlo simulation code'
	hdu.header['NPHOT'] = nphot
	if verbose:
		print '   saving ...'
	hdu.writeto(prefix + "normalisation.fits", clobber=True)
	if verbose:
		print '   saving ... done'
	return normalisation

if not os.path.exists(prefix + "normalisation.fits"):
	compute_normalisation(prefix, binmapfunction=binmapfunction, verbose=True)

def run(prefix, nphot, nmu, n_nh_bins, geometry, binmapfunction, verbose=False):
	rdata_transmit = numpy.zeros((nbins, nbins, nmu*n_nh_bins))
	rdata_reflect = numpy.zeros((nbins, nbins, nmu*n_nh_bins))
	#rdata = [0] * nbins
	energy_lo, energy_hi = bin2energy(range(nbins))
	energy = (energy_hi + energy_lo)/2.
	deltae = energy_hi - energy_lo
	
	pbar = progressbar.ProgressBar(widgets=[
		progressbar.Percentage(), progressbar.Counter('%5d'), 
		progressbar.Bar(), progressbar.ETA()], maxval=nbins).start()

	binrange = [range(nbins+1), range(nmu*n_nh_bins+1)]
	for i in range(nbins):
		photons = PhotonBunch(i=i, nphot=nphot, verbose=verbose, geometry=geometry)
		for n_interactions in range(1000):
			emission, more = photons.pump()
			if emission is None and not more:
				break
			if emission is None:
				continue
			if len(emission['energy']) == 0:
				if not more:
					break
				continue
			if verbose: print ' received %d emitted photons (after %d interactions)' % (len(emission['energy']), n_interactions)
			beta = emission['beta']
			alpha = emission['alpha']
			assert (beta <= pi).all(), beta
			assert (beta >= 0).all(), beta
		
			# vertical bin, i.e. which viewing angle can see this photon?
			mbin = numpy.asarray(binmapfunction(beta=beta, alpha=alpha)).astype(numpy.uint)
			# highest bin exceeded due to rounding
			mbin[mbin == nmu] = nmu - 1
		
			# bin in NH
			nh = geometry.compute_los_nh(beta, alpha)
			nh[nh<1e-2] = 1e-2
			kbin = ((log10(nh) + 2) * n_nh_bins / (4 + 2)).astype(int)
			kbin[kbin == n_nh_bins] = n_nh_bins - 1
		
			mkbin = kbin * nmu + mbin
		
			bins = emission['bin']
			# produce unique array bins, mbin which contains counts
			counts, xedges, yedges = numpy.histogram2d(bins, mkbin, bins=binrange)
			# record into histogram if it landed within relevant range
			if n_interactions < 1:
				rdata_transmit[i] += counts
			else:
				rdata_reflect[i] += counts
			del counts, emission, bins
			if not more:
				break
		del photons
		pbar.update(pbar.currval + 1)
	pbar.finish()

	return (rdata_transmit, rdata_reflect), nphot

(rdata_transmit, rdata_reflect), nphot = run(prefix, nphot = args.nevents, nmu = nmu, n_nh_bins = n_nh_bins, geometry=geometry, 
	binmapfunction = binmapfunction, verbose=args.verbose)

montecarlo.store(prefix + 'transmit', nphot, rdata_transmit, nmu*n_nh_bins, plot=False)
montecarlo.store(prefix + 'reflect', nphot, rdata_reflect, nmu*n_nh_bins, plot=False)
rdata_both = rdata_transmit + rdata_reflect
del rdata_transmit, rdata_reflect
montecarlo.store(prefix, nphot, rdata_both, nmu*n_nh_bins, plot=True)

