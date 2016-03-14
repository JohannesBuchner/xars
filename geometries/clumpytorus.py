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
import copy
import matplotlib.pyplot as plt

energy_lo, energy_hi = bin2energy(range(nbins))
energy = (energy_hi + energy_lo)/2.
deltae = energy_hi - energy_lo

rng = scipy.random

from wedgetorus import WedgeTorusGeometry
import montecarlo

rng.seed(0)

import argparse
import sys, os

parser = argparse.ArgumentParser(
	description="""Monte-Carlo simulator for X-ray obscurer geometries""",
	epilog="""(C) Johannes Buchner, 2013. Based on work by Murray Brightman & Kirpal Nandra (see 2011 publication)""")

parser.add_argument('--log10nh', type=float, required=True, help='column density (10^X cm^-1)')
parser.add_argument('--blob-size', type=float, required=True, help='length of blobs compared to size of torus region (0-1)')
parser.add_argument('--filling', type=float, required=True, help='filling of volume by blobs (0-1)')
parser.add_argument('--verbose', '-v', default=False, help='Be more talkative, show debug statistics on interactions', action='store_true')
parser.add_argument('--nevents', type=int, default=100000, help='number of input photons per energy bin')
args = parser.parse_args()

nmu = 10     # number of viewing angle bins

nh_in = args.log10nh
nh = 10**(nh_in-22)   # current NH value
print '  NH                : 10^%.1f' % (nh_in)
print '  blob size         : %.1f%%' % (args.blob_size * 100.)
assert 0 < args.blob_size < 1
print '  filling           : %.1f%%' % (args.filling * 100.)
assert 0 < args.filling < 1
prefix = 'clumpytorus_%.1f_%.1f_%.1f_' % (nh_in, args.blob_size * 100, args.filling * 100.)

rdata = [0] * nbins
nphot = args.nevents  # total number of photons to send in

# create geometry
d = dict(Theta_low = [], Theta_high = [],
	r_inner = [], r_outer = [],
	NH = [])

volume = 0.

while True: # volume < 4 * pi / 3 * args.filling
	d2 = copy.deepcopy(d)
	r = rng.uniform(0,1-args.blob_size)
	r_outer = r + args.blob_size
	Theta = rng.uniform(0,1-args.blob_size)
	Theta_outer = Theta + args.blob_size
	d2['r_inner'].append(r)
	d2['r_outer'].append(r_outer)
	d2['Theta_low'].append(Theta * numpy.pi)
	d2['Theta_high'].append(Theta_outer * numpy.pi)
	d2['NH'].append(nh)
	
	valid, reason = WedgeTorusGeometry.is_valid(**d2)
	if valid:
		d = d2
		#volume += (r_outer**2 - r**2) / 2 * (args.blob_size * pi - sin(args.blob_size * pi))
		n_components = len(d['r_outer'])
		volume = n_components * args.blob_size * args.blob_size
		print 'filled %.1f%%' % (volume * 100.), len(d['r_outer']), 'using', r, Theta
		if volume > args.filling:
			break

print
print 'using %d wedges' % n_components
geo = WedgeTorusGeometry(**d)
geo.plot = prefix + "_viz.pdf"
geo.viz()
plt.savefig(geo.plot)
plt.close()
geo.plot = None

header = dict(NH='%f' % nh, LOGNH='%f' % nh_in, BLOBSIZE='%f' % args.blob_size, FILLING='%f' % args.filling)
montecarlo.run(prefix, nbins, nphot, nmu, geometry=geo, 
	binmapfunction = lambda beta: numpy.floor(nmu * beta / pi), 
	extra_fits_header = header)


"""
pbar = progressbar.ProgressBar(widgets=[
	progressbar.Percentage(), progressbar.Counter('%5d'), 
	progressbar.Bar(), progressbar.ETA()], maxval=nbins).start()

for i in range(nbins):
	photons = PhotonBunch(i=i, nphot=nphot, verbose=args.verbose, geometry=geo)
	remainder = [(photons.rad, photons.theta)]
	
	for n_interactions in range(1000):
		emission = photons.pump()
		if emission is None:
			break
		if len(emission['energy']) == 0:
			continue
		if args.verbose: print ' received %d emitted photons (after %d interactions)' % (len(emission['energy']), n_interactions)
		beta = emission['beta']
		beta[beta > pi/2.] = pi/2. - beta[beta > pi/2.]
		beta = numpy.abs(beta)
		#bins = energy2bin(energy)
		# vertical bin, i.e. which viewing angle can see this photon?
		mbin = numpy.array(numpy.floor(nmu * beta / pi), dtype=numpy.uint)
		# highest bin exceeded due to rounding
		mbin[mbin == nmu] = nmu - 1
		
		bins = emission['bin']
		# produce unique array bins, mbin which contains counts
		counts, xedges, yedges = numpy.histogram2d(bins, mbin, bins=[range(nbins+1), range(nmu+1)])
		# record into histogram if it landed within relevant range
		rdata[i] += counts
	pbar.update(pbar.currval + 1)
pbar.finish()

rdata = numpy.array(rdata)

try:
	print 'Loading previous file if exists ...'
	import pyfits
	old_file = pyfits.open(prefix + "rdata.fits")
	rdata_old = old_file[0].data
	prev_nphot = old_file[0].header.get('NPHOT', 0)
	print 'Accumulating onto previous result ...'
	rdata = rdata + rdata_old
except Exception as e:
	print 'updating file failed; writing fresh. Error:', e
	prev_nphot = 0

print 'storing ...'
hdu = pyfits.PrimaryHDU(rdata)
import datetime, time
now = datetime.datetime.fromtimestamp(time.time())
nowstr = now.isoformat()
nowstr = nowstr[:nowstr.rfind('.')]
hdu.header.update('CREATOR', """Johannes Buchner <jbuchner@mpe.mpg.de>""")
hdu.header.update('DATE', nowstr)
hdu.header.update('METHOD', 'Monte-Carlo simulation code')
hdu.header.update('NPHOT', nphot + prev_nphot)
hdu.header.update('NH', '%f' % nh)
hdu.header.update('LOGNH', '%f' % nh_in)
hdu.header.update('BLOBSIZE', '%f' % args.blob_size)
hdu.header.update('FILLING', '%f' % args.filling)
hdu.writeto(prefix + "rdata.fits", clobber=True)

print 'total of %d input / %d output photons across %d bins' % (nphot + prev_nphot, rdata.sum(), nbins)

PhoIndex = 1.9
matrix = rdata
x = energy
total = nphot + prev_nphot
xwidth = deltae
print 'plotting...'
for mu in range(nmu):
	a = energy_lo
	b = energy_hi
	G = PhoIndex
	weights = energy**-PhoIndex * deltae
	y = numpy.multiply(weights, matrix[:,:,mu].transpose()).sum(axis=0) / deltae
	print '%d ... ' % mu
	
	plt.figure(figsize=(7,7))
	plt.plot(energy_lo, y, '+-', drawstyle='steps')
	plt.plot(energy, energy**-PhoIndex, '--')
	plt.gca().set_xscale('log')
	plt.gca().set_yscale('log')
	plt.xlim(0.1, 50)
	plt.savefig(prefix + "_%d.pdf" % mu)
	plt.savefig(prefix + "_%d.png" % mu)
	numpy.savetxt(prefix + "_%d.txt" % mu, numpy.vstack([energy, y]).transpose())
	plt.close()



"""
