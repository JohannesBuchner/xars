"""
Monte-Carlo simulator for X-ray obscurer geometries

Literature:

   * Brightman & Nandra (2011)
   * Leahy & Creighton (1993)
"""

import numpy
import scipy
from numpy import pi, arccos, tan, round, exp, log, log10, sin, cos, logical_and, logical_or, arctan
from binning import nbins, energy2bin, bin2energy
import progressbar
import matplotlib.pyplot as plt
import subprocess

energy_lo, energy_hi = bin2energy(range(nbins))
energy = (energy_hi + energy_lo)/2.
deltae = energy_hi - energy_lo

rng = scipy.random

from wedgetorus import WedgeTorusGeometry
from montecarlo import PhotonBunch

rng.seed(0)

nmu = 15 # number of viewing angle bins

import argparse
import sys, os

class WedgeTorusComponents(object): 
	""" Wedge Torus Segments 
	
	global: PhoIndex, log-normalization
	for each segment:
		log10 NH, r_inner, r_outer, Theta_low / pi, Theta_high / pi
	
	"""
	def __init__(self, n):
		self.n = n
		self.component_length = 5
		self.length = 2 + self.component_length * n
		self.low_bounds = [1.4, -6]
		self.high_bounds = [2.8, -4]
		for i in range(self.n):
			self.low_bounds  += [0, 0, 0, 0, 22]
			self.high_bounds += [1, 1, 1, 1, 26]
	def decode(self, encoding):
		assert len(encoding) == self.length, encoding
		d = dict(Theta_low = [], Theta_high = [],
			r_inner = [], r_outer = [],
			NH = [])
		
		max_r_outer = max([encoding[i+1] for i in range(2, len(encoding), self.component_length)])
		
		for i in range(2, len(encoding), self.component_length):
			d['r_inner'].append(max(encoding[i] / max_r_outer, 0.01))
			d['r_outer'].append(encoding[i+1] / max_r_outer)
			d['Theta_low'].append(encoding[i+2] * numpy.pi)
			d['Theta_high'].append(encoding[i+3] * numpy.pi)
			d['NH'].append(10**(encoding[i+4] - 22))
		parameters = dict(PhoIndex = encoding[0], lognormalization = encoding[1])
		return {'geometry':d, 'parameters':parameters}
	
	def generate(self, random, args, limits):
		while True:
			encoding = []
			for i in range(0, self.length):
				encoding.append(random.uniform(self.low_bounds[i], self.high_bounds[i]))
			
			# normalize
			max_r_outer = max([encoding[i+1] for i in range(2, len(encoding), self.component_length)])
			for i in range(2, len(encoding), self.component_length):
				encoding[i] /= max_r_outer
				encoding[i+1] /= max_r_outer
			
			d = self.decode(encoding)
			
			valid, reason = WedgeTorusGeometry.is_valid(**d['geometry'])
			if valid:
				break
			#else:
			#	print 'generated invalid geometry, will try again', d['geometry'], reason
		return encoding

	def visualize(self, encodings, limits):
		plt.figure()
		encoding = encodings[0]
		d = self.decode(encoding)
		geo = WedgeTorusGeometry(**d['geometry'])
		geo.plot = "wedgetorusviz.pdf"
		geo.viz()
		plt.savefig(geo.plot)
		plt.close()
		geo.plot = None
		return
		
	def parse_seeds(self, generator, seeds, random):
		# convert seeds which may not have length n to good candidates
		candidates = []
		print 'parsing based on %d seeds' % len(seeds)
		for s in seeds:
			if (len(s) - 2) % self.component_length != 0:
				print 'seed skipped because of wrong format', s
				continue
			if len(s) == self.length:
				candidates.append(s)
				continue
			if len(s) < self.length:
				for i in range(10): # fill up with a fresh candidate
					t = generator(random, [])
					for i in range(len(s)):
						t[i] = s[i]
					candidates.append(t)
				continue
			if len(s) > self.length:
				# compute number of choice combinations
				n = scipy.special.binom(((len(s) - 2) / self.component_length), self.n)
				for i in range(min(n, 10)):
					choices = random.sample(range((len(s) - 2) / self.component_length), self.n)
					t = s[:2]
					for j in choices:
						t += s[k+2:k+2+self.component_length]
					candidates.append(t)
		print 'parsing produced %d candidates' % len(seeds)
		return candidates

def run_simulation(geometry, prefix, nphot, accumulate = False, verbose = False):
	"""
	geometry: geometry to use
	prefix: filename prefix
	nphot: total number of photons to send in
	accumulate: if file exists and accumulate is True, the simulation will
		add the new computation onto the existing one. Otherwise, the 
		existing simulation will be returned immediately without any 
		computation.
	"""

	try:
		print 'Loading previous file if exists ...'
		import pyfits
		old_file = pyfits.open(os.environ.get('MAGICDIRECTORY', '') + prefix + "rdata.fits")
		rdata_old = old_file[0].data
		prev_nphot = old_file[0].header.get('NPHOT', 0)
		if not accumulate:
			return prev_nphot, rdata_old
	except Exception as e:
		print 'loading previous file failed; Error:', e
		print 'creating a fresh simulation'
		prev_nphot = 0
		rdata_old = 0
	
	rdata = numpy.zeros((nbins, nbins, nmu))
	rdata = [0] * nbins
	
	pbar = progressbar.ProgressBar(widgets=[
		progressbar.Percentage(), progressbar.Counter('%5d'), 
		progressbar.Bar(), progressbar.ETA()], maxval=nbins).start()

	for i in range(nbins):
		photons = PhotonBunch(i=i, nphot=nphot, verbose=verbose, geometry=geometry)
		remainder = [(photons.rad, photons.theta)]
		
		for n_interactions in range(1000):
			emission = photons.pump()
			if emission is None:
				break
			if len(emission['energy']) == 0:
				continue
			if verbose: print ' received %d emitted photons (after %d interactions)' % (len(emission['energy']), n_interactions)
			beta = emission['beta']
			beta = numpy.abs(beta)
			assert not numpy.any(beta > pi)
			assert not numpy.any(beta < 0)
			#bins = energy2bin(energy)
			# vertical bin, i.e. which viewing angle can see this photon?
			mbin = numpy.array(numpy.round((nmu - 1) * 0.5 * (1 - cos(beta))), dtype=numpy.uint)
			assert not numpy.any(mbin >= nmu)
			assert not numpy.any(mbin < 0)
			mbin[mbin == nmu] = nmu - 1
			bins = emission['bin']
		
			# produce unique array bins, mbin which contains counts
			counts, xedges, yedges = numpy.histogram2d(bins, mbin, bins=[range(nbins+1), range(nmu+1)])
			
			# record into histogram if it landed within relevant range
			rdata[i] += counts
		
		if i % 100 == 0 and verbose:
			plt.figure(figsize=(7,5))
			markers = ['x', '+', 's', 'o', '<', '>', '^', 'v']
		
			for marker, mu in zip(markers + markers[::-1], range(nmu)):
				plt.plot(energy, numpy.log10(1+rdata[i][:,mu]), ls='', marker=marker, label='viewing angle %d' % mu, alpha=0.3)
				#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
				#imshow(counts.T, extent=extent, interpolation='nearest', origin='lower')
			plt.legend(loc='best', prop=dict(size=6))
			plt.xlabel('energy')
			plt.ylabel('counts')
			plt.xlim(0.1, 20)
			plt.savefig(os.environ.get('MAGICDIRECTORY', '') + prefix + "abs_%d.pdf" % (i))
			plt.savefig(os.environ.get('MAGICDIRECTORY', '') + prefix + "abs_%d.png" % (i))
			plt.close()
			
		pbar.update(pbar.currval + 1)
	pbar.finish()
	
	rdata = numpy.array(rdata)
	plt.figure(figsize=(7,5))
	markers = ['x', '+', 's', 'o', '<', '>', '^', 'v']
	for marker, mu in zip(markers + markers[::-1], range(nmu)):
		plt.plot(energy, numpy.log10(1+rdata[:,:,mu].sum(axis=0)), ls='', marker=marker, label='viewing angle %d' % mu, alpha=0.3)
		#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
		#imshow(counts.T, extent=extent, interpolation='nearest', origin='lower')
	plt.legend(loc='best', prop=dict(size=6))
	plt.xlabel('energy')
	plt.ylabel('counts')
	plt.xlim(0.1, 12)
	plt.gca().set_xscale('log')
	plt.savefig(os.environ.get('MAGICDIRECTORY', '') + prefix + "abs.pdf")
	plt.savefig(os.environ.get('MAGICDIRECTORY', '') + prefix + "abs.png")
	plt.close()
	
	print 'Accumulating onto previous result ...'
	rdata_total = rdata + rdata_old
	nphot_total = nphot + prev_nphot

	# store in fits file now
	print 'storing ...'
	hdu = pyfits.PrimaryHDU(rdata_total)
	import datetime, time
	now = datetime.datetime.fromtimestamp(time.time())
	nowstr = now.isoformat()
	nowstr = nowstr[:nowstr.rfind('.')]
	hdu.header.update('CREATOR', """Johannes Buchner <jbuchner@mpe.mpg.de>""")
	hdu.header.update('DATE', nowstr)
	hdu.header.update('METHOD', 'Monte-Carlo simulation code')
	hdu.header.update('PREFIX', prefix)
	hdu.header.update('NPHOT', nphot_total)
	hdu.writeto(os.environ.get('MAGICDIRECTORY', '') + prefix + "rdata.fits", clobber=True)
	numpy.savez(os.environ.get('MAGICDIRECTORY', '') + prefix + "rdata.npz", matrix=rdata_total, x=energy, total=nphot_total, xwidth=deltae)

	print 'total of %d input / %d output photons across %d bins' % (nphot + prev_nphot, rdata.sum(), nbins)
	
	return nphot_total, rdata_total


from subprocess import Popen, PIPE, STDOUT
from StringIO import StringIO
import json
# print ' >>>', prefix, parameters
print 'calling external evaluator ...'
p = subprocess.Popen("./spectrumeval.sh", stdout=subprocess.PIPE, stdin=subprocess.PIPE)
print 'waiting for external evaluator to become ready...'
while True:
	line = p.stdout.readline()
	print ' << "%s"' % line.strip()
	if 'ready for input' in line:
		break

def spectrum_likelihood(prefix, parameters, nphot, rdata):
	# need to call external program
	#stdout, stderr = p.communicate("<< %s %s\n" % (prefix, [str(par) for par in parameters]))
	print '>>', "%s %s\n" % (prefix, " ".join([str(par) for par in parameters.values()]))
	p.stdin.write("%s %s\n" % (prefix, " ".join([str(par) for par in parameters.values()])))
	p.stdin.flush()
	
	while True:
		print 'reading input line ...',
		sys.stdout.flush()
		line = p.stdout.readline().strip()
		print ' << "%s"' % line
		if line.startswith('>>>'):
			obj = json.load(open(os.environ.get('MAGICDIRECTORY', '') + line[3:]))
			if 'error' in obj:
				raise Exception(obj['error'])
			limax = max(obj['stats'])
			limarg = log(exp(numpy.asarray(obj['stats']) - limax).sum()) + limax
			print 'evaluator computed', limarg
			return limarg, obj
	assert False, 'nothing to read!'

class SpectrumSimulationModel(object):
	def __init__(self, encoder, limits, nphot, accumulate, events=[], verbose = False):
		self.events = events
		self.encoder = encoder
		self.limits = limits
		self.accumulate = accumulate
		self.nphot = nphot
		self.verbose = verbose
	
	def likelihood(self, candidate):
		d = self.encoder.decode(candidate)
		geo = WedgeTorusGeometry(**d['geometry'])
		
		prefix = "wedgetorus_" + '_'.join(["%.1f" % c for c in candidate[2:]])
		fullprefix = prefix + "_" + '_'.join(["%.2f" % c for c in candidate[:2]])
		
		# plot geometry
		geo.plot = prefix + "viz.pdf"
		geo.viz()
		plt.savefig(os.environ.get('MAGICDIRECTORY', '') + geo.plot)
		plt.savefig(os.environ.get('MAGICDIRECTORY', '') + geo.plot.replace(".pdf", ".png"))
		plt.close()
		geo.plot = None
		
		# run simulation.
		nphot, rdata = run_simulation(geo, prefix, accumulate = self.accumulate, 
			nphot = self.nphot, verbose = self.verbose)
		# call evaluator
		margli, dataobj = spectrum_likelihood(prefix, d['parameters'], nphot, rdata)
		
		# plot spectra
		plt.figure()
		xlo, xhi, y = dataobj['data_binned']
		xmid  = (numpy.asarray(xlo) + numpy.asarray(xhi))/2.
		xdiff = (-numpy.asarray(xlo) + numpy.asarray(xhi))/2.
		plt.errorbar(x=xmid, xerr=xdiff, y=y, marker='+', color='black', ls='')
		for i, (li, (xlo, xhi, y)) in enumerate(zip(dataobj['stats'], dataobj['model'])):
			xmid  = (numpy.asarray(xlo) + numpy.asarray(xhi))/2.
			xdiff = (-numpy.asarray(xlo) + numpy.asarray(xhi))/2.
			plt.plot(xmid, y, '-', color='blue', 
				alpha=max(exp(li - margli), 0.2), lw=1, 
				label='bin %d, viewing angle: %.1f. stat=%.1f' % (
					i, 180 / pi * arccos(1 - i / 0.5 / (nmu - 1)), li))
		plt.xlabel('Energy [keV]')
		plt.ylabel('Flux [cts/cm^2/s]')
		plt.gca().set_xscale('log')
		plt.gca().set_yscale('log')
		plt.legend(loc='best', prop=dict(size=6))
		plt.savefig(os.environ.get('MAGICDIRECTORY', '') + fullprefix + "spec.pdf")
		plt.savefig(os.environ.get('MAGICDIRECTORY', '') + fullprefix + "spec.png")
		plt.close()
		
		plt.figure()
		for i, (li, (xlo, xhi, y)) in enumerate(zip(dataobj['stats'], dataobj['source'])):
			xmid  = (numpy.asarray(xlo) + numpy.asarray(xhi))/2.
			xdiff = (-numpy.asarray(xlo) + numpy.asarray(xhi))/2.
			plt.plot(xmid, y, '-', color='green', 
				alpha=max(exp(li - margli), 0.2), lw=1, 
				label='bin %d, viewing angle: %.1f. stat=%.1f' % (
					i, arccos(1 - i / 0.5 / (nmu - 1)), li))
		plt.xlabel('Energy [keV]')
		plt.ylabel('Flux [cts/cm^2/s]')
		plt.gca().set_xscale('log')
		plt.gca().set_yscale('log')
		plt.legend(loc='best', prop=dict(size=6))
		plt.savefig(os.environ.get('MAGICDIRECTORY', '') + fullprefix + "src.pdf")
		plt.savefig(os.environ.get('MAGICDIRECTORY', '') + fullprefix + "src.png")
		plt.close()
		return margli
	
	def generate_events(self, candidate, random, k):
		raise NotImplementedError

import evolve

evolve.prng.seed(0)
limits = [[0.1, 9.7]]

evolve.parser.add_argument('--nphot', type=int, default=100000, help='number of input photons per energy bin')
evolve.parser.add_argument('--accumulate', type=bool, default=False, help='continue working on existing spectrum')
evolve.parser.add_argument('--verbose', type=bool, default=False, help='be talkative')

evolve.parse_args()

encoder = WedgeTorusComponents(evolve.args.n_components)
ssm = SpectrumSimulationModel(encoder=encoder, limits=limits, 
	accumulate=evolve.args.accumulate,
	nphot=evolve.args.nphot,
	verbose = evolve.args.verbose
	)

ea, evolve_args = evolve.setup(ssm)
evolve.evolve(ea, evolve_args)
	

