import numpy
import scipy
from numpy import pi, arccos as acos, tan, round, log, log10, sin, cos, logical_and, logical_or, arctan as atan, arctan2 as atan2, exp
from binning import nbins, energy2bin, bin2energy
from xsect import xscatt, xphot, xkfe, xboth, absorption_ratio, fek_ratio, electmass
import progressbar
from photons import PhotonBunch

rng = scipy.random

def plot_interaction(nphot, n_interactions, rad, theta, beta, **kwargs):
	import matplotlib.pyplot as plt
	mask = rad < 1
	xi=rad*sin(theta)
	zi=rad*cos(theta)
	plt.figure(figsize=(5,5))
	plt.ylim(0,1)
	plt.xlim(0,1)
	plt.plot(xi[mask], zi[mask],   '.', ms=1, alpha=0.5, 
		label='%02.2f%% (%02.2f%% out) at interaction # %d' % (
		len(mask) * 100. / nphot, 
		mask.sum() * 100. / nphot,
		n_interactions ), color='red')
	plt.plot(xi[-mask], zi[-mask], '.', ms=1, alpha=0.5,
		label='%02.2f%% (%02.2f%% through) at interaction # %d' % (
		len(mask) * 100. / nphot, 
		(-mask).sum() * 100. / nphot,
		n_interactions ), color='blue')
	plt.vlines(xi[mask].mean(), 0, 1, linestyle='--', color='red')
	plt.vlines(xi[-mask].mean(), 0, 1, linestyle='--', color='blue')
	xv=sin(beta)
	zv=cos(beta)
	plt.plot(xv[mask],  zv[mask],  '+', color='red', ms=1)
	plt.plot(xv[-mask], zv[-mask], '+', color='blue', ms=1)
	
	plt.plot(sin(beta[mask].mean()),  cos(beta[mask].mean()),  'o', color='red')
	plt.plot(sin(beta[-mask].mean()), cos(beta[-mask].mean()), 'o', color='blue')
	plt.legend(loc='best', prop=dict(size=6))

def plot_path(rad_paths, theta_paths, **kwargs):
	import matplotlib.pyplot as plt
	for rad, theta in zip(rad_paths, theta_paths):
		xi=rad*sin(theta)
		zi=rad*cos(theta)
		print theta, rad, xi, zi
		plt.plot(xi, zi, 'o:', **kwargs)
	print

def run(prefix, nphot, nmu, geometry, 
	binmapfunction = lambda beta: numpy.floor(nmu * beta / pi),
	plot_paths = False, plot_interactions = False, verbose = False):
	
	if plot_paths or plot_interactions:
		import matplotlib.pyplot as plt
	
	rdata_transmit = numpy.zeros((nbins, nbins, nmu))
	rdata_reflect = numpy.zeros((nbins, nbins, nmu))
	#rdata = [0] * nbins
	energy_lo, energy_hi = bin2energy(range(nbins))
	energy = (energy_hi + energy_lo)/2.
	deltae = energy_hi - energy_lo
	
	pbar = progressbar.ProgressBar(widgets=[
		progressbar.Percentage(), progressbar.Counter('%5d'), 
		progressbar.Bar(), progressbar.ETA()], maxval=nbins).start()

	binrange = [range(nbins+1), range(nmu+1)]
	for i in range(nbins):
		photons = PhotonBunch(i=i, nphot=nphot, verbose=verbose, geometry=geometry)
		remainder = [(photons.rad, photons.theta)]
		if plot_paths:
			plt.figure("paths", figsize=(4, 4))
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
			#beta = numpy.abs(beta)
			#beta[beta > pi] = pi - beta[beta > pi]
			#beta[beta > pi/2.] = pi/2. - beta[beta > pi/2.]
			assert (beta <= pi).all(), beta
			assert (beta >= 0).all(), beta
			#bins = energy2bin(energy)
			# vertical bin, i.e. which viewing angle can see this photon?
			mbin = numpy.asarray(binmapfunction(beta=beta, alpha=alpha)).astype(numpy.uint)
			# highest bin exceeded due to rounding
			mbin[mbin == nmu] = nmu - 1
			
			bins = emission['bin']
			# produce unique array bins, mbin which contains counts
			counts, xedges, yedges = numpy.histogram2d(bins, mbin, bins=binrange)
			# record into histogram if it landed within relevant range
			if n_interactions < 1:
				rdata_transmit[i] += counts
			else:
				rdata_reflect[i] += counts
			#if (emission['energy'] == 6.40).any():
			#	print ' %f: %d/%d are lines' % (energy[i], (emission['energy'] == 7.06).sum(), (emission['energy'] == 6.40).sum())
			#	linebin = set(bins[emission['energy'] == 6.40])
			#	print linebin, rdata[i][724,:], rdata[900][724,4] if i > 900 else ''
			
			# remove the emitted photons from the remainder
			if plot_paths:
				mask = emission['mask']
				path = [(prev_rad[mask], prev_theta[mask]) for prev_rad, prev_theta in remainder]
				
				remainder = [(prev_rad[-mask], prev_theta[-mask]) 
					for prev_rad, prev_theta in remainder] + [(photons.rad, photons.theta)]
				
				rad_paths   = numpy.transpose([path_rad   for path_rad, path_theta in path] + [2*numpy.ones(mask.sum())])
				theta_paths = numpy.transpose([path_theta for path_rad, path_theta in path] + [beta])
				plt.figure("paths")
				plot_path(rad_paths[:100], theta_paths[:100], 
					color=(['b','r','g','y','k','m','w']*5)[n_interactions],
					alpha=1 - 0.75*numpy.exp(-n_interactions/5.))
			
			if plot_interactions:
				print 'plotting %d photons ...' % len(beta)
				plot_interaction(nphot=nphot, n_interactions=n_interactions, **emission)
				#plt.savefig(prefix + "rdata_%d_%d.pdf" % (i, n_interactions))
				plt.savefig(prefix + "rdata_%d_%d.png" % (i, n_interactions))
				plt.close()
				print 'plotting ... done'
			if not more:
				break
		if plot_paths:
			plt.figure("paths")
			plt.xlim(-1, 1)
			plt.ylim(-1, 1)
			plt.title('Energy: %.2f, %d interactions' % (energy[i], n_interactions))
			#plt.savefig(prefix + "paths_%d.pdf" % (i))
			plt.savefig(prefix + "paths_%d.png" % (i))
			plt.close()
			
		pbar.update(pbar.currval + 1)
	pbar.finish()

	return (rdata_transmit, rdata_reflect), nphot

def store(prefix, nphot, rdata, nmu, extra_fits_header = {}, plot=False):
	
	energy_lo, energy_hi = bin2energy(range(nbins))
	energy = (energy_hi + energy_lo)/2.
	deltae = energy_hi - energy_lo
	try:
		import pyfits
	except ImportError:
		import astropy.io.fits as pyfits
	try:
		print 'Loading previous file if exists ...'
		old_file = pyfits.open(prefix + "rdata.fits")
		rdata_old = old_file[0].data
		prev_nphot = old_file[0].header.get('NPHOT', 0)
		print 'Accumulating onto previous result ...'
		rdata = rdata + rdata_old
	except Exception as e:
		print 'updating file failed; writing fresh. Error:', e
		prev_nphot = 0
	nphot_total = nphot + prev_nphot
	print 'storing ...'
	hdu = pyfits.PrimaryHDU(rdata)
	import datetime, time
	now = datetime.datetime.fromtimestamp(time.time())
	nowstr = now.isoformat()
	nowstr = nowstr[:nowstr.rfind('.')]
	hdu.header['CREATOR'] = """Johannes Buchner <johannes.buchner.acad@gmx.com>"""
	hdu.header['DATE'] = nowstr
	hdu.header['METHOD'] = 'Monte-Carlo simulation code'
	hdu.header['NPHOT'] = nphot_total
	for k,v in extra_fits_header.iteritems():
		hdu.header[k] = v
	hdu.writeto(prefix + "rdata.fits", clobber=True)
	print 'total of %d input / %d output photons across %d bins' % (nphot_total, rdata.sum(), nbins)
	
	if not plot:
		return nphot_total, rdata
	import matplotlib.pyplot as plt
	PhoIndex = 2
	matrix = rdata
	print matrix[600,500:600,3]
	print matrix[500,100:500,3]
	#for i in range(len(energy)):
	#	eselect = energy.reshape((-1,1)) * (matrix[i] > 0)
	#	eselect = eselect[eselect > 0]
	#	if len(eselect) == 0: continue
	#	print '%.2f %.2f' % (energy[i], eselect.min())
	#print energy.reshape((1,-1,1)) * (matrix > 0)
	#print numpy.argmax(energy.reshape((1,-1,1)) * (matrix > 0), axis=0)
	x = energy
	total = nphot_total
	xwidth = deltae
	print 'plotting...'
	NH = 1e24/1e22
	weights = (energy**-PhoIndex * deltae).reshape((-1,1))
	yall = (weights * matrix.sum(axis=2)).sum(axis=0) / deltae
	for mu in range(nmu):
		y = (weights * matrix[:,:,mu]).sum(axis=0) / deltae
		print '%d ... ' % mu
		
		plt.figure(figsize=(10,10))
		plt.plot(energy, exp(-xphot*NH) * energy**-PhoIndex / nmu, '-', color='red', linewidth=1)
		plt.plot(energy, exp(-xscatt*NH) * energy**-PhoIndex / nmu, '-', color='pink')
		plt.plot(energy, exp(-xkfe*NH) * energy**-PhoIndex / nmu, '-', color='orange')
		plt.plot(energy, energy**-PhoIndex / nmu, '--', color='gray')
		plt.plot(energy_lo, y / total, '-', color='k') #, drawstyle='steps')
		plt.plot(energy_lo, yall / total / nmu, '-', color='gray', alpha=0.3, linewidth=3) #, drawstyle='steps')
		#plt.plot(energy, exp(-xboth) * energy**-PhoIndex, '-', color='yellow')
		plt.gca().set_xscale('log')
		plt.gca().set_yscale('log')
		#plt.xlim(0.1, 10 * (1 + 10))
		plt.xlim(3, 40)
		lo, hi = 1e-8, 1
		plt.vlines(6.40, lo, hi, linestyles=[':'], color='grey', alpha=0.5)
		plt.vlines(7.06, lo, hi, linestyles=[':'], color='grey', alpha=0.5)
		plt.ylim(lo, hi)
		plt.show()
		plt.savefig(prefix + "_%d.pdf" % mu)
		plt.savefig(prefix + "_%d.png" % mu)
		numpy.savetxt(prefix + "_%d.txt" % mu, numpy.vstack([energy, y]).transpose())
		plt.close()
	return nphot_total, rdata


