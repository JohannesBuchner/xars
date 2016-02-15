import numpy
import scipy
from numpy import pi, arccos as acos, tan, round, log, log10, sin, cos, logical_and, logical_or, arctan as atan
from binning import nbins, energy2bin, bin2energy
from xsect import xscatt, xphot, xkfe, xboth, absorption_ratio, fek_ratio, electmass
import progressbar

rng = scipy.random

class PhotonBunch(object):
	def __init__(self, i, # input energy bin
		nphot, # number of photons
		geometry,
		verbose = False):
		self.verbose = verbose
		self.rad = numpy.zeros(nphot)
		self.phi = numpy.zeros(nphot) # initial left-right angle for position
		self.theta = numpy.zeros(nphot) # initial up-down angle for position
		self.alpha = numpy.zeros(nphot) # initial left-right direction for direction -- value does not matter due to symmetry
		mu = numpy.linspace(-1, 1, nphot) # uniform distribution
		self.beta = acos(mu) # up-down angle for direction
		# self.beta = acos(numpy.linspace(-cone_in, cone_in, nphot)) # up-down angle for direction
		self.geometry = geometry
		energy_lo, energy_hi = bin2energy(i)
		e = (energy_lo + energy_hi) / 2.
		if self.verbose: print 'PhotonBunch of size %d with energy %.2f keV' % (nphot, e)
		self.energy = e * numpy.ones(nphot)
		self.bin = i * numpy.ones(nphot, dtype=numpy.uint)
	
	def cut(self, mask):
		self.phi, self.theta, self.rad = self.phi[mask], self.theta[mask], self.rad[mask]
		self.alpha, self.beta = self.alpha[mask], self.beta[mask]
		self.energy, self.bin = self.energy[mask], self.bin[mask]
	
	def set(self, phi, theta, rad, alpha, beta, energy, bin):
		self.phi, self.theta, self.rad, self.alpha, self.beta, self.energy, self.bin = phi, theta, rad, alpha, beta, energy, bin
	
	def get(self):
		return self.phi, self.theta, self.rad, self.alpha, self.beta, self.energy, self.bin
	
	def pump(self):
		phi, theta, rad, alpha, beta, energy, bin = self.get()
		
		if self.verbose: print 'photon iteration of %d photons %s' % (len(energy), '_'*20)
		if len(energy) == 0:
			return None
		
		if self.verbose: print '  .. distance travelled'
		nphot = len(energy)
		r1 = rng.uniform(size=nphot)
		taur = - log(1. - r1) # effective tau
		dist = taur / xboth[bin] # distance travelled
		#if self.verbose: print '  .. .. tau min: %f mean: %f max: %f ' % (taur.min(), taur.mean(), taur.max())
		if self.verbose: print '  .. .. dist min: %f mean: %f max: %f ' % (dist.min(), dist.mean(), dist.max())
		
		phi0 = phi
		theta0 = theta
		rad0 = rad
	  	
	  	#if self.verbose: print '  .. computing position'
		# compute position
		xi=rad0*sin(theta0)*cos(phi0)
		yi=rad0*sin(theta0)*sin(phi0)
		zi=rad0*cos(theta0)
		
		inside, (xf,yf,zf), (rad, phi, theta) = self.geometry.compute_next_point((xi, yi, zi), (dist, beta, alpha))
		outside = -inside

		# emit
		if self.verbose: print '  .. emitting %d, %d left' % ((-inside).sum(), inside.sum())
		self.set(phi, theta, rad, alpha, beta, energy, bin)
		self.cut(inside)
		emit = dict(phi=phi0[outside], theta=theta0[outside], rad=rad0[outside], 
			beta=beta[outside], energy=energy[outside], bin=bin[outside], mask=outside)
		#print '   ', rad.shape, theta.shape, len(inside)
		xf, yf, zf = xf[inside], yf[inside], zf[inside]
		phi, theta, rad, alpha, beta, energy, bin = self.get()
		nphot = len(energy)
		
	  	if self.verbose: print '  .. checking if absorbed'
		# what effect are we dealing with
		r2 = rng.uniform(size=nphot)
		
		# Photon-absorbed
		q = absorption_ratio[bin]
		photabsorbed = r2 < q
		#if self.verbose: print '  .. .. photon-abs prob min: %f mean: %f max: %f ' % (q.min(), q.mean(), q.max())
		
		# Iron to photon-absorption effectiveness
		omega = fek_ratio[bin[photabsorbed]] # Importance of lines compared to photon-absorption
		#print '  ..  .. omega:', omega
		r3 = rng.uniform(size=photabsorbed.sum())
		
		# Are we coming out as a line
		is_line = r3 < omega
		photabsorbed_line = photabsorbed.copy()
		photabsorbed_notline = photabsorbed.copy()
		# set the absorbed ones (where we have true) to the criterion
		photabsorbed_line[photabsorbed] = is_line
		photabsorbed_notline[photabsorbed] = -is_line
		
		r4 = rng.uniform(size=photabsorbed_line.sum())
		
		is_Kbeta = r4 > 0.866
		is_lineKbeta = photabsorbed_line.copy()
		is_lineKalpha = photabsorbed_line.copy()
		is_lineKbeta[photabsorbed_line] = is_Kbeta
		is_lineKalpha[photabsorbed_line] = -is_Kbeta
		
		# emit Fe-Kbeta line  # 7.0580
		#e = energy2bin(energy[photabsorbed_line][is_lineKbeta])
		if is_lineKbeta.any():
			energy[is_lineKbeta] = 7.06
			e = energy2bin(energy[is_lineKbeta])
			bin[is_lineKbeta] = e
			e = None
		# emit Fe-Kalpha line # 6.4038 & 6.3908
		if is_lineKalpha.any():
			energy[is_lineKalpha] = 6.40
			e = energy2bin(energy[is_lineKalpha])
			bin[is_lineKalpha] = e
			e = None
		# no, we are just being swallowed (photon-absorption).
		# photabsorbed & -is_line
		# remove photabsorbed_line 
		
		absorbed = photabsorbed_notline
		if self.verbose: print '  .. .. line emission: %3d (%3d Kalpha, %3d Kbeta)' % (is_line.sum(), (-is_Kbeta).sum(), is_Kbeta.sum())
		if self.verbose: print '  .. .. absorbed: %d (%d lost)' % (photabsorbed.sum(), absorbed.sum())
		
	  	if self.verbose: print '  .. checking if scattered'
		# Were we compton-scattered?
		photscattered = (-photabsorbed)
		nphotscattered = photscattered.sum()
		#assert nphotscattered.shape == energy.shape
		if nphotscattered > 0:
			if self.verbose: print '  .. .. scattered: %d' % nphotscattered
			a = rng.uniform(size=nphotscattered)
			# compute new direction:
			alpha[photscattered] = a * 2. * pi # left-right angle uniform randomly

			# compute up-down angle
			beta0 = beta[photscattered]
			alpha0 = alpha[photscattered]
			r5 = rng.uniform(size=nphotscattered)
			r5a = rng.uniform(size=nphotscattered)
			
			x = 2. * r5a - 1
			mus = numpy.where(r5 > 0.75, 
				numpy.sign(x) * numpy.abs(x)**(1/3.), 
				x)
			betas = acos(mus)
			mu = cos(beta0)*cos(betas)+sin(beta0)*sin(betas)*cos(alpha0)
			beta[photscattered] = acos(mu)
			
			# new energy
			#if self.verbose: print '  .. .. mus: %.2f' % (mus.mean())
			loss = (1. + (1. - mus) * energy[photscattered] / electmass)
			if self.verbose: print '  .. .. energy loss: %.3f, %.3f, %.3f' % (
				loss.min(), loss.mean(), loss.max())
			energy[photscattered] /= loss
			if self.verbose: print '  .. .. new energy: %.3f, %.3f, %.3f' % (
				energy[photscattered].min(), energy[photscattered].mean(), energy[photscattered].max())
			bin[photscattered] = energy2bin(energy[photscattered])
		
	  	if self.verbose: print '  .. checking if outside of energy range'
		# if in relevant energy range, find right bin
		energy_outside = numpy.logical_or(energy < 0.1, energy > 1100)
		dropouts = numpy.logical_or(absorbed, energy_outside)
		remainders = -dropouts
		if self.verbose: print '  .. .. outside of energy range: %d' % (energy_outside.sum())
		if self.verbose: print '  .. .. %d left' % (remainders.sum())
		
		# cut to remainders (so that we only have to compute new positions for few)
		self.set(phi, theta, rad, alpha, beta, energy, bin)
		self.cut(remainders)
		phi, theta, rad, alpha, beta, energy, bin = self.get()
		
		# compute new positions for remainder
		xf, yf, zf = xf[remainders], yf[remainders], zf[remainders]
		rad = (xf**2+yf**2+zf**2)**0.5
		phi = numpy.zeros_like(xf)
		phi_mask = xf != 0
		phi[phi_mask] = atan(yf[phi_mask] / xf[phi_mask])
		theta = numpy.where(rad == 0, 0., acos(zf / rad))
		self.set(phi, theta, rad, alpha, beta, energy, bin)
		
		# next round
		return emit

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
	
	#rdata = numpy.zeros((nbins, nbins, nmu))
	rdata = [0] * nbins
	energy_lo, energy_hi = bin2energy(range(nbins))
	energy = (energy_hi + energy_lo)/2.
	deltae = energy_hi - energy_lo
	
	pbar = progressbar.ProgressBar(widgets=[
		progressbar.Percentage(), progressbar.Counter('%5d'), 
		progressbar.Bar(), progressbar.ETA()], maxval=nbins).start()

	binrange = numpy.array([range(nbins+1), range(nmu+1)])
	for i in range(nbins):
		photons = PhotonBunch(i=i, nphot=nphot, verbose=verbose, geometry=geometry)
		remainder = [(photons.rad, photons.theta)]
		if plot_paths:
			plt.figure("paths", figsize=(4, 4))
		for n_interactions in range(1000):
			emission = photons.pump()
			if emission is None:
				break
			if len(emission['energy']) == 0:
				continue
			if verbose: print ' received %d emitted photons (after %d interactions)' % (len(emission['energy']), n_interactions)
			beta = emission['beta']
			#beta = numpy.abs(beta)
			#beta[beta > pi] = pi - beta[beta > pi]
			#beta[beta > pi/2.] = pi/2. - beta[beta > pi/2.]
			assert (beta <= pi).all(), beta
			assert (beta >= 0).all(), beta
			#bins = energy2bin(energy)
			# vertical bin, i.e. which viewing angle can see this photon?
			mbin = numpy.asarray(binmapfunction(beta)).astype(numpy.uint)
			# highest bin exceeded due to rounding
			mbin[mbin == nmu] = nmu - 1
			
			bins = emission['bin']
			# produce unique array bins, mbin which contains counts
			counts, xedges, yedges = numpy.histogram2d(bins, mbin, bins=binrange)
			# record into histogram if it landed within relevant range
			rdata[i] += counts
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

	rdata = numpy.array(rdata)
	return rdata, nphot

def store(prefix, nphot, rdata, nmu, extra_fits_header = {}):
	
	energy_lo, energy_hi = bin2energy(range(nbins))
	energy = (energy_hi + energy_lo)/2.
	deltae = energy_hi - energy_lo
	
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
	nphot_total = nphot + prev_nphot
	print 'storing ...'
	hdu = pyfits.PrimaryHDU(rdata)
	import datetime, time
	now = datetime.datetime.fromtimestamp(time.time())
	nowstr = now.isoformat()
	nowstr = nowstr[:nowstr.rfind('.')]
	hdu.header.update('CREATOR', """Johannes Buchner <jbuchner@mpe.mpg.de>""")
	hdu.header.update('DATE', nowstr)
	hdu.header.update('METHOD', 'Monte-Carlo simulation code')
	hdu.header.update('NPHOT', nphot_total)
	for k,v in extra_fits_header.iteritems():
		hdu.header.update(k, v)
	hdu.writeto(prefix + "rdata.fits", clobber=True)
	print 'total of %d input / %d output photons across %d bins' % (nphot_total, rdata.sum(), nbins)

	import matplotlib.pyplot as plt
	PhoIndex = 1.9
	matrix = rdata
	print matrix[600,500:600,3]
	print matrix[500,100:500,3]
	x = energy
	total = nphot_total
	xwidth = deltae
	print 'plotting...'
	for mu in range(nmu):
		a = energy_lo
		b = energy_hi
		G = PhoIndex
		weights = energy**-PhoIndex * deltae
		# .transpose()
		y = (weights * matrix[:,:,mu]).sum(axis=0) / deltae
		print '%d ... ' % mu
		
		plt.figure(figsize=(7,7))
		plt.plot(energy_lo, y / total, '+ ') #, drawstyle='steps')
		plt.plot(energy, energy**-PhoIndex, '--')
		plt.gca().set_xscale('log')
		plt.gca().set_yscale('log')
		plt.xlim(0.1, 10 * (1 + 6))
		lo, hi = plt.ylim()
		plt.vlines(6.40, lo, hi, linestyles=[':'], color='grey', alpha=0.5)
		plt.vlines(7.06, lo, hi, linestyles=[':'], color='grey', alpha=0.5)
		plt.ylim(lo, hi)
		plt.savefig(prefix + "_%d.pdf" % mu)
		plt.savefig(prefix + "_%d.png" % mu)
		numpy.savetxt(prefix + "_%d.txt" % mu, numpy.vstack([energy, y]).transpose())
		plt.close()
	return nphot_total, rdata


