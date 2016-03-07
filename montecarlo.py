import numpy
import scipy
from numpy import pi, arccos as acos, tan, round, log, log10, sin, cos, logical_and, logical_or, arctan as atan, arctan2 as atan2, exp
from coordtrans import to_spherical, to_cartesian
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
		#mu = numpy.linspace(-1, 1, nphot) # uniform distribution
		mu = numpy.random.uniform(-1, 1, size=nphot)
		self.beta = acos(mu) # up-down angle for direction
		# self.beta = acos(numpy.linspace(-cone_in, cone_in, nphot)) # up-down angle for direction
		self.geometry = geometry
		energy_lo, energy_hi = bin2energy(i)
		e = (energy_lo + energy_hi) / 2.
		if self.verbose: print 'PhotonBunch of size %d with energy %.2f keV' % (nphot, e)
		self.energy = e * numpy.ones(nphot)
		self.bin = i * numpy.ones(nphot, dtype=numpy.uint)
		self.stuck = self.rad != 0 # False
	
	def cut(self, mask):
		# TODO: we could also refill with new photons
		self.phi, self.theta, self.rad = self.phi[mask], self.theta[mask], self.rad[mask]
		self.alpha, self.beta = self.alpha[mask], self.beta[mask]
		self.energy, self.bin = self.energy[mask], self.bin[mask]
		self.stuck = self.stuck[mask]
	
	def cut_free(self, free_mask):
		free = -self.stuck
		free[free] = free_mask
		self.cut(free)
	
	def set(self, phi, theta, rad, alpha, beta, energy, bin):
		self.phi, self.theta, self.rad, self.alpha, self.beta, self.energy, self.bin = phi, theta, rad, alpha, beta, energy, bin
	
	def get(self):
		return self.phi, self.theta, self.rad, self.alpha, self.beta, self.energy, self.bin
	
	def get_free(self):
		free = -self.stuck
		return self.phi[free], self.theta[free], self.rad[free], self.alpha[free], self.beta[free], self.energy[free], self.bin[free]

	def update_free(self, phi, theta, rad, alpha, beta, energy, bin):
		free = -self.stuck
		self.phi[free], self.theta[free], self.rad[free], self.alpha[free], self.beta[free], self.energy[free], self.bin[free] =  phi, theta, rad, alpha, beta, energy, bin
	
	def get_stuck(self):
		return self.alpha[self.stuck], self.beta[self.stuck], self.energy[self.stuck], self.bin[self.stuck]

	def update_and_free_stuck(self, alpha, beta, energy, bin, freed):
		mask = self.stuck.copy()
		mask[self.stuck] = freed
		self.alpha[mask], self.beta[mask], self.energy[mask], self.bin[mask] =  alpha, beta, energy, bin
		# not stuck any more
		self.stuck[mask] = False
		#self.stuck[self.stuck] = freed 
	
	def pump(self):
		phi, theta, rad, alpha, beta, energy, bin = self.get_free()
		
		if self.verbose: print 'photon iteration: %d free photons, %d scattering %s' % (len(energy), self.stuck.sum(), '_'*20)
		
		# first half deals with free photons
		if len(energy) > 0:
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
			if self.verbose: print '  .. emitting %d to outside, %d inside material' % ((-inside).sum(), inside.sum())
			self.update_free(phi, theta, rad, alpha, beta, energy, bin)
			self.cut_free(inside)
			emit = dict(phi=phi0[outside], theta=theta0[outside], rad=rad0[outside], 
				beta=beta[outside], energy=energy[outside], bin=bin[outside], mask=outside)
			#print '   ', rad.shape, theta.shape, len(inside)
			xf, yf, zf = xf[inside], yf[inside], zf[inside]
			phi, theta, rad, alpha, beta, energy, bin = self.get_free()
			nphot = len(energy)
			assert nphot == inside.sum()
		
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
			
			# emit line in random direction
			if is_line.any():
				nline = photabsorbed_line.sum()
				alpha_random = numpy.random.uniform(nline)
				beta_random = acos(numpy.random.uniform(-1, 1, size=nline))
				alpha[photabsorbed_line] = alpha_random
				beta[photabsorbed_line] = beta_random
			self.update_free(phi, theta, rad, alpha, beta, energy, bin)

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
			
			#self.energy[free] = energy
			#self.bin[free] = bin
			
		  	if self.verbose: print '  .. checking if outside of energy range'
			# if in relevant energy range, find right bin
			#energy = self.energy
			# for unstuck and photabsorbed_line we have to check the energy
			# just do it for all
			energy_outside = numpy.logical_or(energy < 0.1, energy > 1100)
			dropouts = numpy.logical_or(energy_outside, absorbed)
			remainders = -dropouts
			if self.verbose: print '  .. .. outside of energy range: %d' % (energy_outside.sum())
			if self.verbose: print '  .. .. %d left' % (remainders.sum())
		
			# cut to remainders (so that we only have to compute new positions for few)
			#self.set(phi, theta, rad, alpha, beta, energy, bin)
			self.cut_free(remainders)
			phi, theta, rad, alpha, beta, energy, bin = self.get_free()
		
			# compute new positions for remainder
			xf, yf, zf = xf[remainders], yf[remainders], zf[remainders]
			rad = (xf**2+yf**2+zf**2)**0.5
			phi = atan2(yf, xf)
			theta = numpy.where(rad == 0, 0., acos(zf / rad))
			self.update_free(phi, theta, rad, alpha, beta, energy, bin)
			self.stuck[-self.stuck] = photscattered[remainders]
			if self.verbose: print '  .. .. %d stuck in scattering' % (self.stuck.sum())
		else:
			emit = None
		
		# now deal with stuck photons, i.e. those undergoing compton
		# scattering
		alpha, beta, energy, bin = self.get_stuck()
		nphotscattered = len(energy)
		if nphotscattered > 0:
			if self.verbose: print '  scattering: %d' % nphotscattered
			u = rng.uniform(size=nphotscattered)
			
			Efrac = energy / electmass
			
			mu = -1.0 + 2.0*u
			term1 = (1.0 + mu**2) / 2.0
			term2 = 1.0 / (1.0 + Efrac * (1.0 - mu))
			term3 = 2.0 * Efrac**2 * (1 - mu)**2 * term2 / term1
			
			ang_diff_scat = term1 * (term2**2) * (1.0 + term3)
			
			coin = rng.uniform(size=nphotscattered)
			processed = coin < ang_diff_scat
			if self.verbose: print '  .. passed: %d (%.1f%%)' % (processed.sum(), processed.sum()*100.0/nphotscattered)
			if processed.any():
				# we can now update the direction and energy of the 
				# stuck ones that were processed
				mu = mu[processed]
				E = energy[processed]
			
				# new energy
				#if self.verbose: print '  .. .. mus: %.2f' % (mus.mean())
				loss = (1. + (1. - mu) * E / electmass)
				if self.verbose: print '  .. .. energy loss: %.3f, %.3f, %.3f' % (
					loss.min(), loss.mean(), loss.max())
				E /= loss
				if self.verbose: print '  .. .. new energy: %.3f, %.3f, %.3f' % (
					E.min(), E.mean(), E.max())
				
				### update angle
				# beta ~ theta, alpha ~ phi in normal notation
				alpha, beta = alpha[processed], beta[processed]
				orig_vector = numpy.transpose([sin(beta)*cos(alpha), 
					sin(beta)*sin(alpha), cos(beta)])
				
				# need random vector to find rotation plane
				random_vector = numpy.random.normal(size=(mu.size, 3))
				# vector of rotation axis
				cross_vector = numpy.cross(random_vector, orig_vector)
				cross_vector /= (cross_vector**2).sum(axis=0)**0.5
				
				# we rotate by acos(mu)
				# delta_theta = acos(mu)
				delta_alpha = rng.uniform(0, 2*pi, size=mu.size)
				zz = mu.reshape((-1,1))
				sinT = (1 - mu**2)**0.5
				xx = (cos(delta_alpha) * sinT).reshape((-1,1))
				yy = (sin(delta_alpha) * sinT).reshape((-1,1))
				perturbed_vector = random_vector * xx + cross_vector * yy + orig_vector * zz
				# transform back to spherical coordinates
				xy = perturbed_vector[:,0]**2 + perturbed_vector[:,1]**2
				zz = perturbed_vector[:,2]**2
				#beta  = np.arctan2(xy**0.5, zz)
				#alpha = np.arctan2(perturbed_vector[:,1], perturbed_vector[:,0])
				_, beta_new, alpha_new = to_spherical(perturbed_vector.transpose())
				self.update_and_free_stuck(alpha_new, beta_new, E, energy2bin(E), processed)
				
	  	if self.verbose: print '  .. finally checking all, if outside of energy range'

		phi, theta, rad, alpha, beta, energy, bin = self.get()
		energy_outside = numpy.logical_or(energy < 0.1, energy > 1800)
		dropouts = energy_outside
		remainders = -dropouts
		if self.verbose: print '  .. .. outside of energy range: %d' % (energy_outside.sum())
		if self.verbose: print '  .. .. %d left' % (remainders.sum())
		# cut out
		if dropouts.any():
			self.cut(remainders)
		
		# next round
		return emit, len(self.energy) > 0

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
	
	rdata = numpy.zeros((nbins, nbins, nmu))
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
			#if n_interactions < 1:
			#	continue
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

	rdata = numpy.array(rdata)
	return rdata, nphot

def store(prefix, nphot, rdata, nmu, extra_fits_header = {}):
	
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
	for mu in range(nmu):
		a = energy_lo
		b = energy_hi
		G = PhoIndex
		weights = energy**-PhoIndex * deltae
		# .transpose()
		y = (weights * matrix[:,:,mu]).sum(axis=0) / deltae
		print '%d ... ' % mu
		
		plt.figure(figsize=(10,10))
		plt.plot(energy, 0.1 * exp(-xphot*NH) * energy**-PhoIndex, '-', color='red', linewidth=1)
		plt.plot(energy, 0.1 * exp(-xscatt*NH) * energy**-PhoIndex, '-', color='pink')
		plt.plot(energy, 0.1 * exp(-xkfe*NH) * energy**-PhoIndex, '-', color='orange')
		plt.plot(energy, 0.1 * energy**-PhoIndex, '--', color='gray')
		plt.plot(energy_lo, y / total, '-', color='k') #, drawstyle='steps')
		#plt.plot(energy, exp(-xboth) * energy**-PhoIndex, '-', color='yellow')
		plt.gca().set_xscale('log')
		plt.gca().set_yscale('log')
		#plt.xlim(0.1, 10 * (1 + 10))
		plt.xlim(3, 40)
		lo, hi = 1e-8, 1
		plt.vlines(6.40, lo, hi, linestyles=[':'], color='grey', alpha=0.5)
		plt.vlines(7.06, lo, hi, linestyles=[':'], color='grey', alpha=0.5)
		plt.ylim(lo, hi)
		plt.savefig(prefix + "_%d.pdf" % mu)
		plt.savefig(prefix + "_%d.png" % mu)
		numpy.savetxt(prefix + "_%d.txt" % mu, numpy.vstack([energy, y]).transpose())
		plt.close()
	return nphot_total, rdata


