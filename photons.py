from __future__ import print_function, division
import numpy
import scipy
from numpy import pi, arccos as acos, tan, round, log, log10, sin, cos, logical_and, logical_or, arctan as atan, arctan2 as atan2, exp
from coordtrans import to_spherical, to_cartesian
from binning import nbins, energy2bin, bin2energy
from xsects import xscatt, xphot, xlines_cumulative, xboth, absorption_ratio, xlines_energies, electmass

rng = scipy.random

def projected_distance(x, y, z, alpha, beta):
	"""
	Project point (rad, theta, phi) onto line from (0, 0, 0) in direction (alpha, beta)
	and return the distance (which may be negative).
	"""
	vx, vy, vz = to_cartesian((1, alpha, beta))
	#x, y, z = to_cartesian((rad, theta, phi))
	return vx * x + vy * y + vz * z


class PhotonBunch(object):
	def __init__(self, i, # input energy bin
		nphot, # number of photons
		geometry,
		verbose = False):
		self.verbose = verbose
		self.rad = numpy.zeros(nphot)
		self.phi = numpy.zeros(nphot) # initial left-right angle for position
		self.theta = numpy.zeros(nphot) # initial up-down angle for position
		self.alpha = rng.uniform(0, 2*pi, size=nphot) # initial left-right direction for direction -- value does not matter due to symmetry
		self.distance = numpy.zeros(nphot) # total distance travelled
		#mu = numpy.linspace(-1, 1, nphot) # uniform distribution
		mu = rng.uniform(-1, 1, size=nphot)
		self.beta = acos(mu) # up-down angle for direction
		# self.beta = acos(numpy.linspace(-cone_in, cone_in, nphot)) # up-down angle for direction
		self.geometry = geometry
		energy_lo, energy_hi = bin2energy(i)
		e = (energy_lo + energy_hi) / 2.
		if self.verbose: print('PhotonBunch of size %d with energy %.2f keV' % (nphot, e))
		#self.energy = e * numpy.ones(nphot)
		self.energy = rng.uniform(low=energy_lo, high=energy_hi, size=nphot)
		#bin = energy2bin(self.energy)
		#assert (bin == i).all(), (bin.min(), bin.max(), self.energy.min(), self.energy.max(), energy_lo, energy_hi, self.energy[bin!=i])
		self.bin = i * numpy.ones(nphot, dtype=numpy.uint)
		self.stuck = self.rad != 0 # False
	
	def cut(self, mask):
		# TODO: we could also refill with new photons
		self.phi, self.theta, self.rad = self.phi[mask], self.theta[mask], self.rad[mask]
		self.alpha, self.beta = self.alpha[mask], self.beta[mask]
		self.energy, self.bin = self.energy[mask], self.bin[mask]
		self.distance = self.distance[mask]
		self.stuck = self.stuck[mask]
	
	def cut_free(self, free_mask):
		"""cut the free photons according to free_mask """
		free = ~self.stuck
		free[free] = free_mask
		self.cut(free)
		return free
	
	def set(self, phi, theta, rad, alpha, beta, energy, bin):
		self.phi, self.theta, self.rad, self.alpha, self.beta, self.energy, self.bin = phi, theta, rad, alpha, beta, energy, bin
	
	def get(self):
		return self.phi, self.theta, self.rad, self.alpha, self.beta, self.energy, self.bin
	
	def get_free(self):
		free = ~self.stuck
		return self.phi[free], self.theta[free], self.rad[free], self.alpha[free], self.beta[free], self.energy[free], self.bin[free]

	def update_free(self, phi, theta, rad, alpha, beta, energy, bin):
		free = ~self.stuck
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
		if self.verbose: print('photon iteration: %d free photons, %d scattering %s' % ((~self.stuck).sum(), self.stuck.sum(), '_'*20))
		phi, theta, rad, alpha, beta, energy, bin = self.get_free()
		
		
		# first half deals with free photons
		if len(energy) > 0:
			if self.verbose: print('  .. distance travelled')
			nphot = len(energy)
			r1 = rng.uniform(size=nphot)
			taur = - log(1. - r1) # effective tau
			dist = taur / xboth[bin] # distance travelled
			#if self.verbose: print '  .. .. tau min: %f mean: %f max: %f ' % (taur.min(), taur.mean(), taur.max())
			if self.verbose: print('  .. .. dist min: %f mean: %f max: %f ' % (dist.min(), dist.mean(), dist.max()))
		
			phi0 = phi
			theta0 = theta
			rad0 = rad
			xi, yi, zi = to_cartesian((rad0, theta0, phi0))
			
			#if self.verbose: print '  .. computing position'
			# compute position
			inside, (xf,yf,zf), (rad, phi, theta) = self.geometry.compute_next_point((xi, yi, zi), (dist, beta, alpha))
			distance = ((xi - xf)**2 + (yi - yf)**2 + (zi - zf)**2)**0.5
			outside = ~inside

			# emit
			if self.verbose: print('  .. emitting %d to outside, %d inside material' % ((~inside).sum(), inside.sum()))
			self.update_free(phi, theta, rad, alpha, beta, energy, bin)
			mask = ~self.stuck
			mask[mask] = inside
			self.distance[mask] += distance[inside]
			if self.verbose: print('  .. distance travelled: %f +- %f' % (self.distance.mean(), self.distance.std()))
			assert (self.distance >= 0).all(), self.distance[~(self.distance >= 0)]
			mask = ~self.stuck
			mask[mask] = outside
			# normalise distance measure relative to photons
			# directly from origin in direction of alpha, beta
			self.distance[mask] -= projected_distance(
				xi[outside], yi[outside], zi[outside], 
				alpha[outside], beta[outside])
			emit = dict(phi=phi0[outside], theta=theta0[outside], rad=rad0[outside], 
				beta=beta[outside], alpha=alpha[outside],
				energy=energy[outside], bin=bin[outside], 
				x=xi[outside], y=yi[outside], z=zi[outside],
				mask=mask, distance=self.distance[mask])
			self.cut_free(inside)
			#print '   ', rad.shape, theta.shape, len(inside)
			xf, yf, zf = xf[inside], yf[inside], zf[inside]
			phi, theta, rad, alpha, beta, energy, bin = self.get_free()
			nphot = len(energy)
			assert nphot == inside.sum()
		
			if self.verbose: print('  .. checking if absorbed')
			# what effect are we dealing with
			r2 = rng.uniform(size=nphot)
		
			# Photon-absorbed
			q = absorption_ratio[bin]
			photabsorbed = r2 < q
			#if self.verbose: print '  .. .. photon-abs prob min: %f mean: %f max: %f ' % (q.min(), q.mean(), q.max())
		
			# Iron to photon-absorption effectiveness
			omega = xlines_cumulative[bin[photabsorbed]] # Importance of lines compared to photon-absorption
			# omega *= 0 # disable fluorescent lines
			#print '  ..  .. omega:', omega
			r3 = rng.uniform(size=photabsorbed.sum())
			
			# Are we coming out as a line
			# This handles all lines (lines loaded in xsects)
			line_mask = omega > r3[:,None]
			iline = numpy.where(line_mask.any(axis=1),
				line_mask.argmax(axis=1), -1)
			is_line = iline >= 0
			photabsorbed_line = photabsorbed.copy()
			photabsorbed_notline = photabsorbed.copy()
			# set the absorbed ones (where we have true) to the criterion
			photabsorbed_line[photabsorbed] = is_line
			photabsorbed_notline[photabsorbed] = ~is_line
		
			r4 = rng.uniform(size=photabsorbed_line.sum())
			
			if is_line.any():
				e2 = xlines_energies[iline[is_line]]
				e = energy2bin(e2)
				energy[photabsorbed_line] = e2
				bin[photabsorbed_line] = e
				del e2, e
			
			# emit line in random direction
			if photabsorbed_line.any():
				nline = photabsorbed_line.sum()
				alpha_random = rng.uniform(0, 2*pi, size=nline)
				beta_random = acos(rng.uniform(-1, 1, size=nline))
				alpha[photabsorbed_line] = alpha_random
				beta[photabsorbed_line] = beta_random
			self.update_free(phi, theta, rad, alpha, beta, energy, bin)

			# no, we are just being swallowed (photon-absorption).
			# photabsorbed & ~is_line
			# remove photabsorbed_line 
		
			absorbed = photabsorbed_notline
			if self.verbose: print('  .. .. line emission: %3d' % (is_line.sum()))
			if self.verbose: print('  .. .. absorbed: %d (%d lost)' % (photabsorbed.sum(), absorbed.sum()))
		
			if self.verbose: print('  .. checking if scattered')
			# Were we compton-scattered?
			photscattered = (~photabsorbed)
			nphotscattered = photscattered.sum()
			#assert nphotscattered.shape == energy.shape
			if nphotscattered > 0:
				if self.verbose: print('  .. .. scattered: %d' % nphotscattered)
			
			#self.energy[free] = energy
			#self.bin[free] = bin
			
			if self.verbose: print('  .. checking if outside of energy range')
			# if in relevant energy range, find right bin
			#energy = self.energy
			# for unstuck and photabsorbed_line we have to check the energy
			# just do it for all
			energy_outside = bin <= 0
			dropouts = numpy.logical_or(energy_outside, absorbed)
			remainders = ~dropouts
			if self.verbose: print('  .. .. outside of energy range: %d' % (energy_outside.sum()))
			if self.verbose: print('  .. .. %d left' % (remainders.sum()))
		
			# cut to remainders (so that we only have to compute new positions for few)
			#self.set(phi, theta, rad, alpha, beta, energy, bin)
			self.cut_free(remainders)
			phi, theta, rad, alpha, beta, energy, bin = self.get_free()
		
			# compute new positions for remainder
			xf, yf, zf = xf[remainders], yf[remainders], zf[remainders]
			rad, theta, phi = to_spherical((xf, yf, zf))
			#rad = (xf**2+yf**2+zf**2)**0.5
			#phi = atan2(yf, xf)
			#theta = numpy.where(rad == 0, 0., acos(zf / rad))
			self.update_free(phi, theta, rad, alpha, beta, energy, bin)
			self.stuck[~self.stuck] = photscattered[remainders]
			if self.verbose: print('  .. .. %d stuck in scattering' % (self.stuck.sum()))
		else:
			emit = None
		
		# now deal with stuck photons, i.e. those undergoing compton
		# scattering
		alpha, beta, energy, bin = self.get_stuck()
		nphotscattered = len(energy)
		if nphotscattered > 0:
			if self.verbose: print('  scattering: %d' % nphotscattered)
			a = rng.uniform(size=nphotscattered)
			# compute new direction:
			alpha_new = a * 2. * pi # left-right angle uniform randomly

			# compute up-down angle
			beta0 = beta
			r5 = rng.uniform(size=nphotscattered)
			r5a = rng.uniform(size=nphotscattered)

			x = 2. * r5a - 1
			mus = numpy.where(r5 > 0.75, 
				numpy.sign(x) * numpy.abs(x)**(1/3.), 
				x)
			betas = acos(mus)
			mu = cos(beta0)*cos(betas)+sin(beta0)*sin(betas)*cos(alpha_new - alpha)
			beta_new = acos(mu)

			# new energy
			#if self.verbose: print '  .. .. mus: %.2f' % (mus.mean())
			loss = (1. + (1. - mus) * energy / electmass)
			if self.verbose: print('  .. .. energy loss: %.3f, %.3f, %.3f' % (
				loss.min(), loss.mean(), loss.max()))
			E = energy / loss
			if self.verbose: print('  .. .. new energy: %.3f, %.3f, %.3f' % (
				E.min(), E.mean(), E.max()))
			processed = E >= 0
			self.update_and_free_stuck(alpha_new, beta_new, E, energy2bin(E), processed)
		
		if self.verbose: print('  .. finally checking all, if outside of energy range')

		phi, theta, rad, alpha, beta, energy, bin = self.get()
		energy_outside = numpy.logical_or(energy < 0.1, energy > 1800)
		dropouts = energy_outside
		remainders = ~dropouts
		if self.verbose: print('  .. .. outside of energy range: %d' % (energy_outside.sum()))
		if self.verbose: print('  .. .. %d left' % (remainders.sum()))
		# cut out
		if dropouts.any():
			self.cut(remainders)
		
		# next round
		return emit, len(self.energy) > 0

