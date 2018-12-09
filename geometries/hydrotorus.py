from __future__ import print_function, division
import numpy
import scipy
from numpy import pi, arccos as acos, tan, round, log, log10, sin, cos, logical_and, logical_or, arctan as atan
from coordtrans import to_spherical, to_cartesian
import progressbar
import h5py
import matplotlib.pyplot as plt
from raytrace import grid_raytrace_finite_flat, grid_raytrace_flat

class HydroTorusGeometry(object):
	"""
	A uniform grid of densities is provided as a input file
	"""
	def __init__(self, filename, verbose = False):
		f = h5py.File(filename, 'r')
		rho = numpy.array(f['rho'].value, dtype=numpy.float64)
		pc_cm = 3.0856776e+18
		nH_Msun = 1.18803e+57
	
		# Compute total mass by multiplying with
		# convert from pc to grid cells
		# convert from pc to cm
		# density in Msun/pc^3
		# 32 pc side length of grid
		# convert density to n(H)/cm^3
		# so multiply by m(H)/Msun, and by (pc/cm)^3
		# in units of 1e22
		self.rho = rho * nH_Msun / pc_cm**3 * (32 * pc_cm / 256) / 1e22
		self.rho_flat = numpy.array(self.rho.flatten())
		self.rho = self.rho_flat.reshape(self.rho.shape)
		assert (self.rho >= 0).all()
		#print self.rho.sum(axis=1).max()
		self.center = f['center'].value
		self.verbose = verbose
	
	def compute_next_point(self, location, direction):
		(xi, yi, zi) = location
		(dist, beta, alpha) = direction
		a, b, c = to_cartesian((1, beta, alpha))
		x = xi + self.center[0]
		y = yi + self.center[1]
		z = zi + self.center[2]
		t = grid_raytrace_finite_flat(self.rho_flat, self.rho.shape[0], x, y, z, a, b, c, dist)
		
		# are we still inside?
		inside = t >= 0
		
		# compute new position
		t = numpy.where(t < 0, dist, t)
		xf, yf, zf = xi + a*t, yi + b*t, zi + c*t
		
		# compute spherical coordinates
		rad, phi, theta = to_spherical((xf, yf, zf))
		
		return inside, (xf,yf,zf), (rad, phi, theta)
	
	def compute_los_nh(self, beta, alpha):
		a, b, c = to_cartesian((1, beta, alpha))
		x = a*0 + self.center[0]
		y = b*0 + self.center[1]
		z = c*0 + self.center[2]
		NH = grid_raytrace_flat(self.rho_flat, self.rho.shape[0], x, y, z, a, b, c)
		return NH
	
	def viz(self): 
		""" Visualize the current geometry """
		import matplotlib.pyplot as plt
		
		fig = plt.figure(figsize=(5,5), frameon=False)
		plt.axis('off')
		ax = plt.gca()
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)
		cmap = plt.cm.gray_r
		cmap = 'Greens'
		logrho = log10(self.rho[:,127,:] + 1e-3)
		plt.imshow(logrho.transpose(), cmap=cmap, vmin=-3, vmax=+3)
		plt.plot(self.center[2], self.center[0], 'x', color='r', ms=4, mew=0.5)
		plt.xlim(0, len(logrho))
		plt.ylim(0, len(logrho))
		"""
		for i in numpy.linspace(1, 179, 100):
			zero = numpy.zeros(100) + 0.
			dist = 100 + zero
			beta = i / 180. * pi + zero
			alpha = numpy.linspace(0, pi, 100)
			inside1, _, _ = self.compute_next_point((zero, zero, zero), (dist, beta, alpha))
			inside2, _, _ = self.compute_next_point((zero, zero, zero), (dist, pi-beta, alpha))
			NH = self.compute_los_nh(beta, alpha)
			print '%.1f %.1f%% %.1f%% NH=%.2f' % (i, inside1.mean()*100, inside2.mean()*100, NH.mean())
		"""
		


