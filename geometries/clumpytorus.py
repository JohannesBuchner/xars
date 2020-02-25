from __future__ import print_function, division
import numpy
import scipy
from numpy import pi, arccos as acos, tan, round, log, log10, sin, cos, logical_and, logical_or, arctan as atan
from coordtrans import to_spherical, to_cartesian
import h5py
import matplotlib.pyplot as plt
from raytrace import sphere_raytrace_finite, sphere_raytrace

def sigma_convert(sigma):
	return round((sigma * 1.4) / 5) * 5

class ClumpyTorusGeometry(object):
	def __init__(self, filename, verbose = False):
		f = h5py.File(filename, 'r')
		self.sigma = f.get('sigma')
		# load spheres
		self.x, self.y, self.z = f['x'][()], f['y'][()], f['z'][()]
		self.r = f['radius'][()]
		self.NH = 10**(f['NH'][()] - 22)
		# density = NH of crossing / distance of crossing
		self.rho = self.NH / (2 * self.r)
		self.verbose = verbose
	
	def compute_next_point(self, location, direction):
		(xi, yi, zi) = location
		(dist, beta, alpha) = direction
		a, b, c = to_cartesian((1, beta, alpha))
		t = sphere_raytrace_finite(self.x, self.y, self.z, self.r, self.rho, xi, yi, zi, a, b, c, dist)
		
		# are we still inside?
		inside = t >= 0
		
		# compute new position
		xf, yf, zf = xi + a*t, yi + b*t, zi + c*t
		
		# compute spherical coordinates
		rad, phi, theta = to_spherical((xf, yf, zf))
		
		return inside, (xf,yf,zf), (rad, phi, theta)
	
	def compute_los_nh(self, beta, alpha):
		a, b, c = to_cartesian((1, beta, alpha))
		mindistances = numpy.zeros(1)
		NH = sphere_raytrace(self.x, self.y, self.z, self.r, self.rho, a, b, c, mindistances)[0,:]
		return NH
	
	def viz(self): 
		""" Visualize the current geometry """
		
		import matplotlib
		import matplotlib.pyplot as plt
		from matplotlib.collections import PatchCollection
		import matplotlib.patches as mpatches
		import matplotlib.lines as mlines
		
		x, y, z = self.x, self.y, self.z
		NH = log10(self.NH) + 22
		r = self.r
		intersecting = numpy.abs(y) < r

		plt.figure(figsize=(5,5), frameon=False)
		plt.axis('off')
		if self.sigma is not None:
			plt.title(r'$\sigma=%d^\circ$' % (sigma_convert(self.sigma[()])))
		cmap = plt.cm.gray_r
		patches = []
		colors = []
		for x, y, z, r, NH in zip(x[intersecting], y[intersecting], z[intersecting], r[intersecting], NH[intersecting]):
			r2 = r - numpy.abs(y)
			circle = plt.Circle((x,z),r2)
			patches.append(circle)
			#colors.append((min(26, max(20, NH))-20)/(26-20.))
			colors.append((min(26, max(20, NH))))
		
		collection = PatchCollection(patches, cmap=plt.cm.gray_r, edgecolors="none")
		collection.set_array(numpy.array(colors))
		collection.set_clim(20, 26)
		coll = plt.gca().add_collection(collection)
		
		plt.plot(0, 0, 'x ', color='r', ms=4, mew=2)
		plt.ylim(-1, 1)
		plt.xlim(-1, 1)
		
		if self.sigma is None or self.sigma[()] < 30:
			# add colorbar
			ax = plt.axes([0.1, 0.1, 0.8, 0.02], frameon=False)
			cbar = plt.colorbar(coll, cax=ax, ticks=[20, 21, 22, 23, 24, 25, 26],
				cmap=plt.cm.gray_r, orientation='horizontal')
			cbar.solids.set_edgecolor("face")
			cbar.outline.set_linewidth(0)
			cbar.set_label('Cloud column density')
			
