from __future__ import print_function, division
from jax import numpy
import scipy
from jax.numpy import pi, arccos as acos, tan, round, log, log10, sin, cos, logical_and, logical_or, arctan as atan, arccos
from coordtrans import to_spherical, to_cartesian

import matplotlib.pyplot as plt

class BNTorusGeometry(object):
	def __init__(self, Theta_tor, NH, verbose = False):
		self.Theta_tor = Theta_tor
		self.NH = NH
		self.verbose = verbose
	
	def compute_next_point(self, location, direction):
		(xi, yi, zi) = location
		(dist, beta, alpha) = direction
		d = dist / self.NH # distance in units of nH
		
		if self.verbose: print('  .. .. mean in nH units: ', d.mean())

		# compute relative vector traveled
		xv, yv, zv = to_cartesian((d, beta, alpha))
		
		# compute new position
		xf, yf, zf = xi + xv, yi + yv, zi + zv
		
		# compute intersection with cone border
		a = zv**2 - (xv**2 + yv**2)*tan(0.5*pi - self.Theta_tor)**2
		b = 2.*zi*zv - (2.*xi*xv + 2.*yi*yv)*tan(0.5*pi - self.Theta_tor)**2
		c = zi**2 - (xi**2 + yi**2)*tan(0.5*pi - self.Theta_tor)**2
		quad = b**2 - 4.*a*c
		
		# compute the two solutions
		e1=(-b-quad**0.5)/(2.*a)
		e2=(-b+quad**0.5)/(2.*a)
		# if both are positive and e1<1
		twosolmask = logical_and(logical_and(quad > 0., e1 > 0.),
			logical_and(e1 < 1., e2 > 0.))
	  	#if self.verbose: print '  .. %d of %d have 2 solutions' % (twosolmask.sum(), len(twosolmask))
		# compute the two possible new positions
		#print 'twosol:', twosolmask
		
		x1=xi[twosolmask]+e1[twosolmask]*xv[twosolmask]
		x2=xi[twosolmask]+e2[twosolmask]*xv[twosolmask]
		y1=yi[twosolmask]+e1[twosolmask]*yv[twosolmask]
		y2=yi[twosolmask]+e2[twosolmask]*yv[twosolmask]
		z1=zi[twosolmask]+e1[twosolmask]*zv[twosolmask]
		z2=zi[twosolmask]+e2[twosolmask]*zv[twosolmask]
		
		#print 'e2', e2
		ltsol = e2[twosolmask] < 1.
		gtsol = ~ltsol
		
		asol = twosolmask.copy()
		asol[twosolmask] = ltsol
		bsol = twosolmask.copy()
		bsol[twosolmask] = gtsol
		xf[asol] = x2[ltsol]+xv[asol]-(x1[ltsol]-xi[asol])
		yf[asol] = y2[ltsol]+yv[asol]-(y1[ltsol]-yi[asol])
		zf[asol] = z2[ltsol]+zv[asol]-(z1[ltsol]-zi[asol])
		
		xf[bsol] += (x2[gtsol]-x1[gtsol])
		yf[bsol] += (y2[gtsol]-y1[gtsol])
		zf[bsol] += (z2[gtsol]-z1[gtsol])
		
	  	#print '  .. using symmetries'
		# use symmetries
		# bring to upper side of torus
		zf = numpy.abs(zf)
		# compute spherical coordinates
		rad, theta, phi = to_spherical((xf, yf, zf))
		assert not numpy.isnan(rad).any()
		assert not numpy.isnan(theta).any()
		assert not numpy.isnan(phi).any()
	  	#if self.verbose: print '  .. checking if left cone'
		# are we inside the cone?
		inside = numpy.logical_and(rad < 1., theta > self.Theta_tor)
		return inside, (xf,yf,zf), (rad, phi, theta)

	def viz(self): 
		""" Visualize the current geometry """
		Theta_tor = self.Theta_tor * 180 / pi
		nh = log10(self.NH) + 22
		
		import matplotlib
		import matplotlib.pyplot as plt
		from matplotlib.collections import PatchCollection
		import matplotlib.path as mpath
		import matplotlib.patches as mpatches
		import matplotlib.lines as mlines
		plt.figure(figsize=(5,5))
		font = 'sans-serif'
		ax = plt.axes([0,0,1,1])
		
		thickness = max(0, min(1, (nh - 20.) / 5))
		plt.text(0.35, 0.5,  "nH=%2.1f" % nh, ha="right", va='center',
			family=font, size=14)
		#print thickness
		ax.add_line(mlines.Line2D([0,0.9], [0.5,0.5], lw=1.,alpha=0.4, ls='dashed', color='grey'))
		ax.add_line(mlines.Line2D([0.4,0.4], [0.5,0.9], lw=1.,alpha=0.4, ls='dashed', color='grey'))
		ax.add_patch(  mpatches.Arc((0.4,0.5), 0.5, 0.5, theta2=90, theta1=90 - Theta_tor, 
			color='black', linewidth=1, alpha=1, fill=False, ls='dashed') )
		plt.text(0.4 + 0.02, 0.5 + 0.25 + 0.02, "%2.0f" % Theta_tor, ha="left", va='bottom',
			family=font, size=14)
		ax.add_patch( mpatches.Wedge((0.4,0.5), 0.3, -90 + Theta_tor, 90 - Theta_tor, color='black', 
			linewidth=0, alpha=thickness, fill=True) )
		ax.add_patch( mpatches.Wedge((0.4,0.5), 0.3, 90 + Theta_tor, -90 - Theta_tor, color='black', 
			linewidth=0, alpha=thickness, fill=True) )
	
		ax.add_patch( mpatches.Circle((0.4,0.5), 0.02, color='red', 
			linewidth=0, alpha=1, fill=True) )
		
		ax.set_xticks([])
		ax.set_yticks([])




