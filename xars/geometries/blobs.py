from __future__ import print_function, division
import numpy
from numpy import sin, cos

class MultiBlobGeometry(object):
	def __init__(self, blobs, verbose = False):
		self.blobs = blobs
		self.verbose = verbose
	
	def compute_next_point(self, location, direction):
		# unit direction vector
		(xi, yi, zi) = location
		(dist, beta, alpha) = direction
		xv=sin(beta)*cos(alpha)
		yv=sin(beta)*sin(alpha)
		zv=cos(beta)
		
		nh_contributions = [] # total nh in this direction
		# go through each blob
		for blob in self.blobs:
			nh = blob['NH']
			xb,yb,zb = blob['pos']
			rb = blob['r']
			
			#   compute intersections
			l = numpy.transpose([xv, yv, zv])
			l2 = numpy.square(l).sum(axis=1)
			assert len(a) == len(xv)
			oc = numpy.transpose([xi - xb, yi - yb, zi - zb])
			loc = (l * oc).sum(axis=1)
			ococrb = (oc * oc).sum(axis=1) - rb**2
			
			det = loc**2 - l2 * (ococrb)
			nosol  = det < 0
			twosol = det >= 0
			
			d1 = -loc + numpy.sqrt(det)
			d2 = -loc - numpy.sqrt(det)
			
			# points are at i + v * d12
			
			# if d1 or d2 are negative, we should ignore it (would be going backwards)
			# then the "intersection" point is the starting point
			d1[(d1 < 0) | nosol] = 0
			d2[(d2 < 0) | nosol] = 0
			
			dist_within = numpy.abs(d1 - d2)
			
			nh_contributions.append(dist_within * nh / (2 * rb))
			
		#   compute distance inside the blob
		#   compute distance outside the blob
		# compute average NH
		# compute distance traveled in units of average NH
		# 
		
		d = dist / self.NH # distance in units of nH
		
		if self.verbose: print('  .. .. mean in nH units: ', d.mean())
		
		# compute relative vector traveled
		xv=d*sin(beta)*cos(alpha)
		yv=d*sin(beta)*sin(alpha)
		zv=d*cos(beta)
		
		# compute new position
		xf=xi+xv
		yf=yi+yv
		zf=zi+zv
		
		if self.verbose: print('  .. computing crossing with cone')
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
		if self.verbose: print('  .. %d of %d have 2 solutions' % (twosolmask.sum(), len(twosolmask)))
		# compute the two possible new positions
		x1=xi[twosolmask]+e1[twosolmask]*xv[twosolmask]
		x2=xi[twosolmask]+e2[twosolmask]*xv[twosolmask]
		y1=yi[twosolmask]+e1[twosolmask]*yv[twosolmask]
		y2=yi[twosolmask]+e2[twosolmask]*yv[twosolmask]
		z1=zi[twosolmask]+e1[twosolmask]*zv[twosolmask]
		z2=zi[twosolmask]+e2[twosolmask]*zv[twosolmask]
		
		ltsol = e2[twosolmask] < 1.
		
		xf[twosolmask][ltsol] = x2[ltsol]+xv[ltsol]-(x1[ltsol]-xi[ltsol])
		yf[twosolmask][ltsol] = y2[ltsol]+yv[ltsol]-(y1[ltsol]-yi[ltsol])
		zf[twosolmask][ltsol] = z2[ltsol]+zv[ltsol]-(z1[ltsol]-zi[ltsol])
		
		xf[twosolmask][~ltsol] += (x2[~ltsol]-x1[~ltsol])
		yf[twosolmask][~ltsol] += (y2[~ltsol]-y1[~ltsol])
		zf[twosolmask][~ltsol] += (z2[~ltsol]-z1[~ltsol])
		
		#print '  .. using symmetries'
		# use symmetries
		# bring to upper side of torus
		zf = numpy.abs(zf)
		# compute spherical coordinates
		rad = (xf**2+yf**2+zf**2)**0.5
		phi = numpy.where(xf == 0., 0, atan(yf / xf))
		theta = acos(zf / rad)
		
		if self.verbose: print('  .. checking if left cone')
		# are we inside the cone?
		inside = numpy.logical_and(rad < 1., theta > self.Theta_tor)
		return inside, (xf,yf,zf), (rad, phi, theta)
	
	def viz(self): 
		""" Visualize the current geometry """
		Theta_tor = self.Theta_tor * 180 / pi
		import matplotlib
		import matplotlib.pyplot as plt
		from matplotlib.collections import PatchCollection
		import matplotlib.path as mpath
		import matplotlib.patches as mpatches
		import matplotlib.lines as mlines
		plt.figure(figsize=(5,5))
		font = 'sans-serif'
		ax = plt.axes([0,0,1,1])
	
		thickness = max(0, min(1, (self.NH - 20.) / 5))
		plt.text(0.35, 0.5,  "nH=%2.1f" % self.NH, ha="right", va='center',
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

