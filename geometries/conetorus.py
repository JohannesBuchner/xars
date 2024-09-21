from __future__ import print_function, division
from jax import numpy
import scipy
from jax.numpy import pi, arccos as acos, tan, round, log, log10, sin, cos, logical_and, logical_or, arctan as atan, arccos
from coordtrans import to_spherical, to_cartesian

import matplotlib.pyplot as plt

class ConeTorusGeometry(object):
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
		xv, yv, zv = to_cartesian((1, beta, alpha))
		
		#print xv, yv, zv
		# compute intersection with cone border
		a = zv**2 - (xv**2 + yv**2)*tan(0.5*pi - self.Theta_tor)**2
		b = 2.*zi*zv - (2.*xi*xv + 2.*yi*yv)*tan(0.5*pi - self.Theta_tor)**2
		c = zi**2 - (xi**2 + yi**2)*tan(0.5*pi - self.Theta_tor)**2
		quad = b**2 - 4.*a*c
		
		# compute the two solutions
		with numpy.errstate(invalid='ignore'):
			e1=(-b-quad**0.5)/(2.*a)
			e2=(-b+quad**0.5)/(2.*a)
		# if both are positive and e1<1
		#assert (quad >= 0).all()
		x1 = xi + e1 * xv
		x2 = xi + e2 * xv
		y1 = yi + e1 * yv
		y2 = yi + e2 * yv
		z1 = zi + e1 * zv
		z2 = zi + e2 * zv
		
		# check if those points are outside the sphere
		with numpy.errstate(invalid='ignore'):
			ri1 = x1**2 + y1**2 + z1**2 <= 1
			ri2 = x2**2 + y2**2 + z2**2 <= 1
		#print 'cone points radius', ri1, ri2
		
		# compute intersection with sphere
		#  a=xD2+yD2-zD2, b=2xExD+2yEyD-2zEzD, and c=xE2+yE2-zE2.
		#  a=xD2+yD2+zD2, b=2xExD+2yEyD+2zEzD, and c=xE2+yE2+zE2-1
		sa = 1
		sb = 2.*zi*zv + 2.*xi*xv + 2.*yi*yv
		sc = zi**2 + xi**2 + yi**2 - 1
		squad = sb**2 - 4.*sa*sc
		# compute the two solutions
		se1 = (-sb-squad**0.5)/(2.*sa)
		se2 = (-sb+squad**0.5)/(2.*sa)

		sx1 = xi + se1 * xv
		sx2 = xi + se2 * xv
		sy1 = yi + se1 * yv
		sy2 = yi + se2 * yv
		sz1 = zi + se1 * zv
		sz2 = zi + se2 * zv

		# check if inside cone
		with numpy.errstate(invalid='ignore'):
			r1 = sx1**2 + sy1**2 + sz1**2
			r2 = sx2**2 + sy2**2 + sz2**2
		theta1 = numpy.where(r1 == 0, 0., arccos(sz1 / r1))
		theta2 = numpy.where(r2 == 0, 0., arccos(sz2 / r2))
		theta1[theta1 > pi/2] = pi - theta1[theta1 > pi/2]
		theta2[theta2 > pi/2] = pi - theta2[theta2 > pi/2]
		si1 = theta1 >= self.Theta_tor
		si2 = theta2 >= self.Theta_tor
		#print 'circle points theta', theta1, theta2, si1, si2 #, (sx1, sy1, sz1), (sx2, sy2, sz2)
		
		# now we should go through them
		es = numpy.transpose([e1, e2, se1, se2])
		#print es, 'es'
		good = numpy.transpose([ri1, ri2, si1, si2])
		# they will be ignored, because we go in the pos direction
		with numpy.errstate(invalid='ignore'):
			good[~(es > 0)] = False
		es[~good] = -1
		nsol = good.sum(axis=1)
		#assert nsol.shape == xi.shape, (nsol.shape, xi.shape)
		es.sort(axis=1)
		# now handle the different cases
		#print nsol, good
		#assert (nsol >= 0).all(), nsol.min()
		#assert (nsol <= 3).all(), nsol.min()
		#assert ((nsol == 0) + (nsol == 1) + (nsol == 3)).all(), nsol
		
		# default:
		
		xf, yf, zf = xi + d * xv, yi + d * yv, zi + d * zv
		# no solution, we are on the brink of stepping out
		# set to False
		inside = numpy.zeros(xf.shape, dtype=bool)
		
		es_1 = nsol == 1
		t = es[es_1,-1]
		#    go from where we are to there
		inside[es_1] = d[es_1] < t # if true we never left
		#print 't:', t, d[es_1], es_1
		
		# 3 solutions
		es_2 = nsol == 3
		if es_2.any():
			#print 'BRANCH2'
			t1 = es[es_2,-3]
			t2 = es[es_2,-2]
			t3 = es[es_2,-1]
			# two cases: either we never make it to t1
			inside_2 = d[es_2] < t1 # if true we are still inside
			#print inside.shape, inside_2.shape
			inside[es_2] = inside_2
			#print 'es_2', es_2
			#print 'inside_2', inside_2
			if not inside_2.all():
				#print 'BRANCH3'
				# alternatively, we make it to t1:
				# have to add to the distance we are allowed to travel
				# the distance between t1 and t2 
				es_3 = es_2
				es_3[es_3] = ~inside_2 # select those cases
				#print 'd', d, es_3, d[es_3].shape, t2.shape, t1.shape, inside_2.shape
				#print 't:', t1, t2, t3
				#assert es_3.shape == (len(xi),)
				#assert es_3.sum() == (~inside_2).sum(), (es_3, ~inside_2)
				#assert d[es_3].size == (~inside_2).sum(), (es_3, ~inside_2)
				dmod = d[es_3] + t2[~inside_2] - t1[~inside_2]
				xf[es_3] = xi[es_3] + dmod * xv[es_3]
				yf[es_3] = yi[es_3] + dmod * yv[es_3]
				zf[es_3] = zi[es_3] + dmod * zv[es_3]
				inside[es_3] = dmod < t3[~inside_2]
		
		# bring to upper side of torus
		#zf = numpy.abs(zf)
		# compute spherical coordinates
		rad, theta, phi = to_spherical((xf, yf, zf))
		
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

def test_origin_x():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi, yi, zi = numpy.zeros((3,3))
	# x-direction
	beta  = numpy.zeros((3,)) + pi/2
	alpha  = numpy.zeros((3,))# + pi/2
	dist = numpy.array([0.01, 0.99, 1.01])
	inside, (xf,yf,zf), (rad, phi, theta) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	assert numpy.allclose(zf, 0), zf
	assert numpy.allclose(yf, 0), yf
	assert numpy.allclose(xf, dist), xf
	assert numpy.all(inside == [True, True, False]), inside
	
def test_origin_y():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi, yi, zi = numpy.zeros((3,3))
	# y-direction
	beta  = numpy.zeros((3,)) + pi/2
	alpha  = numpy.zeros((3,)) + pi/2
	dist = numpy.array([0.01, 0.99, 1.01])
	inside, (xf,yf,zf), (rad, phi, theta) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	assert numpy.allclose(zf, 0), zf
	assert numpy.allclose(yf, dist), yf
	assert numpy.allclose(xf, 0), xf
	assert numpy.all(inside == [True, True, False]), inside
	
def test_origin_z():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi, yi, zi = numpy.zeros((3,3))
	beta  = numpy.zeros((3,))
	alpha  = numpy.zeros((3,))
	dist = numpy.array([0.01, 0.99, 1.01])
	inside, (xf,yf,zf), (rad, phi, theta) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	assert numpy.allclose(zf, dist), zf
	assert numpy.allclose(yf, 0), yf
	assert numpy.allclose(xf, 0), xf
	assert numpy.all(inside == [False, False, False]), inside
	
def test_origin_xyz():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi, yi, zi = numpy.zeros((3,3))
	beta  = pi/4 + pi/8 + numpy.zeros(3)
	alpha  = pi/8 + numpy.zeros(3)
	dist = numpy.array([0.01, 0.99, 1.01])
	inside, (xf,yf,zf), (rad, phi, theta) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	assert numpy.allclose(xf, 0.85355339*dist), xf
	assert numpy.allclose(yf, 0.35355339*dist), yf
	assert numpy.allclose(zf, 0.38268343*dist), zf
	assert numpy.all(inside == [True, True, False]), inside
	assert numpy.allclose(rad, dist), rad
	assert numpy.allclose(phi, alpha), alpha
	assert numpy.allclose(theta, beta), theta

def test_vertical_x():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi, yi, zi = numpy.zeros((3,3))
	xi += 0.5
	beta  = numpy.zeros(3) # go up
	alpha = numpy.zeros(3)
	dist  = numpy.array([0.01, 0.99, 1.01])
	inside, (xf,yf,zf), (rad, phi, theta) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	assert numpy.allclose(xf, 0.5), xf
	assert numpy.allclose(yf, 0.0), yf
	assert numpy.allclose(zf, dist), zf
	assert numpy.all(inside == [True, False, False]), inside

def test_vertical_y():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi, yi, zi = numpy.zeros((3,3))
	xi += 0.5
	beta  = numpy.zeros(3) # go up
	alpha = numpy.zeros(3)
	dist  = numpy.array([0.01, 0.99, 1.01])
	inside, (xf,yf,zf), (rad, phi, theta) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	assert numpy.allclose(xf, 0.5), xf
	assert numpy.allclose(yf, 0.0), yf
	assert numpy.allclose(zf, dist), zf
	assert numpy.all(inside == [True, False, False]), inside

def test_vertical_xyz():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi = 0.4 + numpy.zeros(3)
	yi = 0.4 + numpy.zeros(3)
	zi = 0.1 + numpy.zeros(3)
	beta  = numpy.zeros(3) # go up
	alpha = numpy.zeros(3)
	dist  = numpy.array([0.01, 0.99, 1.01])
	inside, (xf,yf,zf), (rad, phi, theta) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	assert numpy.allclose(xf, 0.4), xf
	assert numpy.allclose(yf, 0.4), yf
	assert numpy.allclose(zf, 0.1 + dist), zf
	assert numpy.all(inside == [True, False, False]), inside

def test_vertical_xyz_far():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	# above that position we hit the sphere first, then the cone
	xi = 0.7 + numpy.zeros(3)
	yi = 0.7 + numpy.zeros(3)
	zi = 0.1 + numpy.zeros(3)
	beta  = numpy.zeros(3) # go up
	alpha = numpy.zeros(3)
	dist  = numpy.array([0.01, 0.8, 1.01])
	inside, (xf,yf,zf), (rad, phi, theta) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	assert numpy.allclose(xf, 0.7), xf
	assert numpy.allclose(yf, 0.7), yf
	assert numpy.allclose(zf, 0.1 + dist), zf
	assert numpy.all(inside == [True, False, False]), inside


def test_horizontal_x_shortside():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi = -0.99 + numpy.zeros(3)
	yi = 0 + numpy.zeros(3)
	zi = 0.001 + numpy.zeros(3)
	beta  = numpy.zeros(3) + pi/2 # go horizontal
	alpha = numpy.array([0, 0, pi]) # go in x
	dist  = numpy.array([0.01, 0.98, 0.02])
	inside, (xf,yf,zf), (rad, phi, theta) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	assert numpy.allclose(xf, [-0.98, -0.01, -1.01]), xf
	assert numpy.allclose(yf, 0), yf
	assert numpy.allclose(zf, 0.001), zf
	assert numpy.all(inside == [True, True, False]), inside

def test_horizontal_x_farside():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi = -0.99 + numpy.zeros(3)
	yi = 0 + numpy.zeros(3)
	zi = 0.001 + numpy.zeros(3)
	beta  = numpy.zeros(3) + pi/2 # go horizontal
	alpha = numpy.zeros(3) # go in x
	dist  = numpy.array([1.0, 1.5, 3])
	xskipped = zi * 2 # in 45 degree torus x=z
	inside, (xf,yf,zf), (rad, phi, theta) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	print('expectation:')
	print('naive:', xi + dist)
	print('skipped:', xskipped)
	print(xi + dist + xskipped)
	print('testing:')
	assert numpy.allclose(xf, xi + dist + xskipped), xf
	assert numpy.allclose(yf, 0), yf
	assert numpy.allclose(zf, 0.001), zf
	assert numpy.all(inside == [True, True, False]), inside

	


