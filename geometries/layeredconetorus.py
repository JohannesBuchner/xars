import numpy
import scipy
from numpy import pi, arccos as acos, tan, round, log, log10, sin, cos, logical_and, logical_or, arctan as atan, arccos
from coordtrans import to_spherical, to_cartesian
from raytrace import cone_raytrace_finite
import progressbar
import matplotlib.pyplot as plt

class LayeredConeTorusGeometry(object):
	def __init__(self, Theta_tors, NHs, verbose = False):
		self.Theta_tors = numpy.array(Theta_tors, dtype=float)
		self.NHs = numpy.array(NHs, dtype=float)
		self.verbose = verbose
	
	def compute_next_point(self, (xi, yi, zi), (dist, beta, alpha)):
		a, b, c = to_cartesian((1, beta, alpha))
		
		t = cone_raytrace_finite(self.Theta_tors, self.NHs, xi, yi, zi, a, b, c, dist)
		# are we still inside?
		inside = t >= 0
		
		t = numpy.where(t < 0, dist, t)
		# compute new position
		xf, yf, zf = xi + a*t, yi + b*t, zi + c*t
		
		# compute spherical coordinates
		rad, theta, phi = to_spherical((xf, yf, zf))
		
		return inside, (xf,yf,zf), (rad, theta, phi)
	
	def compute_los_nh(self, beta, alpha):
		NH_los = beta*0
		beta[beta > pi/2] = pi - beta[beta > pi/2]
		for theta, NH in zip(self.Theta_tors, self.NHs):
			NH_los[beta > theta] = NH
		return NH_los

	def viz(self): 
		""" Visualize the current geometry """
		import matplotlib
		import matplotlib.pyplot as plt
		from matplotlib.collections import PatchCollection
		import matplotlib.path as mpath
		import matplotlib.patches as mpatches
		import matplotlib.lines as mlines
		plt.figure(figsize=(5,5))
		font = 'sans-serif'
		ax = plt.axes([0,0,1,1])
		
		describe_all = len(self.Theta_tors) <= 3
		ax.add_line(mlines.Line2D([0,0.9], [0.5,0.5], lw=1.,alpha=0.4, ls='dashed', color='grey'))
		ax.add_line(mlines.Line2D([0.4,0.4], [0.5,0.9], lw=1.,alpha=0.4, ls='dashed', color='grey'))
		
		for i, (Theta_tor, NH) in enumerate(zip(self.Theta_tors, self.NHs)):
			Theta_tor = Theta_tor * 180 / pi
			nh = log10(NH) + 22
			thickness = max(0, min(1, (nh - 20.) / 5))
			plt.text(0.35, 0.5,  "nH=%2.1f" % nh, ha="right", va='center',
				family=font, size=14)
			#print thickness
			if describe_all or i == 0:
				ax.add_patch(  mpatches.Arc((0.4,0.5), 0.5, 0.5, theta2=90, theta1=90 - Theta_tor, 
					color='black', linewidth=1, alpha=1, fill=False, ls='dashed') )
			ax.add_patch( mpatches.Wedge((0.4,0.5), 0.3, -90 + Theta_tor, 90 - Theta_tor, color=(1-thickness, 1-thickness, 1-thickness), 
				linewidth=0, fill=True) )
			ax.add_patch( mpatches.Wedge((0.4,0.5), 0.3, 90 + Theta_tor, -90 - Theta_tor, color=(1-thickness, 1-thickness, 1-thickness), 
				linewidth=0, fill=True) )
		
		if describe_all:
			text = ', '.join(["%2.0f" % (Theta_tor * 180 / pi) for Theta_tor in self.Theta_tors])
		else:
			text = ', '.join(["%2.0f" % (Theta_tor * 180 / pi) for Theta_tor in self.Theta_tors]) + '...'
		plt.text(0.4 + 0.02, 0.5 + 0.25 + 0.02, text, ha="left", va='bottom',
			family=font, size=14)
		ax.add_patch( mpatches.Circle((0.4,0.5), 0.02, color='red', 
			linewidth=0, alpha=1, fill=True) )
		
		ax.set_xticks([])
		ax.set_yticks([])

class ConeTorusGeometry(LayeredConeTorusGeometry):
	def __init__(self, Theta_tor, NH, verbose = False):
		self.Theta_tor = Theta_tor
		self.NH = NH
		self.verbose = verbose
		super(ConeTorusGeometry, self).__init__(Theta_tors=[Theta_tor], NHs=[NH], verbose=verbose)
	
def test_origin_x():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi, yi, zi = numpy.zeros((3,3))
	# x-direction
	beta  = numpy.zeros((3,)) + pi/2
	alpha  = numpy.zeros((3,))# + pi/2
	dist = numpy.array([0.01, 0.99, 1.01])
	inside, (xf,yf,zf), (rad, theta, phi) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	assert numpy.allclose(zf, 0), zf
	assert numpy.allclose(yf, 0), yf
	assert numpy.allclose(xf, dist), (xf, dist)
	assert numpy.all(inside == [True, True, False]), inside
	
def test_origin_y():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi, yi, zi = numpy.zeros((3,3))
	# y-direction
	beta  = numpy.zeros((3,)) + pi/2
	alpha  = numpy.zeros((3,)) + pi/2
	dist = numpy.array([0.01, 0.99, 1.01])
	inside, (xf,yf,zf), (rad, theta, phi) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
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
	inside, (xf,yf,zf), (rad, theta, phi) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
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
	inside, (xf,yf,zf), (rad, theta, phi) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	assert numpy.allclose(xf, 0.85355339*dist), xf
	assert numpy.allclose(yf, 0.35355339*dist), yf
	assert numpy.allclose(zf, 0.38268343*dist), zf
	assert numpy.all(inside == [True, True, False]), inside
	print rad, theta, phi
	assert numpy.allclose(rad, dist), rad
	assert numpy.allclose(theta, beta), (beta, theta)
	assert numpy.allclose(phi, alpha), (alpha, phi)

def test_vertical_x():
	torus = ConeTorusGeometry(Theta_tor=45 * pi / 180, NH=1, verbose=True)
	xi, yi, zi = numpy.zeros((3,3))
	xi += 0.5
	beta  = numpy.zeros(3) # go up
	alpha = numpy.zeros(3)
	dist  = numpy.array([0.01, 0.99, 1.01])
	inside, (xf,yf,zf), (rad, theta, phi) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
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
	inside, (xf,yf,zf), (rad, theta, phi) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
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
	inside, (xf,yf,zf), (rad, theta, phi) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
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
	inside, (xf,yf,zf), (rad, theta, phi) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
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
	inside, (xf,yf,zf), (rad, theta, phi) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
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
	inside, (xf,yf,zf), (rad, theta, phi) = torus.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	print 'expectation:'
	print 'naive:', xi + dist
	print 'skipped:', xskipped
	print xi + dist + xskipped
	print 'testing:'
	assert numpy.all(inside == [True, True, False]), inside
	assert numpy.allclose(yf, 0), yf
	assert numpy.allclose(zf, 0.001), zf
	assert numpy.allclose(xf[inside], (xi + dist + xskipped)[inside]), xf[inside]

if __name__ == '__main__':
	import matplotlib.pyplot as plt
	torus = LayeredConeTorusGeometry(Theta_tors=numpy.array([30,45,60,80]) * pi / 180, NHs=[0.1,1,10,100], verbose=True)
	torus.viz()
	plt.savefig('layeredconetorus.pdf', bbox_inches='tight')
	plt.close()
	
	


