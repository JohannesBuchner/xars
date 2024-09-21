import numpy
from numpy import pi, sin, cos, arccos, arctan2

def to_spherical(cartesian_location):
	(xf, yf, zf) = cartesian_location
	xf = numpy.asarray(xf).reshape((-1,))
	yf = numpy.asarray(yf).reshape((-1,))
	zf = numpy.asarray(zf).reshape((-1,))
	rad = (xf**2+yf**2+zf**2)**0.5
	#phi = numpy.fmod(numpy.arctan2(yf, xf) + 2*pi, pi)
	phi = arctan2(yf, xf)
	mask = ~(rad == 0)
	theta = numpy.zeros_like(rad)
	theta[mask] = arccos(zf[mask] / rad[mask])
	return (rad, theta, phi)

def to_cartesian(spherical_location):
	(rad, theta, phi) = spherical_location
	xv = rad * sin(theta) * cos(phi)
	yv = rad * sin(theta) * sin(phi)
	zv = rad * cos(theta)
	return (xv, yv, zv)


def test_random():
	for i in range(100):
		x, y, z = numpy.random.uniform(-1, 1, size=3)
		rad, theta, phi = to_spherical((x, y, z))
		assert rad <= 3**0.5, rad
		assert rad >= 0
		assert theta <= pi
		assert theta >= 0
		assert phi <= pi, phi
		assert phi >= -pi, phi
		xv, yv, zv = to_cartesian((rad, theta, phi))
	
		assert numpy.allclose(xv, x)
		assert numpy.allclose(yv, y)
		assert numpy.allclose(zv, z)
	

