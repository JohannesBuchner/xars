import numpy
from numpy import pi, tan, round, log, log10, sin, cos, logical_and, logical_or, arccos, arctan, arctan2

def to_spherical((xf, yf, zf)):
	xf = numpy.asarray(xf).reshape((-1,))
	yf = numpy.asarray(yf).reshape((-1,))
	zf = numpy.asarray(zf).reshape((-1,))
	rad = (xf**2+yf**2+zf**2)**0.5
	#phi = numpy.fmod(numpy.arctan2(yf, xf) + 2*pi, pi)
	phi = numpy.arctan2(yf, xf)
	theta = numpy.where(rad == 0, 0., arccos(zf / rad))
	return (rad, theta, phi)

def to_cartesian((rad, theta, phi)):
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
	

