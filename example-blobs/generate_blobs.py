import h5py
import numpy

# make unit spheres at distance 3 with various densities

for nh in numpy.linspace(22, 26, 21):
	with h5py.File('torusblob%.1f.hdf5' % nh, 'w') as f:
		x = numpy.array([0.3])
		y = numpy.array([0.])
		z = numpy.array([0.])
		R = numpy.array([0.1])
		NH = numpy.array([nh])
		f['sigma'] = 0
		f['x'] = x
		f['y'] = y
		f['z'] = z
		f['radius'] = R
		f['NH'] = NH


