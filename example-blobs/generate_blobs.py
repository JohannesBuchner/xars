from __future__ import print_function, division
import h5py
import numpy
from numpy import log10, log, pi, exp, arctan2, arcsin, cos
import matplotlib.pyplot as plt

# make unit spheres at distance 3 with various densities

for nh in numpy.arange(21, 26.5, 0.5):
	plt.figure(figsize=(5,5))
	for fcov in numpy.logspace(-5, 0, 21):
		R = 0.1
		phi = 2 * pi * fcov
		# compute distance
		# assuming angle subtended is :
		# 2pi * (1 - sqrt(d^2 - R^2) / d)
		# solve(phi=2pi * (1 - sqrt(d^2 - R^2) / d), d)
		# https://www.wolframalpha.com/input/?i=solve(phi%3D2pi+*+(1+-+sqrt(d%5E2+-+R%5E2)+%2F+d),+d)
		d = 2 * pi * R / (4 * pi * phi - phi**2)**0.5
		
		if d < 2 and int(log10(fcov)) == log10(fcov):
			circle = plt.Circle((d, 0), R)
			circle.set_facecolor('none')
			circle.set_edgecolor('k')
			plt.gca().add_artist(circle)
			plt.text(d, 0.12, '$f_\mathrm{cov}=%d$' % log10(fcov), 
				va='bottom', ha='center', rotation=90)
		
		phimax = arcsin(R / d)
		area = cos(pi/2 - phimax) * (2 * phimax)
		area_all = 1 * 2 * pi
		if nh == 22:
			print('fcov=', fcov, 'd=', d, 'angular diameter=', 2*arctan2(R, 2*d)/pi*180)
			print('   irradiated area: %.5f/%.5f = %.5f%% (phi=%.5f)' % (area, area_all, area / area_all * 100, phimax))
		
		with h5py.File('blob_nh%.1f_fcov%.2f.hdf5' % (nh, numpy.log10(fcov)), 'w') as f:
			x = numpy.array([d])
			y = numpy.array([0.])
			z = numpy.array([0.])
			R = numpy.array([R])
			NH = numpy.array([nh])
			f['sigma'] = 0
			f['x'] = x
			f['y'] = y
			f['z'] = z
			f['radius'] = R
			f['NH'] = NH
			f['fcov'] = fcov
			f['phi'] = phimax
			f['radfrac'] = area / area_all
	plt.plot(0, 0, 'x', color='r', ms=4, mew=2)
	plt.xlim(-1, 1)
	plt.ylim(-1, 1)
	plt.savefig('geometry.pdf', bbox_inches='tight')
	plt.savefig('geometry.png', bbox_inches='tight')
	plt.close()
