import numpy
from numpy import pi, log10
import matplotlib.pyplot as plt
from xars.geometries.clumpytorus import ClumpyTorusGeometry

plt.figure(figsize=(10, 8))
for sigma in [0, 5, 20, 60]:
	for f, fCT in [(0, 0), (8, 0.25), (6, 0.3), (4, 0.45), (3, 0.6), ]:
		Theta_tor = sigma * 1.4
		if f == 0 and sigma == 0:
			sigmas = 'dring-output/dring_empty'
		elif sigma == 0:
			sigmas = 'dring-output/dring_25.5_i0.5_n%d' % f
		else:
			sigmas = '10000_%d_f%d_gexp5core' % (sigma, f)

		geometryfile = '/mnt/data/daten/PostDoc/research/agn/torus/clumpy/%s.hdf5' % (sigmas)
		geometry = ClumpyTorusGeometry(geometryfile, verbose=True)
		#photons = PhotonBunch(i=100, nphot=10000, verbose=True, geometry=geometryfile)
		#print(photons.beta.mean(), photons.beta.min(), photons.beta.max())
		beta = pi / 2 + numpy.zeros(10000)
		alpha = numpy.linspace(0, 2 * pi, 10000 + 1)[:-1]
		nh = geometry.compute_los_nh(beta, alpha)
		lognh = log10(nh + 0.00001) + 22
		lognh[lognh<20] = 20
		plt.hist(
			lognh, cumulative=True, histtype='step', density=True,
			label='TORsigma=%d CTKcover=%.1f' % (Theta_tor, fCT),
			bins=numpy.linspace(20, 26, 40)
		)
		print('%2d %.1f %.2e %.2f %.2f %.2f' % (Theta_tor, fCT, nh.mean() * 1e22, lognh.mean(), (lognh>22).mean(), (lognh>24).mean()))
		#print(lognh.
plt.legend(loc='upper left')
plt.savefig("doc/eqcovfracs.pdf", bbox_inches='tight')
plt.savefig("doc/eqcovfracs.png", bbox_inches='tight')
plt.close()
