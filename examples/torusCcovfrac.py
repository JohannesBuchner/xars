import h5py
import numpy
from numpy import cos, log10
import matplotlib.pyplot as plt
import sys

plt.figure(figsize=(10, 8))
for i, sigma in enumerate([0, 5, 20, 60][::-1]):
	for j, (f, fCT) in enumerate([(0, 0), (8, 0.25), (6, 0.3), (4, 0.45), (3, 0.6), ]):
		Theta_tor = sigma * 1.4
		if f == 0 and sigma == 0:
			sigmas = 'dring-output/dring_empty'
		elif sigma == 0:
			sigmas = 'dring-output/dring_25.5_i0.5_n%d' % f
		else:
			sigmas = '10000_%d_f%d_gexp5core' % (sigma, f)

		geometryfile = '/mnt/data/daten/PostDoc/research/agn/torus/clumpy/%s.hdf5' % (sigmas)
		f = h5py.File(geometryfile, 'r')
		lognh = numpy.log10(f['NH_samples'][:])
		lognh[lognh<20] = 20
		#lognh = 22 + log10(nh)
		plt.hist(
			lognh, cumulative=True, histtype='step', density=True,
			label='TORsigma=%d CTKcover=%.1f' % (Theta_tor, fCT),
			bins=numpy.linspace(20, 26, 40)
		)
		print('%d %.1f %.2f %.2f' % (Theta_tor, fCT, (lognh>22).mean(), (lognh>24).mean()))
		#print(lognh.
plt.legend(loc='lower right')
plt.savefig("doc/covfracs.pdf", bbox_inches='tight')
plt.savefig("doc/covfracs.png", bbox_inches='tight')
plt.close()
