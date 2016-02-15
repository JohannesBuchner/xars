import numpy, scipy
import matplotlib.pyplot as plt
from numpy import exp, log, log10

nbins = 1250 # number of energy bins

def energy2bin(energy):
	energy = numpy.asarray(energy).reshape((-1,))
	n = 70. * numpy.log10((energy + 23.3)/23.4)
	dn = 0.01 * numpy.ones_like(n)
	ibin = numpy.round(n/dn-0.5)
	m1 = n > 8.34
	dn[m1] = 0.022
	ibin[m1] = numpy.round((n[m1]-8.34)/dn[m1]+834.-0.5)
	m2 = n > 9.308
	dn[m2] = 0.05
	ibin[m2] = numpy.round((n[m2]-9.308)/dn[m2]+878.-0.5)
	m3 = n > 10.258
	dn[m3] = 0.1
	ibin[m3] = numpy.round((n[m3]-10.258)/dn[m3]+897.-0.5)
	m4 = n > 11.158
	dn[m4] = 0.35
	ibin[m4] = numpy.round((n[m4]-11.158)/dn[m4]+906.-0.5)
	
	return numpy.array(ibin, dtype=numpy.uint)

def bin2energy(i):
	i = numpy.asarray(i).reshape((-1,))
	dn = 0.01 * numpy.ones_like(i)
	n = i * dn
	m1 = n > 8.34
	dn[m1] = 0.022
	n[m1] = 8.34 + (i[m1] - 834.0) * dn[m1]
	m2 = n > 9.308
	dn[m2] = 0.05
	n[m2] = 9.308+(i[m2]-878.0)*dn[m2]
	m3 = n > 10.258
	dn[m3] = 0.1
	n[m3] = 10.258 + (i[m3]-897.0)*dn[m3]
	m4 = n > 11.158
	dn[m4] = 0.35
	n[m4] = 11.158+(i[m4] - 906.0)*dn[m4]
	
	e=23.4*10**(n/70.0)-23.3
	f=23.4*10.**((n+dn)/70.0)-23.3
	deltae = f-e
	
	return e, f

if __name__ == '__main__':
	energy = numpy.logspace(log10(0.1), log10(1100))
	results = energy2bin(energy)
	plt.figure()
	plt.plot(energy, results, '+-')
	#plt.plot(energy, numpy.round(400*numpy.log10(energy)), '+-')
	bins = range(1250)
	plt.plot(bin2energy(bins)[0], bins, '+-')
	plt.gca().set_xscale('log')
	#plt.gca().set_xscale('log')
	plt.savefig('weirdfunc.pdf')
	
	for b in bins:
		print b, bin2energy(b)[0], energy2bin(bin2energy(b)[0])
	
