import numpy, scipy
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

def bin2energy_wrong(i):
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

def bin2energy_lo(i):
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
	
	return 23.4 * 10**(n/70.0)-23.3

def bin2energy_hi(i):
	return bin2energy_lo(i+1)

def bin2energy(i):
	j = numpy.asarray(i)
	return bin2energy_lo(j), bin2energy_hi(j)
	
def test_bin2energy():
	allbins = numpy.arange(nbins)
	elo, ehi = bin2energy(allbins)
	edgediff = ehi[:-1] - elo[1:]
	assert numpy.allclose(0, edgediff, atol=1e-5, rtol=1), (allbins[edgediff > 0], edgediff[edgediff > 0])

def test_energy2bin():
	E = numpy.random.uniform(0.1, 100, size=100000)
	bin = energy2bin(E)
	allbins = numpy.arange(nbins)
	elo, ehi = bin2energy(allbins)
	
	for i in allbins:
		assert (E[bin == i] >= elo[i]).all(), (E[bin == i].min(), E[bin == i].max(), elo[i], ehi[i])
		assert (E[bin == i] < ehi[i]).all(), (E[bin == i].min(), E[bin == i].max(), elo[i], ehi[i])


if __name__ == '__main__':
	import matplotlib.pyplot as plt
	energy = numpy.logspace(log10(0.1), log10(1100))
	results = energy2bin(energy)
	plt.figure()
	plt.plot(energy, results, '+-')
	#plt.plot(energy, numpy.round(400*numpy.log10(energy)), '+-')
	bins = range(nbins)
	plt.plot(bin2energy(bins)[0], bins, '+-')
	plt.gca().set_xscale('log')
	#plt.gca().set_xscale('log')
	plt.savefig('weirdfunc.pdf')
	
	for b in bins:
		print b, bin2energy(b)[0], energy2bin(bin2energy(b)[0])
	
