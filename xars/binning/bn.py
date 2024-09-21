import numpy

"""
nbins = 1000

def energy2bin(energy):
	energy = numpy.asarray(energy).reshape((-1,))
	n = 70. * numpy.log10((energy + 23.3)/23.4)
	dn = 0.01 * numpy.ones_like(n)
	ibin = numpy.array(n/dn, dtype=int)
	m1 = n >= 8.34
	dn[m1] = 0.022
	ibin[m1] = (n[m1]- 8.340)/dn[m1]+834.
	m2 = n >= 9.308
	dn[m2] = 0.05
	ibin[m2] = (n[m2]- 9.308)/dn[m2]+878.
	m3 = n >= 10.258
	dn[m3] = 0.1
	ibin[m3] = (n[m3]-10.258)/dn[m3]+897.
	m4 = n >= 11.158
	dn[m4] = 0.22
	ibin[m4] = (n[m4]-11.158)/dn[m4]+906.
	m5 = n >= 12.918
	dn[m5] = 0.5
	ibin[m5] = (n[m5]-12.918)/dn[m5]+914.
	m6 = n >= 15.418
	dn[m6] = 1.0
	ibin[m6] = (n[m6]-15.418)/dn[m6]+919.
	
	return numpy.array(ibin, dtype=int)

def bin2energy_lo(i):
	i = numpy.asarray(i).reshape((-1,))
	dn = 0.01 * numpy.ones_like(i)
	n = i * dn
	m1 = n >= 8.34
	dn[m1] = 0.022
	n[m1] = 8.34 + (i[m1] - 834.0) * dn[m1]
	m2 = n >= 9.308
	dn[m2] = 0.05
	n[m2] = 9.308 + (i[m2] - 878.0)*dn[m2]
	m3 = n >= 10.258
	dn[m3] = 0.1
	n[m3] = 10.258 + (i[m3] - 897.0)*dn[m3]
	m4 = n >= 11.158
	dn[m4] = 0.22
	n[m4] = 11.158 + (i[m4] - 906.0)*dn[m4]
	m5 = n >= 12.918
	dn[m5] = 0.5
	n[m5] = 12.918 + (i[m5] - 914.0)*dn[m5]
	m6 = n >= 15.418
	dn[m6] = 1.0
	n[m6] = 15.418 + (i[m6] - 919.0)*dn[m6]
	
	return 23.4 * 10**(n/70.0)-23.3

"""

nbins = 1250 # number of energy bins

def energy2bin(energy):
	energy = numpy.asarray(energy).reshape((-1,))
	n = 70. * numpy.log10((energy + 23.3)/23.4)
	dn = 0.01 * numpy.ones_like(n)
	ibin = numpy.array(n/dn, dtype=int)
	m1 = n > 8.34
	dn[m1] = 0.022
	ibin[m1] = (n[m1]-8.34)/dn[m1]+834.
	m2 = n > 9.308
	dn[m2] = 0.05
	ibin[m2] = (n[m2]-9.308)/dn[m2]+878.
	m3 = n > 10.258
	dn[m3] = 0.1
	ibin[m3] = (n[m3]-10.258)/dn[m3]+897.
	m4 = n > 11.158
	dn[m4] = 0.35
	ibin[m4] = (n[m4]-11.158)/dn[m4]+906.
	
	return numpy.array(ibin, dtype=int)

def bin2energy_lo(i):
	i = numpy.asarray(i).reshape((-1,))
	dn = 0.01 * numpy.ones_like(i)
	n = i * dn
	m1 = n > 8.34
	dn[m1] = 0.022
	n[m1] = 8.34 + (i[m1] - 834.0) * dn[m1]
	m2 = n > 9.308
	dn[m2] = 0.05
	n[m2] = 9.308 + (i[m2]-878.0)*dn[m2]
	m3 = n > 10.258
	dn[m3] = 0.1
	n[m3] = 10.258 + (i[m3]-897.0)*dn[m3]
	m4 = n > 11.158
	dn[m4] = 0.35
	n[m4] = 11.158 + (i[m4] - 906.0)*dn[m4]
	
	return 23.4 * 10**(n/70.0)-23.3

