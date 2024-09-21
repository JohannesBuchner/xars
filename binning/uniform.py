from __future__ import print_function, division
from jax import numpy
from jax.numpy import exp, log, log10

nbins = 1100
# uniform bins
def energy2bin(energy):
	b = numpy.array(((energy-2)*80. + 0.00000001), dtype=int)
	b[energy > 10] = ((10-2)*80 + (energy[energy > 10] - 10)*20 + 0.00000001)
	return b
def bin2energy_lo(bin):
	e = numpy.array(bin/80. + 2)
	e[bin > 640] = (bin[bin > 640] - 640)/20. + 10
	return e

"""
# break bins
def energy2bin(energy):
	b = numpy.array(((energy-2.158)*10 + 0.00000001), dtype=int)
	b[energy > 11.158] = ((11.158-2.158)*10 + (energy[energy > 11.158] - 11.158)/0.135 + 0.00000001)
	return b
def bin2energy_lo(bin):
	e = numpy.array(bin/10. + 2.158)
	e[bin > 90] = (bin[bin > 90] - 90)*0.135 + 11.158
	return e
# log-uniform bins
def energy2bin(energy):
	x = log10(energy)
	b = (10**(x / 10.)-1)*20.
	b = numpy.array((loge + 0.00000001), dtype=int)
	b[loge > 1] = (1 + (loge[loge > 1] - 10)*40 + 0.00000001)
	return 10**b
def bin2energy_lo(bin):
	e = numpy.log10(bin/20.+1)*10
	e[bin>400] += log10(bin/20.-400/20+1)*100
	return e
# log-linear bins
def energy2bin(energy):
	a, b, c = 10., 20., 1000.
	z = log10(400./a+1)
	x = numpy.array(energy).reshape((-1,))
	bin = a*(10**(x/b)-1)
	bin[bin>400] = a*(10**((c*z+x[bin>400])/(b+c))-1)
	return bin.astype(int)
def bin2energy_lo(bin):
	x = numpy.array(bin).reshape((-1,))
	a, b, c = 10., 20., 1000.
	return log10(x/a+1)*b + (x>400)*(log10(x/a+1)-log10(400/a+1))*c
"""

