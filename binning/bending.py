from __future__ import print_function, division
import numpy
from numpy import exp, log, log10

nbins = 1000
r = 1.5
A = log((8.1 + 0.015)/8.10)**(-1./r)

def bin2energy_lo(bin):
	bin = numpy.asarray(bin)
	with numpy.errstate(invalid='ignore'):
		return numpy.where(bin < 800, bin * 0.01 + 0.1, 8.1*exp(((bin-800)/A)**r))
def energy2bin(energy):
	x = numpy.asarray(energy)
	with numpy.errstate(invalid='ignore'):
		return (numpy.where(x < 8.1, (x-0.1) / 0.01, A*(log(x/8.1))**(1./r)+800)).astype(int)

