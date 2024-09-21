import numpy

nbins = 1100

# uniform bins
def energy2bin(energy):
	b = numpy.array(((energy-2)*80. + 0.00000001), dtype=int)
	b[energy > 10] = ((10-2)*80 + (energy[energy > 10] - 10)*20 + 0.00000001)
	return b

def bin2energy_lo(binid):
	e = numpy.array(binid/80. + 2)
	e[binid > 640] = (binid[binid > 640] - 640)/20. + 10
	return e

