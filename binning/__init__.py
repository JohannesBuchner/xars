import numpy, scipy
from numpy import exp, log, log10

#from bn import nbins, bin2energy_lo, energy2bin
#from uniform import nbins, bin2energy_lo, energy2bin
from bending import nbins, bin2energy_lo, energy2bin

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
	#E = numpy.random.uniform(6.275, 6.2875, size=1000000)
	E = numpy.random.uniform(0.5, 100, size=100000)
	bin = energy2bin(E)
	allbins = numpy.arange(nbins)
	elo, ehi = bin2energy(allbins)
	
	for i in allbins:
		assert (E[bin == i] >= elo[i]).all(), (i, E[bin == i].min(), E[bin == i].max(), elo[i], ehi[i])
		assert (E[bin == i] <  ehi[i]).all(), (i, E[bin == i].min(), E[bin == i].max(), elo[i], ehi[i])

def test_reversible_bins():
	E = numpy.random.uniform(9.5, 12, size=100000)
	bin = energy2bin(E)
	elo, ehi = bin2energy(bin)
	assert (elo <= E).all(), (elo[~(elo <= E)], E[~(elo <= E)])
	assert (ehi > E).all(), (ehi[~(ehi > E)], E[~(ehi > E)])
	erand = numpy.random.uniform(elo, ehi)
	bin2 = energy2bin(erand)
	assert (bin2 == bin).all()

def test_reversible():
	allbins = numpy.arange(nbins)
	elo, ehi = bin2energy(allbins)
	emid = (elo + ehi) / 2.
	imid = energy2bin(emid)
	assert (imid == allbins).all(), (allbins[(imid != allbins)], imid[(imid != allbins)])
	ilo = energy2bin(elo + 0.0001)
	assert (ilo == allbins).all(), (allbins[(ilo != allbins)], ilo[(ilo != allbins)])
	ihi = energy2bin(ehi - 0.0001)
	assert (ihi == allbins).all(), (allbins[(ihi != allbins)], ihi[(ihi != allbins)])

if __name__ == '__main__':
	import matplotlib.pyplot as plt
	energy = numpy.logspace(log10(0.1), log10(1100), 10000)
	results = energy2bin(energy)
	bins = numpy.arange(nbins)
	numpy.savetxt('binning_current.txt', numpy.transpose([bin2energy_lo(bins), bin2energy_hi(bins)]))
	plt.figure()
	#plt.plot(energy, results, '-')
	#plt.plot(energy, numpy.round(400*numpy.log10(energy)), '+-')
	plt.plot(bin2energy(bins)[0], bins, '+-')
	plt.gca().set_xscale('log')
	#plt.gca().set_yscale('log')
	plt.xlabel('Energy [keV]')
	plt.ylabel('Bin number')
	#plt.ylim(800, 1000)
	#plt.xlim(8, 12)
	plt.savefig('binning.pdf')
	
	for b in bins:
		print b, bin2energy(b)[0], energy2bin(bin2energy(b)[0])
	
