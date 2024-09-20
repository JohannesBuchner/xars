from __future__ import print_function, division
import numpy
from numpy import log10
from . import energy2bin, bin2energy_lo, bin2energy_hi, bin2energy, nbins

import matplotlib.pyplot as plt
energy = numpy.logspace(log10(0.1), log10(1100), 10000)
results = energy2bin(energy)
bins = numpy.arange(nbins)
numpy.savetxt('binning_current.txt', numpy.transpose([bin2energy_lo(bins), bin2energy_hi(bins)]))
plt.figure(figsize=(5, 10))
plt.subplot(2, 1, 1)
#plt.plot(energy, results, '-')
#plt.plot(energy, numpy.round(400*numpy.log10(energy)), '+-')
plt.plot(bin2energy(bins)[0], bins, '+-')
plt.gca().set_xscale('log')
#plt.gca().set_yscale('log')
plt.xlabel('Energy [keV]')
plt.ylabel('Bin number')
#plt.ylim(800, 1000)
#plt.xlim(8, 12)
plt.subplot(2, 1, 2)
plt.plot(bin2energy(bins)[0], 1000 * (bin2energy_hi(bins) - bin2energy_lo(bins)))
plt.xlabel('Energy [keV]')
plt.ylabel('Energy resolution [eV]')
plt.xscale('log')
plt.yscale('log')
plt.savefig('binning.pdf')

for b in bins:
    print(b, bin2energy(b)[0], energy2bin(bin2energy(b)[0]))

