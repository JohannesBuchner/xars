"""
Cross sections
---------------

Loading and computation of neutral, solar-abundance cross-sections.
"""

import matplotlib.pyplot as plt
import numpy
from numpy import exp, log

from . import (absorption_ratio, e1, energy, energy_hi, energy_lo, xboth,
               xlines, xlines_energies, xlines_relative, xphot, xscatt)


def test():
    assert e1.shape == energy.shape, (e1.shape, energy.shape)
    for i in range(len(energy)):
        assert (numpy.isclose(e1[i], energy[i])), (e1[i], energy[i], energy_lo[i], energy_hi[i])


# xscatt[:] = 1e-6
xkfe = xlines.sum(axis=1)

# energy = energy_lo
plt.figure()
plt.plot(energy, label='energy')
plt.plot(e1, label='xphot.dat')
plt.gca().set_yscale('log')
plt.legend(loc='best', ncol=2, prop=dict(size=6))
plt.savefig('binning.pdf')
plt.close()

plt.figure(figsize=(7,18))
plt.subplot(3, 1, 1)
plt.plot(energy, xphot, label='absorption')
for i in range(xlines.shape[1]):
    plt.plot(energy, xlines[:,i], label='Fluorescent Line %.2f keV' % (xlines_energies[i]))
plt.plot(energy, numpy.where(xkfe < 1e-5, 1e-5, xkfe), '--', label=r'sum')
plt.ylim(xscatt.min() / 10000, None)
# print 'absorption cross-section:', xphot
plt.plot(energy, xscatt, label='scattering')
# print 'scattering cross-section:', xscatt
plt.plot(energy, xboth, label='both')
plt.gca().set_yscale('log')
plt.gca().set_xscale('log')
plt.legend(loc='best', ncol=2, prop=dict(size=6))
plt.ylabel('cross section [${10}^{-22}$cm$^2$]')
plt.xlabel('energy [keV]')
plt.subplot(3, 1, 2)
plt.plot(energy, exp(-0.01 * xboth), label='$N_H=10^{20}/cm^2$')
plt.plot(energy, exp(-0.10 * xboth), label='$N_H=10^{21}/cm^2$')
plt.plot(energy, exp(-1.00 * xboth), label='$N_H=10^{22}/cm^2$')
plt.plot(energy, exp(-10.0 * xboth), label='$N_H=10^{23}/cm^2$')
plt.plot(energy, exp(-100. * xboth), label='$N_H=10^{24}/cm^2$')
# def comparison(nH, **kwargs):
#    data = numpy.loadtxt('/tmp/%d.qdp' % nH)
#    plt.plot(data[:,0], data[:,2], '--', **kwargs)
# comparison(20, color='b')
# comparison(21, color='g')
# comparison(22, color='r')
# comparison(23, color='c')
# comparison(24, color='m')
plt.plot(energy, absorption_ratio, label='vs scattering')
plt.plot(energy, xlines_relative[:,0], label='vs scattering')
plt.legend(loc='best', ncol=1, prop=dict(size=6))
plt.ylabel('absorption')
plt.xlabel('energy [keV]')
plt.gca().set_xscale('log')
plt.subplot(3, 1, 3)
# e^(-q*N) = prob
# N = -log(prob) / q
plt.plot(energy, -log(0.5) / xboth, label='mean free path')
plt.legend(loc='best', ncol=2, prop=dict(size=6))
plt.ylabel('distance / $N_H$')
plt.xlabel('energy [keV]')
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.savefig('xsects.pdf', bbox_inches='tight')
plt.close()
