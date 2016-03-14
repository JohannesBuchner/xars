"""
Cross sections
---------------

Loading and computation of neutral, solar-abundance cross-sections.
"""

import numpy
import scipy
from numpy import pi, arccos as acos, tan, round, exp, log, log10, sin, cos, logical_and, logical_or
import binning

electmass = 511. # electron rest mass in keV/c^2

xscatt = numpy.zeros(binning.nbins)

energy_lo, energy_hi = binning.bin2energy(range(binning.nbins))
energy = (energy_hi + energy_lo)/2.
deltae = energy_hi - energy_lo

xthom = 6.7e-4         # Thomson cross section
x = energy / electmass # energy in units of electron rest mass

# compute scattering cross section
# Thomson regime
xscatt_thomson = xthom * (1.-(2.*x)+(26.*(x**2.)/5.))
t1=1.+2.*x
t2=1.+x
term1=t2/(x**3.)
term2=(2.*x*t2)/t1
term3=(1./(2.*x))*log(t1)
term4=(1.+3.*x)/(t1**2)
xscatt_compton = 0.75*xthom*(term1*(term2-log(t1))+term3-term4)
xscatt = numpy.where(x < 0.05, xscatt_thomson, xscatt_compton)

e1, xkfe  = numpy.loadtxt('xkfe.dat').transpose()  # cross section of K-Iron absorption?
e2, xphot = numpy.loadtxt('xphot.dat').transpose() # cross section of photon absorption?

assert (e1 == e2).all()
e70 = e1 > 70.
xphot[e70] = xkfe[e70]
#xphot = xphot * 30
assert (xphot >= 0).all()
assert (xscatt >= 0).all()

# from units 1e-21 to 1e-22
xphot *= 10
xkfe *= 10
xscatt *= 10
xscatt_thomson *= 10
xscatt_compton *= 10
xboth = xphot + xscatt
absorption_ratio = xphot / xboth
fek_ratio = 0.342 * xkfe / xphot

def test():
	assert (e1 == e2).all()
	assert e1.shape == energy.shape, (e1.shape, energy.shape)
	assert (e1 == energy_lo).all(), zip(e1, energy, energy_lo, energy_hi)

if __name__ == '__main__':
	#xscatt[:] = 1e-6
	import matplotlib.pyplot as plt
	energy = energy_lo
	plt.figure()
	plt.plot(energy, label='energy')
	plt.plot(e2, label='xphot.dat')
	plt.gca().set_yscale('log')
	plt.legend(loc='best', ncol=2, prop=dict(size=6))
	plt.savefig('binning.pdf')
	plt.close()
	
	plt.figure(figsize=(5,7))
	plt.subplot(3, 1, 1)
	plt.plot(energy, xphot, label='absorption')
	plt.plot(energy, numpy.where(xkfe < 1e-5, 1e-5, xkfe), label=r'Fe $K\alpha$')
	plt.ylim(xscatt.min() / 10, None)
	print 'absorption cross-section:', xphot
	plt.plot(energy, 1.2 * xscatt, label='scattering')
	print 'scattering cross-section:', xscatt
	plt.plot(energy, xboth, label='both')
	plt.gca().set_yscale('log')
	plt.gca().set_xscale('log')
	plt.legend(loc='best', ncol=2, prop=dict(size=6))
	plt.ylabel('cross section')
	plt.xlabel('energy [keV]')
	plt.subplot(3, 1, 2)
	#def comparison(nH, **kwargs):
	#	data = numpy.loadtxt('/tmp/%d.qdp' % nH)
	#	plt.plot(data[:,0], data[:,2], '--', **kwargs)
	plt.plot(energy, exp(-0.01*xboth), label='$N_H=10^{20}/cm^2$')
	plt.plot(energy, exp(-0.1 *xboth), label='$N_H=10^{21}/cm^2$')
	plt.plot(energy, exp(-1   *xboth), label='$N_H=10^{22}/cm^2$')
	plt.plot(energy, exp(-10  *xboth), label='$N_H=10^{23}/cm^2$')
	plt.plot(energy, exp(-100 *xboth), label='$N_H=10^{24}/cm^2$')
	#comparison(20, color='b')
	#comparison(21, color='g')
	#comparison(22, color='r')
	#comparison(23, color='c')
	#comparison(24, color='m')
	plt.plot(energy, absorption_ratio, label='vs scattering')
	plt.plot(energy, fek_ratio, label='vs scattering')
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
	plt.gca().set_yscale('log')
	plt.gca().set_xscale('log')
	plt.savefig('xsect.pdf', bbox_inches='tight')
	plt.close()


