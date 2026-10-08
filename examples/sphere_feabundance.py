"""
Monte-Carlo simulation of a uniform sphere (central source) with a
variable iron abundance.

The sphere has edge column density NH (distance unit = NH/1e22 cm^-2), and
the cross-sections are set to a given Fe abundance (relative to solar)
before each run. The response matrix of each run is stored; afterwards the
resulting angle-averaged spectra (folded with a powerlaw) are plotted for
comparison over the full energy range and around the Fe K complex.

Usage:
  python examples/sphere_feabundance.py --nh=100 --nevents=1000
"""

import argparse
import os

import numpy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from xars import montecarlo
from xars.binning import bin2energy, nbins
from xars.geometries.spheretorus import SphereTorusGeometry
from xars.xsects import set_fe_abundance

parser = argparse.ArgumentParser(
    description="""Monte-Carlo simulator: sphere with variable Fe abundance""",
    epilog="""(C) Johannes Buchner. Based on work by Murray Brightman & Kirpal Nandra""")
parser.add_argument('--nh', type=float, default=100,
                    help='NH in 1e22 cm^-2 (100 = 1e24 cm^-2)')
parser.add_argument('--nevents', type=int, default=1000,
                    help='number of input photons per energy bin')
parser.add_argument('--fe-abundances', type=float, nargs='+', default=[0.1, 1, 10],
                    help='Fe abundances relative to solar to simulate')
parser.add_argument('--gamma', type=float, default=2,
                    help='photon index of the incident powerlaw used for folding')
parser.add_argument('--output', type=str, default='examples/output-sphere-fe',
                    help='output directory')
parser.add_argument('--seed', type=int, default=42, help='random seed')
args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)
nmu = 3  # number of viewing angle bins (sphere output is isotropic; they are summed for the plot)


def binmapfunction(beta, alpha):
    mu = ((0.5 + nmu * numpy.abs(numpy.cos(beta))) - 1).astype(int)
    mu[mu >= nmu] = nmu - 1
    return mu


prefixes = {}
for z_fe in args.fe_abundances:
    numpy.random.seed(args.seed)
    set_fe_abundance(z_fe)
    prefix = os.path.join(args.output, 'fe%g_' % z_fe)
    prefixes[z_fe] = prefix
    print('simulating Fe abundance %g (NH=%.4ge22, %d photons/bin) ...' % (
        z_fe, args.nh, args.nevents))
    geometry = SphereTorusGeometry(NH=args.nh)
    (rdata_transmit, rdata_reflect), nphot = montecarlo.run(
        prefix, nphot=args.nevents, nmu=nmu, geometry=geometry,
        binmapfunction=binmapfunction)
    rdata_transmit += rdata_reflect
    del rdata_reflect
    montecarlo.store(prefix, nphot, rdata_transmit, nmu,
                     extra_fits_header=dict(NH=args.nh, FE_ABUNDANCE=z_fe),
                     plot=False)
    del rdata_transmit

# fold the response matrices with an incident powerlaw and compare
import h5py

energy_lo, energy_hi = bin2energy(numpy.arange(nbins))
energy = (energy_lo + energy_hi) / 2.
deltae = energy_hi - energy_lo
deltae0 = deltae[energy >= 1][0]
weights = energy**-args.gamma * deltae / deltae0

spectra = {}
for z_fe, prefix in prefixes.items():
    with h5py.File(prefix + 'rdata.hdf5', 'r') as f:
        total = f.attrs['NPHOT']
        dset = f['rdata']
        y = numpy.zeros(dset.shape[1])
        # accumulate chunk-wise to keep memory usage low
        chunk = 256
        for i in range(0, dset.shape[0], chunk):
            block = numpy.asarray(dset[i:i + chunk]).sum(axis=2)
            y += (block * weights[i:i + chunk, None]).sum(axis=0)
    spectra[z_fe] = y / total / deltae * deltae0

colors = ['C0', 'C1', 'C2', 'C3', 'C4']

plt.figure(figsize=(11, 5))
plt.subplot(1, 2, 1)
for z_fe, color in zip(prefixes, colors):
    plt.plot(energy, spectra[z_fe] * energy**2, color=color,
             label='Fe %gx solar' % z_fe)
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.xlim(1, 10202)
plt.xlabel('energy [keV]')
plt.ylabel(r'$E^2\ dN/dE$ [arbitrary]')
plt.legend(loc='best')
plt.title('NH=%.4g cm^-2, full range' % (args.nh * 1e22))

plt.subplot(1, 2, 2)
for z_fe, color in zip(prefixes, colors):
    plt.plot(energy, spectra[z_fe], color=color, label='Fe %gx solar' % z_fe)
plt.xlim(5, 9.5)
plt.gca().set_yscale('log')
plt.xlabel('energy [keV]')
plt.ylabel(r'$dN/dE$ [arbitrary]')
plt.vlines([6.40, 7.06], 1e-8, 1, linestyles=':', color='grey', alpha=0.5)
plt.text(6.40, 0.7, 'Fe K$\\alpha$')
plt.text(7.06, 0.7, 'Fe K$\\beta$')
plt.legend(loc='best')
plt.title('Fe K complex (%.1f eV bins near 6.4 keV)' % 1.)
plt.ylim(1e-8, 1)
plt.tight_layout()
plt.savefig(os.path.join(args.output, 'fe_comparison.pdf'))
plt.savefig(os.path.join(args.output, 'fe_comparison.png'))
print('wrote %s' % os.path.join(args.output, 'fe_comparison.pdf'))
