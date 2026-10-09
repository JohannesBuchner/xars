"""
Cross sections
---------------

Loading and computation of neutral, solar-abundance cross-sections.
"""

import os

import numpy
import scipy
from numpy import arccos as acos
from numpy import (cos, exp, log, log10, logical_and, logical_or, pi, round,
                   sin, tan)

from xars import binning

electmass = 511.  # electron rest mass in keV/c^2

xscatt = numpy.zeros(binning.nbins)

energy_lo, energy_hi = binning.bin2energy(numpy.arange(binning.nbins))
energy = (energy_hi + energy_lo) / 2.
deltae = energy_hi - energy_lo

xthom = 6.7e-4          # Thomson cross section
x = energy / electmass  # energy in units of electron rest mass

# compute scattering cross section
# Thomson regime
xscatt_thomson = xthom * (1. - (2. * x) + (26. * (x**2.) / 5.))
t1 = 1. + 2. * x
t2 = 1. + x
term1 = t2 / (x**3.)
term2 = (2. * x * t2) / t1
term3 = (1. / (2. * x)) * log(t1)
term4 = (1. + 3. * x) / (t1**2)
xscatt_compton = 0.75 * xthom * (term1 * (term2 - log(t1)) + term3 - term4)
xscatt = numpy.where(x < 0.05, xscatt_thomson, xscatt_compton)

# convert to units of 1e-22 cm^2 (from 1e-21)
xscatt *= 10
xscatt_thomson *= 10
xscatt_compton *= 10

# When applied the cross section is 120% larger
xscatt_thomson *= 1.2
xscatt_compton *= 1.2
xscatt *= 1.2

# photoelectric and line cross-sections
# find xsects.dat file next to this file (or as given by the environment,
# e.g. a file written by scripts/write_xsects_feabundance.py)
xsects_filename = os.environ.get(
    'XARS_XSECTS_FILE',
    os.path.join(os.path.dirname(__file__), 'xsects.dat'))
xsectsdata = numpy.loadtxt(xsects_filename)
xlines_energies = xsectsdata[0,2:]
xlines_yields = xsectsdata[1,2:]
xlines_fwhm = xsectsdata[2,2:]
xlines_asymmetries = xsectsdata[3,2:]
xsects = xsectsdata[4:,:]
# convert to units of 1e-22 cm^2 (from 1e-21)
xsects[:,1:] *= 10
e1 = xsects[:,0]
xphot = xsects[:,1]
e70 = e1 > 70.
if e70.any():
    lines_max = numpy.max(xsects[:,2:], axis=1)
    xphot[e70] = lines_max[e70] * xphot[e70][0] / lines_max[e70][0]
assert (xphot >= 0).all()
assert (xscatt >= 0).all()

xlines = xsects[:,2:] * xlines_yields
xlines_relative = xlines / xphot[:,None]
# compute probability to come out as fluorescent line
xlines_cumulative = numpy.cumsum(xlines_relative, axis=1)
assert (xlines >= 0).all()
assert (xlines_relative >= 0).all()
assert (xlines_cumulative >= 0).all()

xboth = xphot + xscatt
absorption_ratio = xphot / xboth

# columns of the cross-section table belonging to iron: parse the column
# names from the table header comment of xsects.dat
# (e.g. "# E XPHOT XKFEa2 XKFEa1 XKFEb XKC XKO ... XKNIa2 XKNIa1")
with open(xsects_filename) as f:
    colnames = [line for line in f if line.startswith('# E XPHOT')][-1].lstrip('#').split()
fe_columns = numpy.array(
    [2 + i for i, name in enumerate(colnames[2:]) if name.startswith('XKFE')],
    dtype=int)
assert len(fe_columns) > 0, 'no iron (XKFE*) columns found in xsects.dat header'
# xlines drops the E and XPHOT columns
fe_xlines_columns = fe_columns - 2

# pristine solar-abundance copies; set_fe_abundance restores these so that
# repeated calls set the abundance absolutely, not relative to the last call
# (xphot is a view into xsects, so restoring xsects restores it as well)
_solar_xsects = xsects.copy()
_solar_xlines = xlines.copy()
_solar_xphot = xphot.copy()


def set_fe_abundance(z_fe):
    """Set the iron abundance relative to the solar value (1 = solar).

    The Fe columns of the cross-section table give the Fe K-shell
    photoionisation cross-section, which drives the fluorescent line
    production; the same scaling is applied to the photoelectric absorption
    above the Fe K edge, where Fe dominates the absorption. Call before
    running the Monte-Carlo simulation; can be called repeatedly.
    """
    z_fe = float(z_fe)
    xsects[:] = _solar_xsects
    xlines[:] = _solar_xlines
    fe_photo = xsects[:, fe_columns[0]].copy()
    xsects[:, fe_columns] *= z_fe
    xlines[:, fe_xlines_columns] *= z_fe
    xphot[:] = _solar_xphot + (z_fe - 1) * fe_photo
    assert (xphot >= 0).all(), 'Fe photoelectric scaling made xphot negative'
    xlines_relative[:] = xlines / xphot[:, None]
    xlines_cumulative[:] = numpy.cumsum(xlines_relative, axis=1)
    xboth[:] = xphot + xscatt
    absorption_ratio[:] = xphot / xboth


def test():
    assert e1.shape == energy.shape, (e1.shape, energy.shape)
    for i in range(len(energy)):
        assert (numpy.isclose(e1[i], energy[i])), (e1[i], energy[i], energy_lo[i], energy_hi[i])
