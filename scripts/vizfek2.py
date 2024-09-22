"""
This script demonstrates tracking photons through the simulation, and
picking up the outcoming FeK emission.
"""
import numpy

from xars.binning import bin2energy, energy2bin, nbins
from xars.coordtrans import to_cartesian
from xars.photons import PhotonBunch


def run(
    nphot, geometry,
    verbose=False, emin=0.5, PhoIndex=2
):

    energy_lo, energy_hi = bin2energy(range(nbins))
    energy = (energy_hi + energy_lo) / 2.
    # deltae = energy_hi - energy_lo

    photons = PhotonBunch(i=100, nphot=nphot, verbose=verbose, geometry=geometry)
    photons.energy = emin / numpy.random.power(PhoIndex, size=nphot)
    photons.energy[photons.energy > energy.max()] = energy.max()
    photons.binid = energy2bin(photons.energy)

    for n_interactions in range(1000):
        prev_location = photons.rad.copy(), photons.theta.copy(), photons.phi.copy()
        emission, more = photons.pump()
        if emission is None and not more:
            break
        if emission is None:
            continue
        if len(emission['energy']) == 0:
            if not more:
                break
            continue
        if verbose:
            print(' received %d emitted photons (after %d interactions)' % (len(emission['energy']), n_interactions))
        beta = emission['beta']
        alpha = emission['alpha']
        # bins = emission['binid']
        energy = emission['energy']
        mask_escaped = emission['mask_escaped']

        # get previous location
        prev_rad, prev_theta, prev_phi = prev_location
        rad, theta, phi = prev_rad[mask_escaped], prev_theta[mask_escaped], prev_phi[mask_escaped]
        assert rad.shape == energy.shape
        if n_interactions == 0:
            continue  # we do not consider direct emission

        # filter to FeKa band
        mask = numpy.logical_and(emission['energy'] > 6.1, emission['energy'] < 6.5)

        if mask.any():
            # convert last scattering position to xyz
            x, y, z = to_cartesian((rad[mask], theta[mask], phi[mask]))

            yield x, y, z, alpha[mask], beta[mask], energy[mask]

        if not more:
            break


if __name__ == '__main__':
    nh = 1e24 / 1e22

    from xars.geometries.conetorus import ConeTorusGeometry
    geometry = ConeTorusGeometry(Theta_tor=80 / 180. * numpy.pi, NH=nh, verbose=False)

    from astropy.table import Table

    events = None
    # numpy.random.seed(101)

    for _ in range(10):
        for events_info in run(1000000, geometry, emin=6, PhoIndex=2):
            new_events = numpy.transpose(events_info)
            if events is None:
                events = new_events
            else:
                events = numpy.concatenate((events, new_events))

        print("events collected:", events.shape)
        # numpy.savetxt('vizfek.txt.gz', events)
        Table(events, names=['x', 'y', 'z', 'alpha', 'beta', 'energy']).write('vizfek.fits', format='fits', overwrite=True)
