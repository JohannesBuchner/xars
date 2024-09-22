import numpy
import tqdm
from numpy import cos, exp, pi, sin

from .binning import bin2energy, energy2bin, nbins
from .photons import PhotonBunch
from .xsects import (absorption_ratio, electmass, xboth, xlines_cumulative,
                     xlines_energies, xphot, xscatt)

rng = numpy.random


def plot_interaction(nphot, n_interactions, rad, theta, beta, **kwargs):
    import matplotlib.pyplot as plt
    mask = rad < 1
    xi = rad * sin(theta)
    zi = rad * cos(theta)
    plt.figure(figsize=(5,5))
    plt.ylim(0,1)
    plt.xlim(0,1)
    plt.plot(
        xi[mask], zi[mask], '.', ms=1, alpha=0.5,
        label='%02.2f%% (%02.2f%% out) at interaction # %d' % (
            len(mask) * 100. / nphot,
            mask.sum() * 100. / nphot,
            n_interactions), color='red')
    plt.plot(
        xi[~mask], zi[~mask], '.', ms=1, alpha=0.5,
        label='%02.2f%% (%02.2f%% through) at interaction # %d' % (
            len(mask) * 100. / nphot,
            (~mask).sum() * 100. / nphot,
            n_interactions), color='blue')
    plt.vlines(xi[mask].mean(), 0, 1, linestyle='--', color='red')
    plt.vlines(xi[~mask].mean(), 0, 1, linestyle='--', color='blue')
    xv = sin(beta)
    zv = cos(beta)
    plt.plot(xv[mask], zv[mask], '+', color='red', ms=1)
    plt.plot(xv[~mask], zv[~mask], '+', color='blue', ms=1)

    plt.plot(sin(beta[mask].mean()), cos(beta[mask].mean()), 'o', color='red')
    plt.plot(sin(beta[~mask].mean()), cos(beta[~mask].mean()), 'o', color='blue')
    plt.legend(loc='best', prop=dict(size=6))


def plot_path(rad_paths, theta_paths, **kwargs):
    import matplotlib.pyplot as plt
    for rad, theta in zip(rad_paths, theta_paths):
        xi = rad * sin(theta)
        zi = rad * cos(theta)
        print(theta, rad, xi, zi)
        plt.plot(xi, zi, 'o:', **kwargs)
    print()


def run(
    prefix, nphot, nmu, geometry,
    binmapfunction,
    plot_paths=False, plot_interactions=False, verbose=False
):

    if plot_paths or plot_interactions:
        import matplotlib.pyplot as plt

    rdata_transmit = numpy.zeros((nbins, nbins, nmu))
    rdata_reflect = numpy.zeros((nbins, nbins, nmu))
    energy_lo, energy_hi = bin2energy(list(range(nbins)))
    energy = (energy_hi + energy_lo) / 2.
    # deltae = energy_hi - energy_lo

    binrange = [list(range(nbins + 1)), list(range(nmu + 1))]
    for i in tqdm.trange(nbins - 1, -1, -1):
        photons = PhotonBunch(i=i, nphot=nphot, verbose=verbose, geometry=geometry)
        remainder = [(photons.rad, photons.theta)]
        if plot_paths:
            plt.figure("paths", figsize=(4, 4))
        for n_interactions in range(1000):
            emission, more = photons.pump()
            if emission is None and not more:
                break
            if emission is None:
                continue
            if len(emission['energy']) == 0:
                if not more:
                    break
                if plot_paths:
                    mask_kept = emission['mask_kept']
                    remainder = [
                        (prev_rad[mask_kept], prev_theta[mask_kept])
                        for prev_rad, prev_theta in remainder] + [(photons.rad, photons.theta)]
                continue
            if verbose:
                print(' received %d emitted photons (after %d interactions)' % (len(emission['energy']), n_interactions))
            beta = emission['beta']
            alpha = emission['alpha']
            assert (beta <= pi).all(), beta
            assert (beta >= 0).all(), beta
            # vertical bin, i.e. which viewing angle can see this photon?
            mbin = numpy.asarray(binmapfunction(beta=beta, alpha=alpha)).astype(numpy.uint)
            # highest bin exceeded due to rounding
            mbin[mbin == nmu] = nmu - 1

            bins = emission['binid']
            # produce unique array bins, mbin which contains counts
            counts, xedges, yedges = numpy.histogram2d(bins, mbin, bins=binrange)
            # record into histogram if it landed within relevant range
            if n_interactions < 1:
                rdata_transmit[i] += counts
            else:
                rdata_reflect[i] += counts
            # if (emission['energy'] == 6.40).any():
            #    print ' %f: %d/%d are lines' % (energy[i], (emission['energy'] == 7.06).sum(), (emission['energy'] == 6.40).sum())
            #    linebin = set(bins[emission['energy'] == 6.40])
            #    print linebin, rdata[i][724,:], rdata[900][724,4] if i > 900 else ''

            # remove the emitted photons from the remainder
            if plot_paths:
                mask_escaped = emission['mask_escaped']
                mask_kept = emission['mask_kept']
                path = [(prev_rad[mask_escaped], prev_theta[mask_escaped]) for prev_rad, prev_theta in remainder]

                remainder = [
                    (prev_rad[mask_kept], prev_theta[mask_kept])
                    for prev_rad, prev_theta in remainder] + [(photons.rad, photons.theta)]

                rad_paths = numpy.transpose([path_rad for path_rad, path_theta in path] + [2 * numpy.ones(mask_escaped.sum())])
                theta_paths = numpy.transpose([path_theta for path_rad, path_theta in path] + [beta])
                plt.figure("paths")
                plot_path(
                    rad_paths[:100], theta_paths[:100],
                    color=['b', 'r', 'g', 'y', 'k', 'm', 'w'][n_interactions % 7],
                    alpha=1 - 0.75 * numpy.exp(-n_interactions / 5.))

            if plot_interactions:
                print('plotting %d photons ...' % len(beta))
                plot_interaction(nphot=nphot, n_interactions=n_interactions, **emission)
                plt.savefig(prefix + "rdata_%d_%d.png" % (i, n_interactions))
                plt.close()
                print('plotting ... done')
            if not more:
                break
        if plot_paths:
            plt.figure("paths")
            plt.xlim(-1, 1)
            plt.ylim(-1, 1)
            plt.title('Energy: %.2f, %d interactions' % (energy[i], n_interactions))
            plt.savefig(prefix + "paths_%d.png" % (i))
            plt.close()

    return (rdata_transmit, rdata_reflect), nphot


def store(prefix, nphot, rdata, nmu, extra_fits_header=None, plot=False):
    energy_lo, energy_hi = bin2energy(list(range(nbins)))
    energy = (energy_hi + energy_lo) / 2.
    deltae = energy_hi - energy_lo
    import h5py
    try:
        print('Loading previous file "%s" if exists ...' % (prefix + "rdata.hdf5"))
        with h5py.File(prefix + "rdata.hdf5", 'r') as old_file:
            print('  reading values')
            rdata_old = old_file['rdata'][()]
            print('  reading header value')
            prev_nphot = old_file.attrs['NPHOT']
            print('Accumulating onto previous result ...')
            rdata = rdata + rdata_old
            del rdata_old
    except Exception as e:
        if "error message = 'no such file or directory'" not in str(e).lower():
            print('updating file failed; writing fresh. Error:', e)
        prev_nphot = 0
    nphot_total = nphot + prev_nphot
    print('  storing ...')
    import datetime
    import time
    now = datetime.datetime.fromtimestamp(time.time())
    nowstr = now.isoformat()
    nowstr = nowstr[:nowstr.rfind('.')]
    with h5py.File(prefix + "rdata.hdf5", 'w') as f:
        f.create_dataset('rdata', data=rdata, compression='gzip', shuffle=True)
        f.create_dataset('energy_lo', data=energy_lo, compression='gzip', shuffle=True)
        f.create_dataset('energy_hi', data=energy_hi, compression='gzip', shuffle=True)

        f.attrs['CREATOR'] = """Johannes Buchner <johannes.buchner.acad@gmx.com>"""
        f.attrs['DATE'] = nowstr
        f.attrs['METHOD'] = 'Monte-Carlo simulation code'
        f.attrs['NPHOT'] = nphot_total
        if extra_fits_header is not None:
            for k,v in extra_fits_header.items():
                f.attrs[k] = v

    print('  total of %d input / %d output photons across %d bins' % (nphot_total, rdata.sum(), nbins))

    if not plot:
        return nphot_total, rdata
    import sys

    import matplotlib.pyplot as plt
    PhoIndex = 2
    matrix = rdata
    total = nphot_total
    deltae0 = deltae[energy >= 1][0]
    NH = 1e24 / 1e22
    weights = (energy**-PhoIndex * deltae / deltae0).reshape((-1,1))  # * deltae.reshape((1, -1)) / deltae.reshape((-1, 1))
    yall = (weights * matrix.sum(axis=2)).sum(axis=0) / deltae * deltae0
    for mu in range(nmu):
        y = (weights * matrix[:,:,mu]).sum(axis=0) / deltae * deltae0
        sys.stdout.write('plotting %d/%d ...\r' % (mu + 1, nmu))
        sys.stdout.flush()

        plt.figure(figsize=(10,10))
        plt.plot(energy, exp(-xphot * NH) * energy**-PhoIndex, '-', color='red', linewidth=1)
        plt.plot(energy, exp(-xscatt * NH) * energy**-PhoIndex, '-', color='pink')
        plt.plot(energy, exp(-xlines_cumulative[:,0] * NH) * energy**-PhoIndex, '-', color='orange')
        plt.plot(energy, energy**-PhoIndex, '--', color='gray')
        plt.plot(energy, y / total * nmu, '-', color='k')  # , drawstyle='steps')
        plt.plot(energy, yall / total, '-', color='gray', alpha=0.3, linewidth=3)  # , drawstyle='steps')
        # plt.plot(energy, exp(-xboth) * energy**-PhoIndex, '-', color='yellow')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        # plt.xlim(0.1, 10 * (1 + 10))
        plt.xlim(1, 40)
        lo, hi = 1e-8, 1
        plt.vlines(6.40, lo, hi, linestyles=[':'], color='grey', alpha=0.5)
        plt.vlines(7.06, lo, hi, linestyles=[':'], color='grey', alpha=0.5)
        # plt.vlines(8.33, lo, hi, linestyles=[':'], color='grey', alpha=0.5)
        # plt.vlines(9.4, lo, hi, linestyles=[':'], color='grey', alpha=0.5)
        # plt.vlines(10.4, lo, hi, linestyles=[':'], color='grey', alpha=0.5)
        plt.ylim(lo, hi)
        # plt.ylim(1e-5, 1e-3)
        # plt.xlim(6, 15)
        # plt.show()
        # plt.savefig(prefix + "_%d.pdf" % mu)
        plt.savefig(prefix + "_%d.png" % mu)
        # numpy.savetxt(prefix + "_%d.txt" % mu, numpy.vstack([energy, y]).transpose())
        plt.close()
    print()
    return nphot_total, rdata
