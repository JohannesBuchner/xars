import numpy
import scipy
from numpy import exp, log, log10

# pick a binning module
# from bn import nbins, bin2energy_lo, energy2bin
# from uniform import nbins, bin2energy_lo, energy2bin
from .bending import bin2energy_lo, energy2bin, nbins


def bin2energy_hi(i):
    return bin2energy_lo(i + 1)


def bin2energy(i):
    j = numpy.asarray(i)
    return bin2energy_lo(j), bin2energy_hi(j)


def test_bin2energy():
    allbins = numpy.arange(nbins)
    elo, ehi = bin2energy(allbins)
    edgediff = ehi[:-1] - elo[1:]
    assert numpy.allclose(0, edgediff, atol=1e-5, rtol=1), (allbins[edgediff > 0], edgediff[edgediff > 0])


def test_energy2bin():
    E = numpy.random.uniform(0.5, 100, size=100000)
    binid = energy2bin(E)
    allbins = numpy.arange(nbins)
    elo, ehi = bin2energy(allbins)

    for i in allbins:
        assert (E[binid == i] >= elo[i]).all(), (i, E[binid == i].min(), E[binid == i].max(), elo[i], ehi[i])
        assert (E[binid == i] < ehi[i]).all(), (i, E[binid == i].min(), E[binid == i].max(), elo[i], ehi[i])


def test_reversible_bins():
    E = numpy.random.uniform(9.5, 12, size=100000)
    binid = energy2bin(E)
    elo, ehi = bin2energy(binid)
    assert (elo <= E).all(), (elo[~(elo <= E)], E[~(elo <= E)])
    assert (ehi > E).all(), (ehi[~(ehi > E)], E[~(ehi > E)])
    erand = numpy.random.uniform(elo, ehi)
    bin2 = energy2bin(erand)
    assert (bin2 == binid).all()


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
    plt.plot(bin2energy(bins)[0], bins, '+-')
    plt.gca().set_xscale('log')
    plt.xlabel('Energy [keV]')
    plt.ylabel('Bin number')
    plt.savefig('binning.pdf')

    for b in bins:
        print(b, bin2energy(b)[0], energy2bin(bin2energy(b)[0]))
