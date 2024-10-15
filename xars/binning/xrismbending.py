import numpy
from numpy import exp, log

nbins = 2800
r = 1.5
A = log((8.1 + 0.015) / 8.10)**(-1. / r)


def bin2energy_lo(bin):
    bin = numpy.asarray(bin)
    with numpy.errstate(invalid='ignore'):
        return numpy.where(
            bin < 600,
            bin * 0.01 + 0.1,  # up to 6.1keV
            numpy.where(
                bin < 2600,
                (bin - 600) * 0.001 + 600 * 0.01 + 0.1,  # 6.1..8.1keV
                8.1 * exp(((bin - 2600) / A)**r)
            )
        )


def energy2bin(energy):
    x = numpy.asarray(energy)
    with numpy.errstate(invalid='ignore'):
        return (numpy.where(
            x < 6.1,
            (x - 0.1) / 0.01,
            numpy.where(
                x < 8.1,
                (x - 6.1) / 0.001 + 600,
                A * (log(x / 8.1))**(1. / r) + 2600))).astype(int)
