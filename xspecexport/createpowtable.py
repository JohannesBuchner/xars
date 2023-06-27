import numpy
from numpy import exp, log
import astropy.io.fits as pyfits

outfilename = 'myatable.fits'

# create some energy grid --- START
# you can replace this with your own grid

nbins = 1000
def bin2energy_lo(i):
	"""Get lower energy edge of bin with index i"""
	# get lower edge of bin, in keV
	r = 1.5
	A = log((8.1 + 0.015)/8.10)**(-1./r)
	with numpy.errstate(invalid='ignore'):
		return numpy.where(i < 800, i * 0.01 + 0.1, 8.1*exp(((i-800)/A)**r))

energy_lo = bin2energy_lo(numpy.arange(nbins))
energy_hi = bin2energy_lo(numpy.arange(nbins) + 1)
# create some energy grid --- END

# central energy of bin
energy = (energy_hi + energy_lo) / 2
# width of bin
deltae = energy_hi - energy_lo

table = []
PhoIndices = [ 1.        ,  1.20000005,  1.39999998,  1.60000002,  1.79999995,
		2.        ,  2.20000005,  2.4000001 ,  2.5999999 ,  2.79999995,
		3. ]
Ecuts = [20.,  30, 40, 60, 100, 140, 200, 400]

# width of bin at 1keV
deltae0 = deltae[energy >= 1][0]
print(deltae, deltae0)

for PhoIndex in PhoIndices:
	for Ecut in Ecuts:
		# compute model: here a powerlaw with an exponential cut-off. 
		# We need to integrate over the energy bin, which we do here by 
		# evaluating at the centre and multiplying by the width
		mymodel = energy**-PhoIndex * exp(-energy / Ecut) * deltae
		# add parameters and spectrum:
		table.append(((PhoIndex, Ecut), mymodel))

hdus = []
hdu = pyfits.PrimaryHDU()
import datetime, time
now = datetime.datetime.fromtimestamp(time.time())
nowstr = now.isoformat()
hdu.header['CREATOR'] = """Johannes Buchner <johannes.buchner.acad@gmx.com>"""
hdu.header['DATE'] = nowstr[:nowstr.rfind('.')]
hdu.header['HDUCLAS1'] = 'XSPEC TABLE MODEL'
hdu.header['MODLNAME'] = 'torus'
hdu.header['ADDMODEL'] = True
hdu.header['MODLUNIT'] = 'photons/cm^2/s'
hdu.header['EXTEND']   = True
hdu.header['REDSHIFT'] = True
hdu.header['SIMPLE']   = True
hdu.header['HDUDOC']   = 'OGIP/92-009'
hdu.header['HDUVERS1'] = '1.0.0'
hdu.header['HDUCLASS'] = 'OGIP'
hdus.append(hdu)

# parameters of this model (must agree with the loops above!)

# NAME, METHOD, INITIAL, DELTA, MINIMUM, BOTTOM, TOP, MAXIMUM, NUMBVALS, VALUE (41)
dtype = [('NAME', 'S12'), ('METHOD', '>i4'), ('INITIAL', '>f4'), ('DELTA', '>f4'), ('MINIMUM', '>f4'), ('BOTTOM', '>f4'), ('TOP', '>f4'), ('MAXIMUM', '>f4'), ('NUMBVALS', '>i4'), ('VALUE', '>f4', (41,))]

parameters = numpy.array([
	('PhoIndex', 0, 2.0, 0.0099999998, 1.0, 1.2, 2.8, 3.0, 11, numpy.array([ 1.        ,  1.20000005,  1.39999998,  1.60000002,  1.79999995,
		2.        ,  2.20000005,  2.4000001 ,  2.5999999 ,  2.79999995,
		3.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,  0.        ])),
	('Ecut', 0, 100.0, 10.0, 20, 20, 400, 400, 8, numpy.array([ 20.        ,  30,  40,  60,  100,
		140        ,  200,  400 ,  0 ,  0,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,  0.        ])),
], dtype=dtype)
# check that parameter definition agrees with computed models
assert numpy.product(parameters['NUMBVALS']) == len(table), ('parameter definition does not match spectra table', parameters['NUMBVALS'], numpy.product(parameters['NUMBVALS']), len(table))
hdu = pyfits.BinTableHDU(data=parameters)
hdu.header['DATE'] = nowstr
hdu.header['EXTNAME'] = 'PARAMETERS'
hdu.header['HDUCLAS1'] = 'XSPEC TABLE MODEL'
hdu.header['HDUVERS1'] = '1.0.0'
hdu.header['NINTPARM'] = len(parameters)
hdu.header['NADDPARM'] = 0
hdu.header['HDUCLAS2'] = 'PARAMETERS'
hdus.append(hdu)

# ENERG_LO, ENERG_HI
dtype = [('ENERG_LO', '>f4'), ('ENERG_HI', '>f4')]
energies = numpy.array(list(zip(energy_lo, energy_hi)), dtype=dtype)
hdu = pyfits.BinTableHDU(data=energies)
hdu.header['DATE'] = nowstr
hdu.header['EXTNAME']  = 'ENERGIES'
hdu.header['TUNIT2']   = 'keV'
hdu.header['TUNIT1']   = 'keV'
hdu.header['HDUCLAS1'] = 'XSPEC TABLE MODEL'
hdu.header['HDUCLAS2'] = 'ENERGIES'
hdu.header['HDUVERS1'] = '1.0.0'
hdus.append(hdu)

# PARAMVAL (4), INTPSPEC
dtype = [('PARAMVAL', '>f4', (len(parameters),)), ('INTPSPEC', '>f4', (nbins,))]
table.sort()
table = numpy.array(table, dtype=dtype)
hdu = pyfits.BinTableHDU(data=table)
hdu.header['DATE'] = nowstr
hdu.header['EXTNAME']  = 'SPECTRA'
hdu.header['TUNIT2']   = 'photons/cm^2/s'
hdu.header['TUNIT1']   = 'none'
hdu.header['HDUCLAS1'] = 'XSPEC TABLE MODEL'
hdu.header['HDUCLAS2'] = 'SPECTRA'
hdu.header['HDUVERS1'] = '1.0.0'

hdus.append(hdu)
hdus = pyfits.HDUList(hdus)

hdus.writeto(outfilename, overwrite=True)
