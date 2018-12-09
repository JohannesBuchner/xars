from __future__ import print_function, division
import numpy
from numpy import pi, exp
import astropy.io.fits as pyfits
import h5py
import sys
from binning import nbins, energy2bin, bin2energy

energy_lo, energy_hi = bin2energy(numpy.arange(nbins))
energy = (energy_hi + energy_lo) / 2
deltae = energy_hi - energy_lo

table = []
PhoIndices = [ 1.        ,  1.20000005,  1.39999998,  1.60000002,  1.79999995,
		2.        ,  2.20000005,  2.4000001 ,  2.5999999 ,  2.79999995,
		3. ]
ThetaIncs = [ 18.20000076,  31.79999924,  41.40000153,  49.5       ,
		56.59999847,  63.29999924,  69.5       ,  75.5       ,
		81.40000153,  87.09999847]
Ecuts = [ 20.,  30, 40, 60, 100, 140, 200, 400 ]
data = {}

outfilename = sys.argv[1]

for filename in sys.argv[2:]:
	print('loading', filename)
	#f = pyfits.open(filename)
	#nphot = int(f[0].header['NPHOT'])
	#matrix = f[0].data
	f = h5py.File(filename)
	nphot = f.attrs['NPHOT']
	matrix = f['rdata']
	a, b, nmu = matrix.shape
	assert a == nbins, f[0].data.shape
	assert b == nbins, f[0].data.shape
	#data[(nh, opening)] = [(nphot, f[0].data)]
	
	for mu, ThetaInc in enumerate(ThetaIncs[::-1]):
		# go through viewing angles
		matrix_mu = matrix[:,:,mu]
		print(ThetaInc)
		for PhoIndex in PhoIndices:
			for Ecut in Ecuts:
				weights = (energy**-PhoIndex * exp(-energy / Ecut) * deltae).reshape((-1,1))
				y = (weights * matrix_mu).sum(axis=0) / (nphot / 10.)
				#print PhoIndex, ThetaInc #, (y/deltae)[energy_lo >= 1][0]
				#print '    ', (weights * matrix[:,:,mu]).sum(axis=0), deltae, (nphot / 1000000.)
				#assert numpy.any(y > 0), y
				table.append(((PhoIndex, Ecut, ThetaInc), y))

hdus = []
hdu = pyfits.PrimaryHDU()
import datetime, time
now = datetime.datetime.fromtimestamp(time.time())
nowstr = now.isoformat()
nowstr = nowstr[:nowstr.rfind('.')]
hdu.header['CREATOR'] = """Johannes Buchner <johannes.buchner.acad@gmx.com>"""
hdu.header['DATE'] = nowstr
hdu.header['HDUCLAS1'] = 'XSPEC TABLE MODEL'
hdu.header['MODLNAME'] = 'disk'
hdu.header['ADDMODEL'] = True
hdu.header['MODLUNIT'] = 'photons/cm^2/s'
hdu.header['EXTEND']   = True
hdu.header['REDSHIFT'] = True
hdu.header['SIMPLE']   = True
hdu.header['HDUDOC']   = 'OGIP/92-009'
hdu.header['HDUVERS1'] = '1.0.0'
hdu.header['HDUCLASS'] = 'OGIP'
hdus.append(hdu)

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
	('Theta_inc', 0, 18.200001, 5.0, 0.0, 18.200001, 87.099998, 90.0, 10, numpy.array([ 18.20000076,  31.79999924,  41.40000153,  49.5       ,
		56.59999847,  63.29999924,  69.5       ,  75.5       ,
		81.40000153,  87.09999847,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,   0.        ])),
], dtype=dtype)
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

hdus.writeto(outfilename, clobber=True)


