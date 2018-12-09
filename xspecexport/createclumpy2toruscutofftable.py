from __future__ import print_function, division
import numpy
from numpy import pi, exp
import h5py
import astropy.io.fits as pyfits
import sys
import progressbar
from binning import nbins, energy2bin, bin2energy

energy_lo, energy_hi = bin2energy(numpy.arange(nbins))
energy = (energy_hi + energy_lo) / 2
deltae = energy_hi - energy_lo

table = []
PhoIndices = [ 1.        ,  1.20000005,  1.39999998,  1.60000002,  1.79999995,
		2.        ,  2.20000005,  2.4000001 ,  2.5999999 ,  2.79999995,
		3. ]
Ecuts = [ 20.,  30, 40, 60, 100, 140, 200, 400 ]

nh_bins = numpy.array([9.99999978e-03, 1.41000003e-02, 1.99999996e-02,
	2.82000005e-02,   3.97999994e-02,   5.62000014e-02,
	7.94000030e-02,   1.12000003e-01,   1.58000007e-01,
	2.24000007e-01,   3.16000015e-01,   4.46999997e-01,
	6.30999982e-01,   8.90999973e-01,   1.25999999e+00,
	1.77999997e+00,   2.50999999e+00,   3.54999995e+00,
	5.01000023e+00,   7.07999992e+00,   1.00000000e+01,
	1.41000004e+01,   2.00000000e+01,   2.82000008e+01,
	3.97999992e+01,   5.62000008e+01,   7.94000015e+01,
	1.12000000e+02,   1.58000000e+02,   2.24000000e+02,
	3.16000000e+02,   4.47000000e+02,   6.31000000e+02,
	8.91000000e+02,   1.26000000e+03,   1.78000000e+03,
	2.51000000e+03,   3.55000000e+03,   5.01000000e+03,
	7.08000000e+03,   1.00000000e+04])

data = {}

outfilename = sys.argv[1]
rdataname = sys.argv[2]
models = sys.argv[3:]
#sigmas  = [10,15,20,30,40,60,'sphere']
#sigmav = [10,15,20,30,40,60,90]
#Theta_tors = [80, 75, 70, 60, 50, 30,  0]
nh_bins_ThetaInc = [(nh, ThetaInc) for nh in nh_bins for ThetaInc in [90,60,0]]
deltae0 = deltae[energy >= 1][0]

widgets = [progressbar.Percentage(), " starting ... ", progressbar.Bar(), progressbar.ETA()]
pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(models)*len(nh_bins_ThetaInc)).start()


for model in models:
	#print 'loading', model
	m = h5py.File(model, 'r')
	sigma = m['sigma'].value * 2**0.5
	Ntot = m['Ntot'].value
	ThetaCloud = (Ntot / 100000.)**-0.5
	Theta_tor = 90 - sigma
	
	norm_filename = '%s_outnormalisation.fits' % model
	normalisations = pyfits.open(norm_filename)[0].data
	filename = rdataname % model # '_outrdata.hdf5'
	#print Ntot, sigma, ThetaCloud, Theta_tor, model
	f = h5py.File(filename, 'r')
	nphot = f.attrs['NPHOT']
	
	matrix = f['rdata']
	a, b, nmu = matrix.shape
	assert a == nbins, matrix.shape
	assert b == nbins, matrix.shape
	
	for mu, ((nh, ThetaInc), norm) in enumerate(zip(nh_bins_ThetaInc, normalisations)):
		# go through viewing angles
		matrix_mu = matrix[:,:,mu] * 1. / nphot / max(norm, 1e-6)
		#print '   ', nh, ThetaInc
		widgets[1] = '| op=%d cloud=%.1f nh=%.3f inc=%02d ' % (Theta_tor, ThetaCloud, nh, ThetaInc)
		pbar.update(pbar.currval + 1)
		for PhoIndex in PhoIndices:
			for Ecut in Ecuts:
				weights = (energy**-PhoIndex * exp(-energy / Ecut) * deltae / deltae0).reshape((-1,1))
				y = (weights * matrix_mu).sum(axis=0)
				table.append(((nh, PhoIndex, Ecut, Theta_tor, ThetaCloud, ThetaInc), y))
pbar.finish()

hdus = []
hdu = pyfits.PrimaryHDU()
import datetime, time
now = datetime.datetime.fromtimestamp(time.time())
nowstr = now.isoformat()
nowstr = nowstr[:nowstr.rfind('.')]
hdu.header['CREATOR'] = """Johannes Buchner <johannes.buchner.acad@gmx.com>"""
hdu.header['DATE'] = nowstr
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

# NAME, METHOD, INITIAL, DELTA, MINIMUM, BOTTOM, TOP, MAXIMUM, NUMBVALS, VALUE (41)
dtype = [('NAME', 'S12'), ('METHOD', '>i4'), ('INITIAL', '>f4'), ('DELTA', '>f4'), ('MINIMUM', '>f4'), ('BOTTOM', '>f4'), ('TOP', '>f4'), ('MAXIMUM', '>f4'), ('NUMBVALS', '>i4'), ('VALUE', '>f4', (41,))]

parameters = numpy.array([
	('nH', 1, 10.0, 1.0, 0.0099999998, 0.0099999998, 10000.0, 10000.0, 41, numpy.array([  9.99999978e-03,   1.41000003e-02,   1.99999996e-02,
		 2.82000005e-02,   3.97999994e-02,   5.62000014e-02,
		 7.94000030e-02,   1.12000003e-01,   1.58000007e-01,
		 2.24000007e-01,   3.16000015e-01,   4.46999997e-01,
		 6.30999982e-01,   8.90999973e-01,   1.25999999e+00,
		 1.77999997e+00,   2.50999999e+00,   3.54999995e+00,
		 5.01000023e+00,   7.07999992e+00,   1.00000000e+01,
		 1.41000004e+01,   2.00000000e+01,   2.82000008e+01,
		 3.97999992e+01,   5.62000008e+01,   7.94000015e+01,
		 1.12000000e+02,   1.58000000e+02,   2.24000000e+02,
		 3.16000000e+02,   4.47000000e+02,   6.31000000e+02,
		 8.91000000e+02,   1.26000000e+03,   1.78000000e+03,
		 2.51000000e+03,   3.55000000e+03,   5.01000000e+03,
		 7.08000000e+03,   1.00000000e+04])),
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
	('Theta_tor', 0, 50.0, 5.0, 0.0, 30, 85, 90.0, 4, numpy.array([ 30,  60,  80, 85, 
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,   0.        ])),
	('Theta_cloud', 0, 1.0, 0.3, 1.0, 1., 5.77350269, 5.77350269, 4, numpy.array([ 1. ,  1.82574186, 3.16227766,  5.77350269  ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,   0.        ])),
	('Theta_inc', 0, 18.200001, 5.0, 0.0, 18.200001, 87.099998, 90.0, 3, numpy.array([ 0,  60,  90,  0.0       ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,   0.        ])),
], dtype=dtype)
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

hdus.writeto(outfilename, clobber=True)


