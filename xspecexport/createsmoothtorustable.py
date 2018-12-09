from __future__ import print_function, division
import numpy
from numpy import pi, exp
import astropy.io.fits as pyfits
import h5py
import sys
import progressbar
from binning import nbins, energy2bin, bin2energy

energy_lo, energy_hi = bin2energy(numpy.arange(nbins))
energy = (energy_hi + energy_lo) / 2
deltae = energy_hi - energy_lo
deltae0 = deltae[energy >= 1][0]

table = []
PhoIndices = [ 1.        ,  1.20000005,  1.39999998,  1.60000002,  1.79999995,
		2.        ,  2.20000005,  2.4000001 ,  2.5999999 ,  2.79999995,
		3. ]
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
prefix = sys.argv[2]
sigmas  = ['5_gexp2','10_gexp2','20_gexp2','30_gexp2','sphere']
sigmav = [5,10,20,30,90]
opening = [90-s for s in sigmav]
filenames = [prefix % o for o in sigmas]

widgets = [progressbar.Percentage(), " starting ... ", progressbar.Bar(), progressbar.ETA()]
pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(filenames)*len(nh_bins)).start()


for Theta_tor, filename in zip(opening, filenames):
	#print 'loading', filename
	f = h5py.File(filename)
	nphot = f.attrs['NPHOT']
	
	geometry = numpy.loadtxt(filename.replace('_rdata.hdf5', '').replace('_transmitrdata.hdf5', '').replace('_reflectrdata.hdf5', '').replace('layered',''))
	
	matrix = f['rdata']
	a, b, nmu = matrix.shape
	assert a == nbins, matrix.shape
	assert b == nbins, matrix.shape
	#data[(nh, opening)] = [(nphot, f[0].data)]
	
	last_angle = 0
	# go through viewing angles
	for mu, nh in enumerate(nh_bins):
		angle = geometry[mu,0]
		# theta = arccos(1-x)
		deltac = numpy.cos(last_angle) - numpy.cos(angle)
		last_angle = angle
		
		matrix_mu = matrix[:,:,mu]
		widgets[1] = '| op=%d nh=%.3f ' % (Theta_tor, nh)
		pbar.update(pbar.currval + 1)
		for PhoIndex in PhoIndices:
			weights = (energy**-PhoIndex * deltae).reshape((-1,1))
			y = (weights * matrix_mu).sum(axis=0) / nphot * 1. / deltac
			#print nh, PhoIndex, Theta_tor #, (y/deltae)[energy_lo >= 1][0]
			#print '    ', (weights * matrix[:,:,mu]).sum(axis=0), deltae, (nphot / 1000000.)
			#assert numpy.any(y > 0), y
			table.append(((nh, PhoIndex, Theta_tor), y))
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
	('Theta_tor', 0, 60.0, 5.0, 0.0, 0.0, 85.0, 85.0, 5, numpy.array([ 0,60,70,80,85,
		             0.        ,   0.        ,   0.        ,
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


