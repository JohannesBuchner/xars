from __future__ import print_function, division
import numpy
from numpy import exp
import h5py
import astropy.io.fits as pyfits
import tqdm
from binning import nbins, bin2energy

energy_lo, energy_hi = bin2energy(numpy.arange(nbins))
energy = (energy_hi + energy_lo) / 2
deltae = energy_hi - energy_lo

table = []
PhoIndices = numpy.arange(1, 3.01, 0.1)
Ecuts = [ 40, 60, 100, 140, 200, 400 ]

data = {}

outfilename = 'blobs.fits'
rdataname = '%s_outreflectrdata.hdf5'
blobnhs = numpy.arange(22, 26.1, 1)
fcovs = numpy.arange(-5, 0.1, 1)
models = ['torusblob%.1f.hdf5' % blobnh for blobnh in blobnhs]
deltae0 = deltae[energy >= 1][0]

pbar = tqdm.tqdm(list(zip(blobnhs, models)))
for NHcloud, model in pbar:
	#print 'loading', model
	#m = h5py.File(model, 'r')
	#NHcloud = m['NH'][()]
	
	filename = rdataname % model
	f = h5py.File(filename, 'r')
	nphot = f.attrs['NPHOT']
	
	matrix = f['rdata']
	a, b, nmu = matrix.shape
	assert a == nbins, matrix.shape
	assert b == nbins, matrix.shape
	# sum over the viewing angles, but not the nh bins
	matrix_noinc = numpy.sum(matrix, axis=2)
	
	# go through viewing angles
	pbar.set_description('| NHblob=%.1f' % (NHcloud))
	matrix_mu = matrix_noinc * 1. / nphot
	for PhoIndex in PhoIndices:
		for Ecut in Ecuts:
			weights = (energy**-PhoIndex * exp(-energy / Ecut) * deltae / deltae0).reshape((-1,1))
			y = (weights * matrix_mu).sum(axis=0)
			table.append(((PhoIndex, Ecut, NHcloud), y))

hdus = []
hdu = pyfits.PrimaryHDU()
import datetime, time
now = datetime.datetime.fromtimestamp(time.time())
nowstr = now.isoformat()
nowstr = nowstr[:nowstr.rfind('.')]
hdu.header['CREATOR'] = """Johannes Buchner <johannes.buchner.acad@gmx.com>"""
hdu.header['DATE'] = nowstr
hdu.header['HDUCLAS1'] = 'XSPEC TABLE MODEL'
hdu.header['MODLNAME'] = 'blobreflect'
hdu.header['ADDMODEL'] = True
hdu.header['MODLUNIT'] = 'photons/cm^2/s'
hdu.header['EXTEND']   = True
hdu.header['REDSHIFT'] = True
hdu.header['SIMPLE']   = True
hdu.header['HDUDOC']   = 'OGIP/92-009'
hdu.header['HDUVERS1'] = '1.0.0'
hdu.header['HDUCLASS'] = 'OGIP'
hdus.append(hdu)

num_values = 41
# NAME, METHOD, INITIAL, DELTA, MINIMUM, BOTTOM, TOP, MAXIMUM, NUMBVALS, VALUE (41)
dtype = [('NAME', 'S12'), ('METHOD', '>i4'), ('INITIAL', '>f4'), ('DELTA', '>f4'), ('MINIMUM', '>f4'), ('BOTTOM', '>f4'), ('TOP', '>f4'), ('MAXIMUM', '>f4'), ('NUMBVALS', '>i4'), ('VALUE', '>f4', (num_values,))]

pyparameters = [
	('PhoIndex', 0, 2.0,   0.0099999998, 1.2, 2.8, PhoIndices),
	('Ecut',     0, 100.0, 10.0,        40, 400, Ecuts),
	('NH_blob',  0, 25.0,  0.5,         22, 26, blobnhs),
]
param_list = []
for paramname, a, default, step, lo, hi, values in pyparameters:
	padded_values = numpy.zeros(num_values)
	padded_values[:len(values)] = values
	param_list.append((paramname, a, default, step, numpy.min(values), lo, hi, numpy.max(values), len(values), padded_values))

parameters = numpy.array(param_list, dtype=dtype)
"""
parameters = numpy.array([
	('PhoIndex', 0, 2.0, 0.0099999998, 1.0, 1.2, 2.8, 3.0, 11, numpy.array([ 1.  ,  1.1, 1.2, 1.3, 1.4,  1.5, 1.6,  1.7, 1.8, 1.9, 2. , 2.1, 2.2,  2.3, 2.4,  2.5, 2.6 , 2.7, 2.8, 2.9, 3.,
                0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ])),
	('Ecut', 0, 100.0, 10.0, 40, 40, 400, 400, 6, numpy.array([ 40,  60,  100,
		140        ,  200,  400 ,  0 ,  0, 0, 0,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
		0.        ,  0.        ,  0.        ,  0.        ,  0.        ,  0.        ])),
	('NH_blob', 0, 25.0, 0.5, 22.0, 22, 26, 26, 21, numpy.array([ 
	         22.0, 22.2, 22.4, 22.6, 22.8, 
	         23.0, 23.2, 23.4, 23.6, 23.8, 
	         24.0, 24.2, 24.4, 24.6, 24.8, 
	         25.0, 25.2, 25.4, 25.6, 25.8, 
		 26.0      ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,
		 0.        ,   0.        ,   0.        ,   0.        ,   0.        ])),
], dtype=dtype)
"""

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


