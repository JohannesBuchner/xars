import numpy
from numpy import exp
import h5py
import progressbar
from xars.tableexport import write_atable

table = []
PhoIndices = [ 1.        ,  1.20000005,  1.39999998,  1.60000002,  1.79999995,
		2.        ,  2.20000005,  2.4000001 ,  2.5999999 ,  2.79999995,
		3., 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5.0 ]
Ecuts = [ 1, 2, 3, 5., 10., 15, 20., 25, 30, 35, 40, 50, 60, 80, 100, 140, 160, 200, 300, 400 ]

outfilename = 'blobs.fits'
rdataname = '%s_outreflectrdata.hdf5'
blobnhs = numpy.round(numpy.arange(22, 26.1, 0.2), 8)
models = ['torusblob%.1f.hdf5' % blobnh for blobnh in blobnhs]

widgets = [progressbar.Percentage(), " starting ... ", progressbar.Bar(), progressbar.ETA()]
pbar = progressbar.ProgressBar(widgets=widgets)

for NHcloud, model in zip(pbar(blobnhs), models):
	filename = rdataname % model
	f = h5py.File(filename, 'r')
	energy_lo = f['energy_lo'][()]
	energy_hi = f['energy_hi'][()]
	nbins = len(energy_lo)
	energy = (energy_hi + energy_lo) / 2
	deltae = energy_hi - energy_lo
	deltae0 = deltae[energy >= 1][0]
	nphot = f.attrs['NPHOT']
	
	matrix = f['rdata']
	a, b, nmu = matrix.shape
	assert a == nbins, matrix.shape
	assert b == nbins, matrix.shape
	# sum over the viewing angles, but not the nh bins
	matrix_noinc = numpy.sum(matrix, axis=2)
	
	# go through viewing angles
	widgets[1] = '| NHblob=%.1f' % (NHcloud)
	matrix_mu = matrix_noinc * 1. / nphot
	for PhoIndex in PhoIndices:
		for Ecut in Ecuts:
			weights = energy**-PhoIndex * exp(-energy / Ecut) * deltae / deltae0
			y = weights @ matrix_mu
			table.append(((PhoIndex, Ecut, NHcloud), y))

energy_info = dict(energy_lo=energy_lo, energy_hi=energy_hi)

write_atable(outfilename,
	parameter_definitions = [
		('PhoIndex', 0, 2.0, 0.0099999998, 1.2, 2.8, PhoIndices),
		('Ecut', 0, 100.0, 10.0, 40, 400, Ecuts),
		('NH_blob', 0, 25.0, 0.5, 22, 26, blobnhs), 
	], table=table, modelname='blobreflect', **energy_info)

