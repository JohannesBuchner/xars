import sys
import h5py

nphot = 0
rdata = 0
for a in sys.argv[1:-1]:
	
	print 'accumulating', a
	with h5py.File(a, 'r') as fa:
		date = fa.attrs['DATE']
		method = fa.attrs['METHOD']
		creator = fa.attrs['CREATOR']
		nphot = nphot + fa.attrs['NPHOT']
		rdata = rdata + fa['rdata'].value
		print 'total of %d input / %d output photons' % (nphot, rdata.sum())


c = sys.argv[-1]
print 'writing to', c
with h5py.File(c, 'w') as f:
	f.create_dataset('rdata', data=rdata, compression='gzip', shuffle=True)
	f.attrs['CREATOR'] = creator
	f.attrs['DATE'] = date
	f.attrs['METHOD'] = method
	f.attrs['NPHOT'] = nphot


