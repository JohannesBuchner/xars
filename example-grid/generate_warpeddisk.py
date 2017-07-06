import numpy
from numpy import exp, pi, log, log10, cos, sin
import sys
import h5py

#  rho:  density [M_sun pc^-3]
#  temp:  gas temperature [K]
#  temp_dust:  dust temperature [K]
#  u,v,w:   three components of velocity [pc/Myr]
#  nxh, nyh, nzh:  Grid size (x,y,z)
#  time:  time in Myr

nxh = 256
nyh = 256
nzh = 256
shape = (nxh, nyh, nzh)
center = numpy.array([nxh/2.+0.5, nyh/2.+0.5, nzh/2.+0.5])

pc_cm = 3.0856776e+18
nH_Msun = 1.18803e+57
CTNH = 10**25
CTrho = CTNH / nH_Msun * pc_cm**3 / (32 * pc_cm / 256)

def gamma(r):
	return r**0.5
def beta(r):
	return sin(gamma(r))/gamma(r)

fcov = 1
for diskfraction in 16, 4, 2, 1:
	rho = numpy.zeros(shape)

	# paint in the geometric shape
	for u in numpy.linspace(0, 2*pi, 401):
		v = numpy.linspace(1e-3, pi, 401)
		v = v[v < pi/diskfraction]
		x = v*( cos(u)*sin(gamma(v)) + sin(u)*cos(gamma(v))*cos(beta(v)))
		y = v*(-cos(u)*cos(gamma(v)) + sin(u)*sin(gamma(v))*cos(beta(v)))
		z = fcov * -v*sin(u)*sin(beta(v))
		
		i = (x / 6.3 * nxh + center[0]).astype(int)
		j = (y / 6.3 * nyh + center[1]).astype(int)
		k = (z / 6.3 * nzh + center[2]).astype(int)
		
		k -= 2
		rho[i,j,k] = CTrho
		k -= 1
		rho[i,j,k] = CTrho
	print diskfraction, rho[tuple(center.astype(int))]


	with h5py.File('warpeddisk_%s.hdf5' % diskfraction, 'w') as fout:
		d = fout.create_dataset('rho', data=rho, compression='gzip', shuffle=True)
		d.attrs['description'] = 'Density'
		d.attrs['unit'] = 'M_sun/pc^3'
		d = fout.create_dataset('nxh', data=nxh)
		d.attrs['description'] = 'number of grid cells in x axis'
		d = fout.create_dataset('nyh', data=nyh)
		d.attrs['description'] = 'number of grid cells in y axis'
		d = fout.create_dataset('nzh', data=nzh)
		d.attrs['description'] = 'number of grid cells in z axis'
		d = fout.create_dataset('center', data=center)
		d.attrs['unit'] = 'Myr'
	
		fout.attrs['CREATOR'] = """Johannes Buchner <johannes.buchner.acad@gmx.com>"""
		fout.attrs['ORIGIN'] = """Torqued disk with z-boost of %s""" % fcov




