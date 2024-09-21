import numpy
import scipy
from numpy import pi, tan, log10, sin, cos, logical_and

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
font = 'sans-serif'

from coordtrans import to_spherical, to_cartesian

def intersection_sort(intersections):
	#print 'sorting intersections ...'
	for i, (oa, da) in enumerate(intersections):
		assert da.dtype == numpy.float, da.dtype
		for j in range(i+1, len(intersections)):
			(ob, db) = intersections[j]
			assert db.dtype == numpy.float, db.dtype
			move_forward = db < da
			if not numpy.any(move_forward):
				continue
			#print '   %d: moving %d of %d forward from %d' % (i, move_forward.sum(), len(move_forward), j), 
			#print '       because:', da, db
			# swap
			tmp = numpy.copy(da[move_forward])
			da[move_forward] = db[move_forward]
			db[move_forward] = tmp
			tmp = numpy.copy(oa[move_forward])
			oa[move_forward] = ob[move_forward]
			ob[move_forward] = tmp
			#print '       now    :', da, db
	return intersections
			

class WedgeTorusGeometry(object):
	"""
	Spherical shell cut with two cones.
	"""
	def __init__(self, Theta_low, Theta_high, r_inner, r_outer, NH, verbose = False):
		validity, reason = WedgeTorusGeometry.is_valid(Theta_low, Theta_high, r_inner, r_outer, NH)
		assert validity, reason
		self.n = len(NH)
		self.Theta_high = Theta_high
		self.Theta_low = Theta_low
		self.NH = NH
		self.verbose = verbose
		self.r_inner = r_inner
		self.r_outer = r_outer
		self.density = numpy.asarray(NH).reshape((-1,))
		self.verbose = verbose
		self.plot = None
		for i in range(self.n):
			# looking through to source radially, the column density should be NH
			self.density[i] = self.NH[i] / (self.r_outer[i] - self.r_inner[i])
	
	@staticmethod
	def is_valid(Theta_low, Theta_high, r_inner, r_outer, NH):
		n = len(NH)
		if len(Theta_low) != n: return False, ['Theta_low length wrong', Theta_low]
		if len(Theta_high) != n: return False, ['Theta_high length wrong', Theta_high]
		if len(r_inner) != n: return False, ['r_inner length wrong', r_inner]
		if len(r_outer) != n: return False, ['r_outer length wrong', r_outer]
		for i in range(n):
			if not (Theta_low[i] >= 0): return False, ['Theta_low[%d] negative' % i, Theta_low[i]]
			if not (Theta_low[i] < Theta_high[i]): return False, ['component %d, Theta_low above Theta_high' % i, Theta_low[i], Theta_high[i]]
			if not (r_inner[i] < r_outer[i]): return False, ['component %d, r_inner above r_outer' % i, r_inner[i], r_outer[i]]
			# looking through to source radially, the column density should be NH
			for j in range(i):
				A = r_inner[i] > r_outer[j]
				B = r_inner[j] > r_outer[i]
				C = Theta_low[i] > Theta_high[j]
				D = Theta_low[j] > Theta_high[i]
				if not (A or B or C or D): return False, ['radial collision' if not (A or B) else '',
					'angle collision' if not (C or D) else '', 
					Theta_low, Theta_high, r_inner, r_outer]
		return True, None
		
	
	def compute_next_point(self, location, direction):
		# d = dist / self.NH # distance in units of nH
		(xi, yi, zi) = location
		(dist, beta, alpha) = direction
		(xi, yi, zi) = (
			numpy.asarray(xi, dtype=numpy.float).reshape((-1,)), 
			numpy.asarray(yi, dtype=numpy.float).reshape((-1,)), 
			numpy.asarray(zi, dtype=numpy.float).reshape((-1,)))
		(dist, beta, alpha) = (
			numpy.asarray(dist, dtype=numpy.float).reshape((-1,)), 
			numpy.asarray(beta, dtype=numpy.float).reshape((-1,)), 
			numpy.asarray(alpha, dtype=numpy.float).reshape((-1,)))
		
		radi, thetai, phii = to_spherical((xi, yi, zi))
		
		xv, yv, zv = to_cartesian((1, beta, alpha))
		if self.verbose: print('ray from', (xi, yi, zi), 'to', (xv, yv, zv))
		
		intersections = []
		
		if self.plot is not None:
			plt.plot(0.4 + 0.3*(xi + xv * numpy.linspace(-0.1, 5)), 
				 0.5 + 0.3*(zi + zv * numpy.linspace(-0.1, 5)), 
				 ':', color='orange')
			plt.plot(0.4 + 0.3*xi, 
				 0.5 + 0.3*zi, '*', color='orange')
		
			plt.savefig(self.plot)
		
		for i in range(self.n):
			# check if already inside
			already_inside = numpy.logical_and(
				numpy.logical_and(thetai >= self.Theta_low[i], thetai <= self.Theta_high[i]),
				numpy.logical_and(radi >= self.r_inner[i], radi <= self.r_outer[i]))
			if self.verbose: print('already inside %d' % i, already_inside)
			intersections.append((i * numpy.ones_like(already_inside, dtype=numpy.uint), numpy.where(already_inside, 0, numpy.nan)))
			
			# sphere intersections
			d1i, d2i = self.line_sphere_intersection((xi, yi, zi), (xv, yv, zv), self.r_inner[i])
			d1o, d2o = self.line_sphere_intersection((xi, yi, zi), (xv, yv, zv), self.r_outer[i])
			
			for inter in d1i, d2i, d1o, d2o:
				rf, thetaf, phif = to_spherical((xi+inter*xv, yi+inter*yv, zi+inter*zv))
				if self.verbose: print('   potential sphere surface intersection with %d:' % i, inter, (rf, thetaf, phif))
				relevant = numpy.logical_and(thetaf <= self.Theta_high[i], thetaf >= self.Theta_low[i])
				inter[~relevant] = numpy.nan
				if self.verbose: print('   sphere surface intersection with %d:' % i, inter, relevant)
				intersections.append((i * numpy.ones_like(inter, dtype=numpy.uint), inter))
				
			# cone intersections
			d1l, d2l = self.line_cone_intersection((xi, yi, zi), (xv, yv, zv), self.Theta_low[i])
			d1h, d2h = self.line_cone_intersection((xi, yi, zi), (xv, yv, zv), self.Theta_high[i])
			# is in sphere?
			for inter in d1l, d2l, d1h, d2h:
				rf, thetaf, phif = to_spherical((xi+inter*xv, yi+inter*yv, zi+inter*zv))
				if self.verbose: print('   potential cone surface intersection with %d:' % i, inter, (rf, thetaf, phif))
				relevant = numpy.logical_and(
					numpy.logical_and(rf <= self.r_outer[i], rf >= self.r_inner[i]),
					numpy.logical_and( thetaf - 1e-6 <= self.Theta_high[i], thetaf + 1e-6 >= self.Theta_low[i]))
				inter[~relevant] = numpy.nan
				if self.verbose: print('   cone surface intersection with %d:' % i, inter, relevant)
				intersections.append((i * numpy.ones_like(inter, dtype=numpy.uint), inter))
		
		# here we need something that goes through the intersections in the correct order, skipping nans
		# we have len(intersections) == 9 * self.n
		intersections = intersection_sort(intersections)
		
		if self.plot is not None:
			plt.plot([0.4 + 0.3 * (xi + d * xv) for i,d in intersections],
				[0.5 + 0.3 * (zi + d * zv)  for i,d in intersections],
				'+', color='orange', ms=5, alpha=0.9)
			plt.savefig(self.plot)
		# compute sequence of densities
		if self.verbose: 
			print('found intersections:')
			for oi, di in intersections:
				print('   ', oi, di)
		
		# find pairs of entry/exit, skipping nans
		#   find entry:
		last_d = numpy.ones_like(intersections[0][1]) * numpy.nan
		last_o = numpy.ones_like(intersections[0][1]) * numpy.nan
		terminate_inside_scene = numpy.ones_like(intersections[0][1]) == 0
		terminated = numpy.ones_like(intersections[0][1]) == 0
		rfinal = dist.copy()
		xf = xi.copy()
		yf = yi.copy()
		zf = zi.copy()
		
		if self.verbose: print('handling intersections; distance to go:', dist)
		for i in range(len(intersections)):
			oi, di = intersections[i]
			# if di is nan, do nothing
			# if di is not nan, and last_d is nan:
			#    this is an entry
			#    set last_d to di
			# if di is not nan, and last_d is not nan:
			#    this is an exit. reduce nh, check if end of trajectory
			#    oi must be last_o
			#    set last_d to nan
			point = -numpy.isnan(di)
			point_entry = numpy.logical_and(point,  numpy.isnan(last_d))
			point_exit  = numpy.logical_and(point, -numpy.isnan(last_d))
			last_d[point_entry] = di[point_entry]
			last_o[point_entry] = oi[point_entry]
			
			if not numpy.any(point_exit):
				continue
			assert numpy.all(last_o[point_exit] == oi[point_exit])
			
			distance = di - last_d
			if self.verbose: print('within %d for %f [%f .. %f]', (oi[point_exit], distance[point_exit], last_d[point_exit], di[point_exit]))
			
			if self.plot:
				plt.plot([0.4 + 0.3 * (xi[point_exit] + last_d[point_exit] * xv[point_exit]), 
					 0.4 + 0.3 * (xi[point_exit] + di * xv[point_exit])],
					[0.5 + 0.3 * (zi[point_exit] + last_d[point_exit] * zv[point_exit]), 
					 0.5 + 0.3 * (zi[point_exit] + di * zv[point_exit])],
					'+-', color='red', lw=4, alpha=0.3)
			
			nh = self.density[oi] * distance
			if self.verbose: print('NH: ', nh)
			
			assert dist.shape == nh.shape, (dist.shape, nh.shape)
			assert point_exit.shape == terminated.shape, (point_exit.shape, terminated.shape)
			assert dist.shape == terminated.shape, (dist.shape, terminated.shape)
			termination = numpy.logical_and(numpy.logical_and(dist <= nh, dist > 0), logical_and(point_exit, -terminated))
			gothrough   = numpy.logical_and(numpy.logical_and(dist > nh,  dist > 0), logical_and(point_exit, -terminated))
			t_dist = dist[termination]
			t_nh = nh[termination]
			t_distance = distance[termination]
			t_last_d = last_d[termination]
			t_this_d = di[termination]
			t_xi = xi[termination]
			t_yi = yi[termination]
			t_zi = zi[termination]
			t_xv = xv[termination]
			t_yv = yv[termination]
			t_zv = zv[termination]
			# so it has come to this: we interact
			assert numpy.all(t_dist > 0)
			assert numpy.all(t_nh > 0)
			assert numpy.all(t_distance > 0)
			assert numpy.all(t_last_d >= 0)
			d_interaction = t_dist / t_nh * t_distance + t_last_d
			assert numpy.all(d_interaction >= t_last_d)
			assert numpy.all(d_interaction <= t_this_d)
			
			if self.verbose: print('terminating %d, at' % termination.sum(), d_interaction)
			if numpy.any(termination):
				xf[termination] = t_xi + d_interaction * t_xv
				yf[termination] = t_yi + d_interaction * t_yv
				zf[termination] = t_zi + d_interaction * t_zv
				(t_rfinal, t_beta, t_alpha) = to_spherical((
					xf[termination], 
					yf[termination], 
					zf[termination]))
				rfinal[termination] = t_rfinal
				beta[termination] = t_beta
				alpha[termination] = t_alpha
				if self.plot: 
					plt.plot(0.4 + 0.3 * xf[termination], 0.5 + 0.3 * zf[termination],
						'p', color='orange', alpha=0.4)
				
				# we actually do terminate inside now
				terminate_inside_scene[termination] = True
				terminated[termination] = True
				if self.verbose: print('termination update:', terminate_inside_scene)
			
			dist[gothrough] -= nh[gothrough]
			if self.verbose: print('distance left:', dist[gothrough])
			last_d[point_exit] = numpy.nan
			last_o[point_exit] = numpy.nan
		
		# we leave the scene
		if self.verbose: print('termination inside:', terminate_inside_scene)
		return terminate_inside_scene, (xf,yf,zf), (rfinal, beta, alpha)
	
	def line_sphere_intersection(self, xxx_todo_changeme2, xxx_todo_changeme3, R):
		# compute sphere intersections
		# (xi + d * xv)**2 + (yi + d * yv)**2 + (zi + d * zv)**2 = R**2
		# 
		# xi**2 + 2 * xi * d * xv + d**2 + d**2 * xv**2 + 
		# yi**2 + 2 * yi * d * yv + d**2 + d**2 * yv**2 + 
		# zi**2 + 2 * zi * d * zv + d**2 + d**2 * zv**2 = R**2
		# 
		# a * x^2 + b * x + c = 0
		(xi, yi, zi) = xxx_todo_changeme2
		(xv, yv, zv) = xxx_todo_changeme3
		c = xi**2 + yi**2 + zi**2 - R**2
		b = 2 * (xi*xv + yi*yv + zi*zv)
		a = (xv**2 + yv**2 + zv**2)
		
		have_soln = b**2 >= 4 * a * c
		d1 = (-b - (b**2 - 4 * a * c)**0.5) / 2.
		d2 = (-b + (b**2 - 4 * a * c)**0.5) / 2.
		#xf = xi + d1 * xv
		#zf = zi + d1 * zv
		##if self.plot: plt.plot(0.4 + 0.3*xf, 0.5 + 0.3*zf, 'o', color='green' if d1 > 0 else 'red', alpha=0.3)
		#xf = xi + d2 * xv
		#zf = zi + d2 * zv
		##if self.plot: plt.plot(0.4 + 0.3*xf, 0.5 + 0.3*zf, 'o', color='green' if d2 > 0 else 'red', alpha=0.3)
		d1[d1 <= 0] = numpy.nan
		d2[d2 <= 0] = numpy.nan
		d1[~have_soln] = numpy.nan
		d2[~have_soln] = numpy.nan
		return d1, d2

	def line_cone_intersection(self, xxx_todo_changeme4, xxx_todo_changeme5, theta):
		# compute cone intersections
		# cos(theta) = (zi + d * zv) / ((xi + d * xv)**2 + (yi + d * yv)**2 + (zi + d * zv)**2)**0.5
		#assert theta >= 0, theta
		(xi, yi, zi) = xxx_todo_changeme4
		(xv, yv, zv) = xxx_todo_changeme5
		a = zv**2 - (xv**2 + yv**2)*tan(pi/2 - theta)**2
		b = 2.*zi*zv - (2.*xi*xv + 2.*yi*yv)*tan(pi/2 - theta)**2
		c = zi**2 - (xi**2 + yi**2)*tan(pi/2 - theta)**2
		quad = b**2 - 4.*a*c
		
		have_soln = quad >= 0
		d1 = (-b - quad**0.5)/(2. * a)
		d2 = (-b + quad**0.5)/(2. * a)
		if self.verbose: print('cone intersection distances', d1, d2)
		
		for di in d1, d2:
			xf = xi + di * xv
			zf = zi + di * zv
			dr, dtheta, dphi = to_spherical( (xi + di * xv, yi + di * yv, zi + di * zv))
			if self.plot: plt.plot(0.4 + 0.3*xf, 0.5 + 0.3*zf, '<', color='green' if di > 0 else 'red', alpha=0.3)
			di[di <= 0] = numpy.nan
			di[(dtheta - theta)**2 > 1e-9] = numpy.nan
		
		d1[~have_soln] = numpy.nan
		d2[~have_soln] = numpy.nan
		return d1, d2
		
	
	def viz_part(self, ax, s_Theta_low, s_Theta_high, NH, r_inner, r_outer):
		# high now refers to higher in z
		Theta_low  = s_Theta_low * 180 / pi
		Theta_high = s_Theta_high * 180 / pi
		nh = log10(NH) + 22
		thickness = max(0, min(1, (nh - 20.) / 5))
		
		ax.add_patch(  mpatches.Arc((0.4,0.5), 0.65, 0.65, theta2=90 - Theta_low, theta1=90 - Theta_high, 
			color='black', linewidth=1, alpha=1, fill=False, ls='dashed') )
		#plt.text(0.4 + 0.02, 0.5 + 0.25 + 0.02, "%2.0f -- %2.0f" % (Theta_low, Theta_high), ha="left", va='bottom',
		#	family=font, size=14)
		#ax.add_patch(  mpatches.Arc((0.4,0.5), 0.5, 0.5, theta1=0, theta2=Theta_low, 
		#	color='black', linewidth=1, alpha=1, fill=False, ls='dashed') )
		plt.text(0.4 + 0.05 + 0.3 * r_outer * cos(pi/2 - (s_Theta_low+s_Theta_high)/2.), 
			 0.5 + 0.02 + 0.3 * r_outer * sin(pi/2 - (s_Theta_low+s_Theta_high)/2.), 
			r"""$r=%2.2f-%2.2f$
$n_H=%2.1f$""" % (r_inner, r_outer, nh), ha="left", va='center',
			family=font, size=14)
		plt.text(0.4 + 0.05 + 0.3 * r_outer * cos(pi/2 - s_Theta_high), 
			 0.5 + 0.02 + 0.3 * r_outer * sin(pi/2 - s_Theta_high), 
			r"$%2.0f^\circ$" % (Theta_high), ha="left", va='bottom',
			family=font, size=14)
		plt.text(0.4 + 0.05 + 0.3 * r_outer * cos(pi/2 - s_Theta_low), 
			 0.5 + 0.02 + 0.3 * r_outer * sin(pi/2 - s_Theta_low), 
			r"$%2.0f^\circ$" % (Theta_low), ha="left", va='top',
			family=font, size=14)
		ax.add_patch( mpatches.Wedge(center=(0.4,0.5), r=0.3 * r_outer, 
			theta2=90 - Theta_low, 
			theta1=90 - Theta_high,
			width=0.3 * (r_outer - r_inner),
			color='black', 
			linewidth=0, alpha=thickness, fill=True) )
		ax.add_patch( mpatches.Wedge(center=(0.4,0.5), r=0.3 * r_outer, 
			theta1=90 + Theta_low, 
			theta2=90 + Theta_high,
			width=0.3 * (r_outer - r_inner),
			color='black', 
			linewidth=0, alpha=thickness, fill=True) )
		
		ax.add_patch( mpatches.Circle((0.4,0.5), 0.02, color='red', 
			linewidth=0, alpha=1, fill=True) )
	
	def viz(self): 
		""" Visualize the current geometry """
		plt.figure(figsize=(5,5))
		ax = plt.axes([0,0,1,1])
		ax.add_line(mlines.Line2D([0,0.9], [0.5,0.5], lw=1.,alpha=0.4, ls='dashed', color='grey'))
		ax.add_line(mlines.Line2D([0.4,0.4], [0.5,0.9], lw=1.,alpha=0.4, ls='dashed', color='grey'))
		ax.set_xticks([])
		ax.set_yticks([])
		
		for i in range(self.n):
			self.viz_part(ax, self.Theta_low[i], self.Theta_high[i], self.NH[i], self.r_inner[i], self.r_outer[i])
		
		plt.xlim(0,1)
		plt.ylim(0,1)

def test_1_wedge():
	geo = WedgeTorusGeometry(Theta_high=[(90--30) / 180. * pi],
		Theta_low=[(90-30) / 180. * pi], 
		r_inner = [0.4], r_outer=[0.9], NH=[10**(24-22)], verbose = True)
	geo.verbose = True
	geo.viz()
	plt.savefig('wedgetest.pdf')
	plt.close()
def test_2_wedges():
	geo = WedgeTorusGeometry(Theta_high=numpy.array([90--30, 90-40]) / 180. * pi,
		Theta_low=numpy.array([90-30, 90-90]) / 180. * pi, 
		r_inner = numpy.array([0.3, 0.1]), r_outer=numpy.array([0.4, 0.9]), NH=numpy.array([10**(24-22), 10**(21-22)]), verbose = True)
	geo.verbose = True
	geo.viz()
	plt.savefig('wedgetest2.pdf')
	plt.close()
def test_random_wedge():
	numpy.random.seed(0)
	for i in range(20):
		while True:
			d = dict(Theta_low=numpy.random.uniform(0,pi,2),
				NH=10**(numpy.random.uniform(20,25,2)-22),
				r_inner = numpy.random.uniform(0,1,2),
				verbose = True)
			d['Theta_high'] = numpy.random.uniform(d['Theta_low'],pi,2)
			d['r_outer'] = numpy.random.uniform(d['r_inner'],1,2)

			try:
				geo = WedgeTorusGeometry(**d)
				break
			except AssertionError as e:
				print(e)
				continue
		geo.viz()
		plt.savefig('wedgetestrand.pdf')
		plt.close()
def test_ray():
	geo = WedgeTorusGeometry(Theta_high=numpy.array([90, 45]) / 180. * pi,
		Theta_low=numpy.array([60, 0]) / 180. * pi, 
		r_inner = numpy.array([0.3, 0.1]), r_outer=numpy.array([0.4, 0.9]), NH=numpy.array([10**(24-22), 10**(21-22)]), verbose = True)
	# send in rays
	geo.plot = "raytest.pdf"
	geo.verbose = True
	geo.viz()
	(xi, yi, zi), (dist, beta, alpha) = (0.29,0.,-0.2), (10, pi/2 - pi/2 - 0.01, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	assert inside, inside
	assert zf < 0.3, zf
	assert yf == 0, yf
	assert xf < 0.3, xf
	assert r < 0.4, r
	assert r > 0.3, r
	
	(xi, yi, zi), (dist, beta, alpha) = (0.29,0.,-0.2), (78, pi/2 - pi/2 - 0.01, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	assert inside
	assert zf < 0.3, zf
	assert yf == 0
	assert xf < 0.3, xf
	assert r < 0.4
	assert r > 0.3

	(xi, yi, zi), (dist, beta, alpha) = (0.29,0.,-0.2), (78.44, pi/2 - pi/2 - 0.01, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	assert inside
	assert r < 0.9, r
	assert r > 0.3, r
	
	(xi, yi, zi), (dist, beta, alpha) = (0.29,0.,-0.2), (90, pi/2 - pi/2 - 0.01, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	assert not inside
	assert theta == beta, [theta, beta]
	assert phi == alpha, [phi, alpha]
	plt.close()
	
def test_ray_2():
	geo = WedgeTorusGeometry(
		Theta_high=numpy.array([95, 80]) / 180. * pi,
		Theta_low= numpy.array([85, 60]) / 180. * pi, 
		r_inner = numpy.array([0.01, 0.01]), r_outer=numpy.array([1, 1]), 
		NH=numpy.array([10**(24-22), 10**(21-22)]), verbose = True)
	# send in rays
	geo.plot = "raytest2.pdf"
	geo.verbose = True
	geo.viz()
	plt.savefig(geo.plot)
	(xi, yi, zi), (dist, beta, alpha) = (0.,0.,0.), (10, pi/2, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	assert inside
	assert r < 0.3
	assert r > 0.
	
	(xi, yi, zi), (dist, beta, alpha) = (0.,0.,0.), (50, pi/2, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	assert inside
	assert r < 0.6, r
	assert r > 0.
	assert theta == beta, [theta, beta]
	assert phi == alpha, [phi, alpha]

	
	(xi, yi, zi), (dist, beta, alpha) = (0.,0.,0.), (101, pi/2, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	assert not inside
	assert theta == beta, [theta, beta]
	assert phi == alpha, [phi, alpha]
	plt.close()
	

def test_ray_3():
	d = {
		'NH': [   3], 
		'r_inner': ([ 0.01]), 
		'r_outer': ([ 1]), 
		'Theta_low': ([ 0.]), 
		'Theta_high': ([ pi/8.])
	}
	geo = WedgeTorusGeometry(verbose = True, **d)
	# send in rays
	geo.plot = "raytest3.pdf"
	geo.verbose = True
	geo.viz()
	plt.savefig(geo.plot)
	(xi, yi, zi), (dist, beta, alpha) = (0.1,0,0.1), (0.01, pi + pi/2, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	plt.close()
	assert inside
	assert xf < xi, xf
	assert xf > 0., xf
	assert yf == 0, yf
	assert numpy.abs(zf-zi)/zi < 1e-6, zf
	
	geo.plot = "raytest3b.pdf"
	geo.viz()
	(xi, yi, zi), (dist, beta, alpha) = (-0.1,0,0.1), (0.01, pi/2, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	plt.close()
	assert inside
	assert xf < 0
	assert xf > xi
	assert yf == 0
	assert zf == zi

def test_ray_4():
	d = {
		'NH': ([ 47.45324503,  58.17402829]), 
		'Theta_low': ([ 2.19167304,  0.1892039 ]), 
		'r_outer': ([ 0.66491791,  0.51659453]), 
		'r_inner': ([ 0.21038256,  0.1289263 ]), 
		'Theta_high': ([ 2.49130462,  1.26301949])
	}
	d = {
		'NH': ([ 47]), 
		'Theta_low': ([ 2.19167304 ]), 
		'r_outer': ([ 0.66491791]), 
		'r_inner': ([ 0.21038256 ]), 
		'Theta_high': ([ 2.49130462])
	}
	
	geo = WedgeTorusGeometry(verbose = True, **d)
	# send in rays
	geo.plot = "raytest4.pdf"
	geo.verbose = True
	geo.viz()
	plt.savefig(geo.plot)
	(xi, yi, zi), (dist, beta, alpha) = (-0.02, 0, -0.2), (1, pi/2+pi/6, 0)
	(xi, yi, zi), (dist, beta, alpha) = (0.5, 0, 0.02), (1, -pi/6, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	plt.close()
	assert not inside

def test_ray_5():
	d = {
		'NH': ([ 47]), 
		'Theta_low': ([ pi/64 ]), 
		'r_outer': ([ 1]), 
		'r_inner': ([ 0.0 ]), 
		'Theta_high': ([ pi/32])
	}
	
	geo = WedgeTorusGeometry(verbose = True, **d)
	# send in rays
	geo.plot = "raytest5.pdf"
	geo.verbose = True
	geo.viz()
	plt.savefig(geo.plot)
	(xi, yi, zi), (dist, beta, alpha) = (-0.02, 0, -0.2), (1, pi/2+pi/6, 0)
	(xi, yi, zi), (dist, beta, alpha) = (0.2, 0, 0.7), (1, 2*pi-pi/4, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	plt.close()
	assert inside
	
	return
	geo.plot = "raytest4b.pdf"
	geo.viz()
	plt.savefig(geo.plot)
	(xi, yi, zi), (dist, beta, alpha) = (-0.02, 0, -0.2), (1, pi/2+pi/4, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	plt.close()
	assert inside

def test_ray_6():
	d = {
		'NH': ([ 100]), 
		'Theta_low': ([ pi/2 - pi/4 ]), 
		'r_outer': ([ 1]), 
		'r_inner': ([ 0.0 ]), 
		'Theta_high': ([ pi/2 + pi/100])
	}
	
	geo = WedgeTorusGeometry(verbose = True, **d)
	# send in rays
	geo.plot = "raytest6.pdf"
	geo.verbose = True
	geo.viz()
	plt.savefig(geo.plot)
	(xi, yi, zi), (dist, beta, alpha) = (0.2, 0, 0.5), (1, pi, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	plt.close()
	assert inside

def test_ray_7():
	d = {
		'NH': ([ 0.1]), 
		'Theta_low': ([ pi/2 - pi/4 ]), 
		'r_outer': ([ 1]), 
		'r_inner': ([ 0.0 ]), 
		'Theta_high': ([ pi/2 + pi/100])
	}
	
	geo = WedgeTorusGeometry(verbose = True, **d)
	# send in rays
	geo.plot = "raytest7.pdf"
	geo.verbose = True
	geo.viz()
	plt.savefig(geo.plot)
	(xi, yi, zi), (dist, beta, alpha) = (0.5, 0, 0.2), (0.1, - pi/2, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	plt.close()
	
	assert inside
	assert xf < -0.5, xf
	assert yf == 0, yf
	assert numpy.allclose(zf, 0.2), zf

def test_ray_8():
	d = {
		'NH': ([ 100 ]), 
		'Theta_low': ([ pi/4 ]), 
		'r_outer': ([ 1]), 
		'r_inner': ([ 0.0 ]), 
		'Theta_high': ([ pi/3])
	}
	
	geo = WedgeTorusGeometry(verbose = True, **d)
	# send in rays
	geo.plot = "raytest8.pdf"
	geo.verbose = True
	geo.viz()
	plt.savefig(geo.plot)
	(xi, yi, zi), (dist, beta, alpha) = (0., 0.2, 0.), (1, 0, 0)
	inside, (xf, yf, zf), (r, theta, phi) = geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
	plt.savefig(geo.plot)
	plt.close()
	
	assert inside
	assert xf == 0, xf
	assert yf == 0.2, yf
	assert zf > 0.01, zf
	assert zf < 0.4, zf


def speed_run():
	numpy.random.seed(0)
	plotting = False
	for i in range(50):
		while True:
			d = dict(Theta_low=numpy.random.uniform(0,pi,2),
				NH=10**(numpy.random.uniform(20,25,2)-22),
				r_inner = numpy.random.uniform(0,1,2),
				verbose = False)
			d['Theta_high'] = numpy.random.uniform(d['Theta_low'] + 1 / 180*pi,pi,2)
			d['r_outer'] = numpy.random.uniform(d['r_inner'] + 0.01,1,2)

			try:
				geo = WedgeTorusGeometry(**d)
				break
			except AssertionError as e:
				print(e)
				continue
		print(' --- setup %02d random geometry using:' % i, d)
		if plotting:
			geo.plot = "raytestrand_%02d.pdf" % i
			geo.viz()
			plt.savefig(geo.plot)
			plt.close()
		N = 100000
		(xi, yi, zi) = [numpy.random.uniform(-0.3, 0.3, N) for i in range(3)]
		dist = 10**numpy.random.uniform(-3, 3, N)
		beta  = numpy.random.uniform(0, pi, N)
		alpha = numpy.random.uniform(0, 2*pi, N)
		#yi *= 0
		#alpha *= 0
		if plotting:
			geo.viz()
		
		geo.compute_next_point((xi, yi, zi), (dist, beta, alpha))
		if plotting:
			plt.savefig(geo.plot)
			plt.close()

	
	
		
