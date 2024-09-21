from numpy import log10
from coordtrans import to_spherical, to_cartesian
import matplotlib.pyplot as plt

class SphereTorusGeometry(object):
	def __init__(self, NH, verbose = False):
		self.NH = NH
		self.verbose = verbose
	
	def compute_next_point(self, location, direction):
		(xi, yi, zi) = location
		(dist, beta, alpha) = direction
		d = dist / self.NH # distance in units of nH
		
		if self.verbose: print('  .. .. mean in nH units: ', d.mean())
		# compute relative vector traveled
		xv, yv, zv = to_cartesian((d, beta, alpha))
		
		# compute new position
		xf, yf, zf = xi + xv, yi + yv, zi + zv
		
		# compute spherical coordinates
		rad, phi, theta = to_spherical((xf, yf, zf))
		
		# are we inside the cone?
		inside = rad < 1.
		return inside, (xf,yf,zf), (rad, phi, theta)

	def viz(self): 
		""" Visualize the current geometry """
		nh = log10(self.NH) + 22
		
		import matplotlib.patches as mpatches
		import matplotlib.lines as mlines
		plt.figure(figsize=(5,5))
		font = 'sans-serif'
		ax = plt.axes([0,0,1,1])
		
		thickness = max(0, min(1, (nh - 20.) / 5))
		plt.text(0.35, 0.5,  "nH=%2.1f" % nh, ha="right", va='center',
			family=font, size=14)
		#print thickness
		ax.add_line(mlines.Line2D([0,0.9], [0.5,0.5], lw=1.,alpha=0.4, ls='dashed', color='grey'))
		ax.add_line(mlines.Line2D([0.4,0.4], [0.5,0.9], lw=1.,alpha=0.4, ls='dashed', color='grey'))
		ax.add_patch(  mpatches.Arc((0.4,0.5), 0.5, 0.5, theta2=90, theta1=90, 
			color='black', linewidth=1, alpha=1, fill=False, ls='dashed') )
		ax.add_patch( mpatches.Wedge((0.4,0.5), 0.3, -90, 90, color='black', 
			linewidth=0, alpha=thickness, fill=True) )
		ax.add_patch( mpatches.Wedge((0.4,0.5), 0.3, 90, -90, color='black', 
			linewidth=0, alpha=thickness, fill=True) )
	
		ax.add_patch( mpatches.Circle((0.4,0.5), 0.02, color='red', 
			linewidth=0, alpha=1, fill=True) )
		
		ax.set_xticks([])
		ax.set_yticks([])
	
