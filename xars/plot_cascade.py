import numpy
import matplotlib.pyplot as plt
import sys
from weirdfunc import energy2bin, bin2energy

a = numpy.load(sys.argv[1])
nmu = a.shape[2]

ticks = [0.1,0.5,1, 2, 5, 10, 20, 50, 100]
ticks_labels = [str(t) for t in ticks]
ticks_values = energy2bin(ticks)

energy_hi, energy_lo = bin2energy(list(range(a.shape[0])))
weights = (energy_hi - energy_lo)
X, Y = numpy.meshgrid(weights, weights)
weights = X * Y
weights = 1

def maxoffdiagonal(m):
	m2 = m.copy()
	numpy.fill_diagonal(m2, 0)
	return m2.max()

view_matrices = [a[:,:,m] / weights for m in range(nmu)]

vmax = numpy.log10(max([maxoffdiagonal(m) for m in view_matrices]))
#vmax = numpy.max(view_matrices)

print([maxoffdiagonal(m) for m in view_matrices])
print([numpy.mean(m) for m in view_matrices])
print(vmax)
for b, m in zip(view_matrices, list(range(nmu))):
	#b =  # / 1e6
	#b /= a.max()
	#numpy.fill_diagonal(b, 0)

	plt.figure(figsize=(10,10))
	plt.title("angle = %.1f, peak: %.0f" % (
		numpy.arccos(m * 1. / nmu) * 180 / numpy.pi, b.max()))
	plt.imshow(numpy.log10(b), cmap=plt.cm.gray_r, vmin=0, vmax=vmax)
	plt.xticks(ticks_values, ticks_labels)
	plt.yticks(ticks_values, ticks_labels)
	plt.savefig(sys.argv[1] + "_%d.png" % m)



