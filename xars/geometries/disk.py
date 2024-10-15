from xars.coordtrans import to_cartesian, to_spherical


class DiskGeometry:
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.NH = 1

    def compute_next_point(self, location, direction):
        (xi, yi, zi) = location
        (dist, beta, alpha) = direction
        d = dist / self.NH  # distance in units of nH

        if self.verbose:
            print('  .. .. mean in nH units: ', d.mean())
        # compute relative vector traveled
        xv, yv, zv = to_cartesian((d, beta, alpha))

        # compute new position
        xf, yf, zf = xi + xv, yi + yv, zi + zv

        # compute spherical coordinates
        rad, phi, theta = to_spherical((xf, yf, zf))

        # are we inside the disk
        inside = zf < 0
        return inside, (xf,yf,zf), (rad, phi, theta)
