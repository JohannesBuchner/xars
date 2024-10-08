#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
from numpy import cos, log

from . import conetorus

"""
      ⎛         za                   ⎞
      ⎜  -NH₀ - ──                   ⎟
      ⎜         zg                   ⎟
      ⎜10         ⋅log(10)⋅cos(θ)    ⎟
zg⋅log⎜────────────────────────── + 1⎟
      ⎝           d⋅zg               ⎠
────────────────────────────────────── = dist
            log(10)⋅cos(θ)

NH = 1/dist
"""


class GradientConeTorusGeometry(conetorus.ConeTorusGeometry):
    def __init__(self, Theta_tor, NH0=21 - 22, z0=0, zg=-1, verbose=False):
        self.NH0 = NH0
        self.z0 = z0
        self.zg = zg
        self.verbose = verbose
        super(GradientConeTorusGeometry, self).__init__(Theta_tor, 1)

    def compute_next_point(self, location, direction):
        # print 'beta, zg', beta, self.zg
        (xi, yi, zi) = location
        (dist, beta, alpha) = direction
        f = log(10) * cos(beta) / self.zg
        # print 'f, dist', f, dist
        logval = f / dist * 10**(-self.NH0 - (zi + self.z0) / self.zg) + 1
        # print 'logval', logval

        self.NH = numpy.where(logval < 0, 1e-100, f / log(logval))
        assert numpy.all(self.NH > 0), [self.NH[~(self.NH > 0)], f, logval]
        return super(GradientConeTorusGeometry, self).compute_next_point(
            (xi, yi, zi), (dist, beta, alpha))
