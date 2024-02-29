# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.analysis.radprof
===========================

Radial profiling.
"""
import numpy as np


def _atv_radplotf(img, x, y):

    # img should be background subtracted !!
    # x and y are the centroid coordinates of the profile

    sz = img.shape

    inrad = 0.5*np.sqrt(2)
    outrad = float(sz[0] // 2)
    drad = 1.0

    nrad = int(np.ceil((outrad - inrad)/drad) + 1)

    # initialize outputs
    radii = np.zeros(nrad, dtype='float32')
    profile = np.zeros(nrad, dtype='float32')

    distsq = np.zeros((sz[0], sz[1]), dtype=float)

    xx = np.arange(sz[1], dtype=float)
    yy = np.arange(sz[0], dtype=float)

    x2 = np.power(xx - x, 2)
    y2 = np.power(yy - y, 2)

    for i in range(sz[0]):
        distsq[i, :] = x2 + y2[i]

    ###import matplotlib.pyplot as plt

    ###plt.imshow(distsq, origin='lower', interpolation='nearest', cmap='gray')

    ###plt.show()

    for i in range(nrad):
        if i == 0:
            rin = 0.0
            rout = inrad
            rin2 = -0.01
        else:
            rin = inrad + drad*(i-1)
            rout = min(rin + drad, outrad)
            rin2 = rin*rin

        rout2 = rout*rout

        mask = (distsq > rin2) & (distsq <= rout2)

        _np = np.sum(mask)

        if _np > 0:
            radii[i] = (rout+rin)/2.0
            profile[i] = np.sum(img[mask])/_np
        else:
            radii[i] = rout

    return radii, profile
