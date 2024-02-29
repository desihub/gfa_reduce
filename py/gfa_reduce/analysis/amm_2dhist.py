# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.analysis.amm_2dhist
==============================

Unknown purpose. Function is used by :mod:`~gfa_reduce.analysis.asterisms`.
"""
import numpy as np


def amm_2dhist(xmin, ymin, nx, ny, dx, dy, xvals, yvals):
    # xmin and ymin are the left and bottom EDGES of the lowest bins

    assert(len(xvals) == len(yvals))

    x_edges_left = xmin + np.arange(nx)*dx
    y_edges_left = ymin + np.arange(ny)*dy

    counts, _, __ = np.histogram2d(xvals, yvals, bins=[x_edges_left, y_edges_left])

    assert(np.sum(__ != y_edges_left) == 0)
    assert(np.sum(_ != x_edges_left) == 0)

    return counts, x_edges_left, y_edges_left
