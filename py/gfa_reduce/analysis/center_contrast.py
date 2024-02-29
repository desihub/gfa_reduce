# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.analysis.center_contrast
===================================

Determine flux ratios.
"""
import numpy as np
import copy

def peak_flux(density_image):
    # probably going to have issues if image isn't square, don't
    # worry about that for now

    sz = density_image.shape
    nx = sz[0]
    ny = sz[1]

    indmax = np.argmax(density_image)

    # this is different from the IDL version because of the way Python
    # indexing works !!!!
    iy = (indmax % ny).astype(int)
    ix = (indmax // ny).astype(int)

    assert(density_image[ix, iy] == np.max(density_image))

    half = 3

    xmin = max(ix - half, 0)
    xmax = min(ix + half + 1, nx)

    ymin = max(iy - half, 0)
    ymax = min(iy + half + 1, ny)

    cutout = density_image[xmin:xmax, ymin:ymax]

    return np.sum(cutout), ix, iy, nx, ny

def center_contrast(density_image):
    peak_flux_best, ix, iy, nx, ny = peak_flux(density_image)

    _density_image = copy.deepcopy(density_image)

    # may need tweaking for GFA's !!!!
    half = 10 # ~2 asec for the GFA's

    xmin = max(ix - half, 0)
    xmax = min(ix + half + 1, nx)
    ymin = max(iy - half, 0)
    ymax = min(iy + half + 1, ny)

    _density_image[xmin:xmax, ymin:ymax] = 0.0
    peak_flux_second_best, _, __, ___, ____ = peak_flux(_density_image)

    # will the denominator ever be zero ???
    contrast = peak_flux_best/peak_flux_second_best

    return contrast
