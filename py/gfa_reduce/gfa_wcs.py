# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.gfa_wcs
==================

WCS-related utilities.
"""
import astropy.io.fits as fits
import numpy as np
from astropy import wcs
import os
import gfa_reduce.common as common
import pickle

def local_tan_wcs(telra, teldec, extname):
    wcs_big = nominal_tan_wcs(telra, teldec, extname)

    crval1, crval2 = ccd_center_radec(wcs_big)

    par = common.gfa_misc_params()

    fname = os.path.join(os.environ[par['meta_env_var']],
                         'dummy_with_headers_local_SIP.zenith.fits.gz')

    h = fits.getheader(fname, extname=extname)

    h['CRVAL1'] = float(crval1)
    h['CRVAL2'] = float(crval2)

    w = wcs.WCS(h)

    return w

def nominal_tan_wcs(telra, teldec, extname):
    # Create a new WCS object.  The number of axes must be set
    # from the start

    par = common.gfa_misc_params()

    fname = os.path.join(os.environ[par['meta_env_var']],
                         par['wcs_templates_filename'])

    templates = pickle.load(open(fname, 'rb'))

    h = templates[extname]

    h['CRVAL1'] = telra
    h['CRVAL2'] = teldec

    w = wcs.WCS(h)

    return w

def ccd_center_radec(_wcs):
    y_pixel_center = 515.5 # need to double check this relative to convention
    x_pixel_center = 1023.5 # need to double check this relative to convention

    ra, dec = _wcs.all_pix2world(x_pixel_center, y_pixel_center, 0)

    return ra, dec
