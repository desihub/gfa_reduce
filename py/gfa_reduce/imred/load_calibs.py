# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.imred.load_calibs
============================

Image calibration utilities: overscan, bias, flat, etc.
"""
import gfa_reduce.common as common
import astropy.io.fits as fits
import os
import numpy as np
from desiutil.log import get_logger


def remove_overscan(image):
    sh = image.shape
    if sh[1] == 2248:
        _image = np.zeros((1032, 2048), dtype=float)
        _image[:, 0:1024] = image[:, 50:1074]
        _image[:, 1024:2048] = image[:, 1174:2198]
    return _image


def read_bias_image(extname):
    log = get_logger()
    assert(common.is_valid_extname(extname))

    par = common.gfa_misc_params()
    bias_fname = os.path.join(os.environ[par['meta_env_var']], \
                              par['master_bias_filename'])

    log.info('Attempting to read master bias : ' + bias_fname +
             ', extension name : ' + extname)

    assert(os.path.exists(bias_fname))

    bias = fits.getdata(bias_fname, extname=extname)

    bias = remove_overscan(bias)
    return bias

def read_flat_image(extname):
    log = get_logger()
    # at some point should add option to return master flat's
    # inverse variance as well
    assert(common.is_valid_extname(extname))

    par = common.gfa_misc_params()
    flat_fname = os.path.join(os.environ[par['meta_env_var']], \
                              par['master_flat_filename'])

    log.info('Attempting to read master flat : ' + flat_fname +
             ', extension name : ' + extname)

    assert(os.path.exists(flat_fname))

    flat = fits.getdata(flat_fname, extname=extname)

    flat = remove_overscan(flat)
    return flat

def read_static_mask_image(extname):
    log = get_logger()
    assert(common.is_valid_extname(extname))

    par = common.gfa_misc_params()
    mask_fname = os.path.join(os.environ[par['meta_env_var']], \
                              par['static_mask_filename'])

    log.info('Attempting to read static bad pixel mask : ' + mask_fname +
             ', extension name : ' + extname)

    assert(os.path.exists(mask_fname))

    mask = fits.getdata(mask_fname, extname=extname)

    mask = remove_overscan(mask)
    return mask
