# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.imred.dq_mask
========================

Image mask utilities.
"""
import gfa_reduce.common as common
import gfa_reduce.imred.load_calibs as load_calibs
import numpy as np


def create_satmask(im, extname):
    # im is just a 2D array of pixels, not a GFA_image object

    par = common.gfa_misc_params()

    gain = common.gfa_camera_gain(extname)

    sat_thresh = par['full_well_electrons']/gain

    satmask = (im >= sat_thresh)

    return satmask

def create_nanmask(im):
    # im is just a 2D array of pixels, not a GFA_image object

    return np.logical_not(np.isfinite(im))

def dq_bitmask(im, extname):
    """
    bit definitions
    2^0 = bad pixels based on master flat (inherit from static bad pixel mask)
    2^1 = questionable based on master flat (inherit from static bad pixel mask)
    2^2 = saturated in image
    2^3 = non-finite values in image
    """

    # im is just a 2D array of pixels, not a GFA_image object

    mask = load_calibs.read_static_mask_image(extname)

    ###maskbits = common.mask_bit_dict()

    ###satmask = create_satmask(im, extname)
    ###mask = np.bitwise_or(mask, (2**maskbits['SATUR'])*satmask)

    ###nanmask = create_nanmask(im)
    ###mask = np.bitwise_or(mask, (2**maskbits['NAN'])*nanmask)

    return mask

def add_dq_bitmask_header_cards(h):
    # h should be a FITS header
    # input gets modified !!

    maskbits = common.mask_bit_dict()

    for k,v in maskbits.items():
        card_name = 'MASKB' + str(v)
        card_value = common.mask_bit_description(k)

        comment = 'Mask bit 2**' + str(v) + '=' + str(2**v) + ' meaning'
        h[card_name] = (card_value, comment)

    return h
