# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.analysis.basic_catalog_stats
=======================================

Determine image summary statistics.
"""
import gfa_reduce.analysis.util as util
import numpy as np


# do a meta-analysis of the per-object FWHM values,
# presumably restricting to relative high s/n sources
# if no good sources are available, then need to return some dummy
# value
def overall_image_fwhm(tab, bad_amps=None, snr_thresh=20):
    # remove bad/dummy FWHM values
    # remove sources that are near image boundaries
    # could remove sources that aren't isolated, but might be fine
    # to leave this as an optimization for the future

    assert(len(np.unique(tab['camera'])) == 1)

    good = util.use_for_fwhm_meas(tab, bad_amps=bad_amps, snr_thresh=snr_thresh)

    # check for case of not enough sources
    ngood = np.sum(good)
    if ngood == 0:
        return -1, -1, -1, -1, 0, good

    sig_to_fwhm = 2.355
    # just take geometric mean of semi-major / semi-minor sigmas ??
    if ngood > 1:
        fwhm_major_pix = np.median(tab[good]['sig_major_pix'])*sig_to_fwhm
        fwhm_minor_pix = np.median(tab[good]['sig_minor_pix'])*sig_to_fwhm
    else:
        fwhm_major_pix = tab[good][0]['sig_major_pix']*sig_to_fwhm
        fwhm_minor_pix = tab[good][0]['sig_minor_pix']*sig_to_fwhm

    fwhm_pix = np.sqrt(fwhm_major_pix*fwhm_minor_pix)

    # this is a temporary HACK -- should really take into account
    # the PA and directionality of the platescale when doing this

    asec_per_pix = util.nominal_pixel_sidelen_arith()

    fwhm_asec= asec_per_pix*fwhm_pix

    return fwhm_major_pix, fwhm_minor_pix, fwhm_pix, fwhm_asec, int(np.sum(good)), good
