# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.analysis.segment
===========================

Segment a 2D numpy array.
"""
import numpy as np
import gfa_reduce.common as common
import gfa_reduce.analysis.util as util
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
try:
    # photutils >= 2.0
    from photutils.segmentation import detect_threshold, detect_sources
except ImportError:
    # photutils < 2.0
    from photutils import detect_threshold, detect_sources


def segmentation_map(image, extname, get_kernel=False):
    """Perform segmentation on `image`.

    Parameters
    ----------
    image : :class:`numpy.ndarray`
        A 2D NumPy array rather than a GFA_image object.
    extname : :class:`str`
        This parameter does not appear to be used.
    get_kernel : :class:`bool`, optional
        If ``True``, return the kernel used in segmentation.

    Returns
    -------
    :class:`~photutils.segmentation.SegmentationImage`
        A segmentation image.
    """

    par = common.gfa_misc_params()

    fwhm_pix = par['nominal_fwhm_asec'] / \
        util.nominal_pixel_sidelen_arith()

    threshold = detect_threshold(image, snr=2.0)

    sigma = fwhm_pix*gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=int(np.round(fwhm_pix)),
                                     y_size=int(np.round(fwhm_pix)))
    kernel.normalize()

    segm = detect_sources(image, threshold, npixels=5, filter_kernel=kernel)

    # add my own dilation of segm.array ?
    # incorporate masking based on master flat/bias in this analysis ?

    if not get_kernel:
        return segm
    else:
        return segm, kernel
