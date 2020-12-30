import gfa_reduce.common as common
import numpy as np
import gfa_reduce.analysis.util as util

def adu_to_surface_brightness(sky_adu_1pixel, acttime, extname):
    """
    convert from ADU (per pixel) to mag per square asec (AB)

    note that this is meant to be applied to an average sky value across
    an entire GFA camera; this function does not take into account
    platescale variations within a camera
    """

    if (sky_adu_1pixel <= 0) or (acttime <= 0):
        return np.nan

    par = common.gfa_misc_params()

    pixel_area_sq_asec = util.nominal_pixel_area_sq_asec(extname)

    sky_adu_per_sq_asec = sky_adu_1pixel/pixel_area_sq_asec

    sky_adu_per_sec_sq_asec = sky_adu_per_sq_asec/acttime

    zp_adu_per_s = util.median_zenith_camera_zeropoint(extname)

    return (zp_adu_per_s - 2.5*np.log10(sky_adu_per_sec_sq_asec))
    
