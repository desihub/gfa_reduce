import ci_reduce.common as common
import numpy as np

def adu_to_surface_brightness(sky_adu_1pixel, acttime, extname):
    """
    convert from ADU (per pixel) to mag per square asec (AB)
    """

    assert(sky_adu_1pixel > 0)
    assert(acttime > 0)

    par = common.ci_misc_params()

    if extname is 'CIC':
        pixel_area_sq_asec = (par['nominal_cen_cd']*3600.0)**2
    else:
        pixel_area_sq_asec = \
            (par['nominal_mer_cd']*3600.0)*(par['nominal_sag_cd']*3600.0)

    sky_adu_per_sq_asec = sky_adu_1pixel/pixel_area_sq_asec

    sky_adu_per_sec_sq_asec = sky_adu_per_sq_asec/acttime

    sky_e_per_sec_sq_asec = sky_adu_per_sec_sq_asec*common.ci_camera_gain(extname)

    return (par['nominal_zeropoint'] - 2.5*np.log10(sky_e_per_sec_sq_asec))
    