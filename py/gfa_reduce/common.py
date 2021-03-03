import numpy as np
import os

# could add utilities re: rotation relative to CS5

def gfa_misc_params():
    """
    repository for various constant values that I don't want to 
    end up hardcoding in various places throughout the code base
    """

    width_pix_native = 2048
    height_pix_native = 1032

    par = {'meta_env_var': 'GFA_REDUCE_META',
           'gaia_env_var': 'GAIA_CAT_DIR',
           'ps1_env_var' : 'PS_CAT_DIR',
           'width_pix_native': width_pix_native,
           'height_pix_native': height_pix_native,
           'width_with_prescan_overscan' : 2248,
           'height_with_prescan_overscan' : 1032,
           'n_cameras': 6,
           'nominal_zeropoint': 27.0621,
           'gfa_exp_extname': 'GFA',
           'guider_exp_extname' : 'GUIDER',
           'master_bias_filename': 'GFA_master_bias-overscan_subtracted.fits', 
           'master_flat_filename': 'flat_LED_5s_20200115-all.fits',
           'master_dark_filename' : 'master_dark_library/master_dark-00026584_00026684.fits',
           'dark_index_filename' : 'master_dark_index-20200128.fits',
           'static_mask_filename': 'GFA_static_badpixels.fits',
           'nominal_sag_cd': (5.0/3.0)*3.55978e-5, 
           'nominal_mer_cd': (5.0/3.0)*3.26627e-5,
           'nominal_cen_cd': (5.0/3.0)*3.70303e-5,
           'full_well_electrons' : 100000.0,
           'sat_thresh_adu' : 40000.0,
           'nominal_fwhm_asec' : 1.25,
           'wcs_templates_filename' : 'dummy_with_headers_gaia.zenith.pkl',
           'reduced_image_flavors' : ['REDUCED', 'INVVAR', 'BITMASK', 'DETMAP'],
           'kterm' : 0.114,
           'zp_filename' : 'dense_field_zeropoints-combined.all_cameras.fits',
           'ephem_filename' : 'gfa_ephemeris.fits',
           'kpno_lat_deg' : 31.9639671,
           'exp_kernel_filename' : 'elg_exponential.fits',
           'devauc_kernel_filename' : 'bgs_devauc.fits',
           'fracflux_nominal_pointsource' : 0.55970293,
           'fracflux_nominal_elg' : 0.40354151,
           'fracflux_nominal_bgs' : 0.18546116}

    return par

def gfa_camera_gain_dict():
    # from 20201123 dome screen measurements

    gains = {'GUIDE0' : 3.632,
             'GUIDE2' : 3.731,
             'GUIDE3' : 3.732,
             'GUIDE5' : 3.748,
             'GUIDE7' : 3.839,
             'GUIDE8' : 3.718}

    return gains

def gfa_camera_gain(extname):
    assert(is_valid_extname(extname))

    gains = gfa_camera_gain_dict()

    return gains[extname]

def gfa_camera_readnoise_dict():
    # units are electrons per pixel

    readnoise_electrons = {'GUIDE0': 19.4077,
                           'GUIDE2': 22.0799,
                           'GUIDE3': 18.9685,
                           'GUIDE5': 18.5241,
                           'GUIDE7': 18.6089,
                           'GUIDE8': 21.0796}

    return readnoise_electrons

def gfa_camera_readnoise(extname):
    assert(is_valid_extname(extname))

    r = gfa_camera_readnoise_dict()
    return r[extname]

def valid_image_extname_list():
    return ['GUIDE0', 'GUIDE2', 'GUIDE3', 'GUIDE5', 'GUIDE7', 'GUIDE8']

def valid_extname_list(fz=True):

    par = gfa_misc_params()
    extnames = valid_image_extname_list() + [par['gfa_exp_extname'], par['guider_exp_extname']]

    return extnames

def valid_gfa_number_list():

    return [0, 2, 3, 5, 7, 8]

def gfa_extname_to_gfa_number_dict():
    return dict(zip(valid_image_extname_list(), valid_gfa_number_list()))

def gfa_extname_to_gfa_number(extname):

    # trim whitespace in case the input somehow has any
    extname = extname.replace(' ', '')

    assert(is_valid_extname(extname))

    d = gfa_extname_to_gfa_number_dict()

    return d[extname]

def is_valid_image_extname(extname):
    return (extname in valid_image_extname_list())

def is_valid_extname(extname, fz=True):
    """
    this is different from the above because when fzipped there's
    a dummy zeroth extension
    """

    return (extname in valid_extname_list(fz=fz))

def mask_bit_dict():
    # keys are strings with shorthand name for each bit, values
    # are the corresponding powers of 2

    d = {'HOTDARK' : 0,
         'SATUR' : 1,
         'FLATBAD' : 2}

    return d

def mask_bit_from_bitname(bitname):

    d = mask_bit_dict()

    assert(bitname in d.keys())

    return d[bitname]

def mask_bit_description_dict():

    d = {'HOTDARK' : 'hot pixel based on master dark',
         'SATUR' : 'saturated pixel',
         'FLATBAD' : 'low flatfield value'}

    return d

def mask_bit_description(bitname):
    # bitname is the shorthand name for 

    d = mask_bit_description_dict()

    assert(bitname in d.keys())

    return d[bitname]

def reduced_image_filename_label(flavor):
    par = gfa_misc_params()

    assert(flavor in par['reduced_image_flavors'])

    label = '_' + flavor.lower()

    return label

def reduced_flavor_to_bunit(flavor):
    par = gfa_misc_params()

    assert(flavor in par['reduced_image_flavors'])

    d = {'REDUCED' : 'ADU', 'INVVAR' : '1/ADU^2', 
         'BITMASK' : 'dimensionless', 'DETMAP' : 'sigma'}

    return d[flavor]

def expid_from_filename(fname):
    # this is meant to be run on an input fname that is a single string
    # as opposed to an array/list of strings

    assert(isinstance(fname, str))

    f = os.path.basename(fname)

    _string = f.split('-')[1]
    
    expid = int(_string[0:8])

    return expid

def valid_amps_list():
    return ['E', 'F', 'G', 'H']

def is_valid_amp(amp):
    return amp in valid_amps_list()

def overscan_bdy_coords(amp):
    assert(is_valid_amp(amp))

    # pixel coordinates meant to be used to slice numpy arrays
    if amp  == 'E':
        return {'x_l': 1074, 'x_u': 1124, 'y_l': 0, 'y_u': 516}
        
    elif amp == 'F':
        return {'x_l': 1124, 'x_u': 1174, 'y_l': 0, 'y_u': 516}
    
    elif amp == 'G':
        return {'x_l': 1124, 'x_u': 1174, 'y_l': 516, 'y_u': 1032}

    elif amp == 'H':
        return {'x_l': 1074, 'x_u': 1124, 'y_l': 516, 'y_u': 1032}

def prescan_bdy_coords(amp):
    assert(is_valid_amp(amp))

    if amp == 'E':
        return {'x_l': 0, 'x_u': 50, 'y_l': 0, 'y_u': 516}
        
    elif amp == 'F':
        return {'x_l': 2198, 'x_u': 2248, 'y_l': 0, 'y_u': 516}

    elif amp == 'G':
        return {'x_l': 2198, 'x_u': 2248, 'y_l': 516, 'y_u': 1032}

    elif amp == 'H':
        return {'x_l': 0, 'x_u': 50, 'y_l': 516, 'y_u': 1032}

def amp_bdy_coords(amp):
    assert(is_valid_amp(amp))

    # assumes prescan/overscan have already been stripped out

    if amp == 'E':
        return {'x_l': 0, 'x_u': 1024, 'y_l': 0, 'y_u': 516}
            
    elif amp == 'F':
        return {'x_l': 1024, 'x_u': 2048, 'y_l': 0, 'y_u': 516}

    elif amp == 'G':
        return {'x_l': 1024, 'x_u': 2048, 'y_l': 516, 'y_u': 1032}
            
    elif amp == 'H':
        return {'x_l': 0, 'x_u': 1024, 'y_l': 516, 'y_u': 1032}
