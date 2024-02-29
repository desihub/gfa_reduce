import ci_reduce.common as common
import ci_reduce.imred.load_calibs as load_calibs
import astropy.io.fits as fits
import glob
import numpy as np
import os
from scipy import ndimage

# in the future it would be good to have access to some relatively long dark 
# frames to investigate what if any additional masking those might suggest

def flatfield_jump_mask(flat):
    medfilt = ndimage.median_filter(flat, 3)
    jump = (np.abs(flat-medfilt) > 0.1).astype('int')

    return jump

def bad_pixels_master_flat(ci_extname):
    flat = load_calibs.read_flat_image(ci_extname)

    # bit 2^0 is pixels that are truly bad based on master flat
    mask = (flat < 0.5).astype('int')

    # bit 2^1 will be questionable pixels based on master flat
    # questionable means that there's a 1-pixel "jump" in the
    # flat field relative to local neighborhood of surrounding pixels
    mask += 2*flatfield_jump_mask(flat)

    return mask

def static_mask_header_cards(hdu, ci_extname):
    h = hdu.header

    h['FLAVOR'] = 'MASK'
    h['EXTNAME'] = ci_extname

    h['MASKB0'] = ('low flat field value', 'Mask bit 2**0=1 meaning')
    h['MASKB1'] = ('single-pixel flat-field jump', 'Mask bit 2**1=2 meaning')

    return h

def write_bad_pixel_mask():
    # bad pixel mask is currently based only on flat field, and flat field
    # can only be constructed from forDK.tar.gz data for CIN, so 
    # for now all extensions will be the same

    par = common.ci_misc_params()
    outname = par['static_mask_filename']
    outname = os.path.join(os.environ[par['etc_env_var']], outname)

    assert(not os.path.exists(outname))

    mask = bad_pixels_master_flat('CIN')

    hdus = []
    for ci_extnum in range(par['n_cameras']):
        ci_extname = common.ci_extnum_to_extname(ci_extnum, fz=False)
        print('Assembling HDU for: ' + ci_extname)
        if len(hdus) == 0:
            hdu = fits.PrimaryHDU(mask)
        else:
            hdu = fits.ImageHDU(mask)
        hdu.header = static_mask_header_cards(hdu, ci_extname)

        hdus.append(hdu)

    hdul = fits.HDUList(hdus)

    hdul.writeto(outname)
