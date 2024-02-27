import ci_reduce.common as common
import astropy.io.fits as fits
import glob
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import scoreatpercentile

# note that is assumed to be constant !!!
dark_exptime = 300.0 # 300 seconds was time used for all long darks Klaus ran

def get_flist():
    dirname = '/project/projectdirs/desi/spectro/data/20190330/'

    flist = glob.glob(dirname + '/*/ci*.fits.fz')

    exptime = []
    flavor = []
    for f in flist:
        h = fits.getheader(f, extname='CIC')

        flavor.append(h['FLAVOR'])
        exptime.append(h['EXPTIME'])

    exptime = np.array(exptime)
    flavor = np.array(flavor)

    good = (exptime == dark_exptime) & (flavor == 'DARK')

    flist = np.array(flist)
    flist = flist[good]

    return flist

def read_dark_frames(ci_extname):
    flist = get_flist()
    imgs = [fits.getdata(f, extname=ci_extname) for f in flist]

    return imgs

def master_dark_1camera(ci_extname):
    # use common to check that CI extention name is valid

    dark_frames = read_dark_frames(ci_extname)
    dark_frames = np.asarray(dark_frames)
    dark_med = np.median(dark_frames, axis=0)

    # we do want the 20190330 version of the bias, since that has the
    # same mislabeling of extensions as the darks

    bias = fits.getdata('/project/projectdirs/desi/users/ameisner/CI/post_install_calibs/CI_master_bias-20190330.fits', extname=ci_extname)

    dark_med -= bias

    return dark_med

def master_dark_header_cards(hdu, ci_extname):
    h = hdu.header
    h['FLAVOR'] = 'DARK'
    h['EXTNAME'] = ci_extname
    h['BUNIT'] = 'ADU'
    h['EXPTIME'] = (1, 'seconds')
    # this is just a placeholder guess; CCDTEMP was not available
    # when the 300 s darks were taken, and as of April 3, 2019 various
    # problems seems to remain with CCDTEMP
    # this nominal temperature should definitely be re-evaluated
    # in the future
    # Ashley has informed me that the CI cameras are operating a mode
    # where they try to cool themselves to 7.5 Celsius, so that's where 7.5 
    # comes from

    h['CCDTEMP'] = (7.5, '[deg] Camera temperature (Celsius)')

    return h

def write_master_dark(outname=None):

    par = common.ci_misc_params()

    if outname is None:
        outname = os.path.join('/project/projectdirs/desi/users/ameisner/CI/post_install_calibs/CI_master_dark-20190330.fits')

    assert(not os.path.exists(outname))

    ci_extnames = common.valid_image_extname_list()

    hdus = []
    for ci_extname in ci_extnames:
        print('Working on master dark for: ' + ci_extname)
        dark_image = master_dark_1camera(ci_extname)
        dark_image = dark_image.astype('float32')

        # convert to counts per second !!!
        dark_image = dark_image/dark_exptime

        dark_image = dark_image.astype('float32')

        if len(hdus) == 0:
            hdu = fits.PrimaryHDU(dark_image)
        else:
            hdu = fits.ImageHDU(dark_image)
        hdu.header = master_dark_header_cards(hdu, ci_extname)

        hdus.append(hdu)

    hdul = fits.HDUList(hdus)

    hdul.writeto(outname)
