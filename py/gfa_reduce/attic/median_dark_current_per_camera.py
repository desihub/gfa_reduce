import astropy.io.fits as fits
import ci_reduce.common as common
import numpy as np

fname = '/project/projectdirs/desi/users/ameisner/CI/post_install_calibs/CI_master_dark.fits'

extnames = common.valid_image_extname_list()

for extname in extnames:
    im = fits.getdata(fname, extname=extname)
    print(extname, np.median(im))
