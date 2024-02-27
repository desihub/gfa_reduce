import astropy.io.fits as fits
import numpy as np
import os

fname_in = '/project/projectdirs/desi/users/ameisner/CI/deprecated_fordk_calibs/CI_master_flat.fits.gz'

hdul = fits.open(fname_in)

for hdu in hdul:
    image = hdu.data
    image.fill(1.0)

    hdu.data = image

outname = '/project/projectdirs/desi/users/ameisner/CI/ci_reduce_etc/fake_uniform_flat.fits'

assert(not os.path.exists(outname))

hdul.writeto(outname)
