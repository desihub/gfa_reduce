import astropy.io.fits as fits
import numpy as np

fname_steve = '/project/projectdirs/desi/users/ameisner/CI/post_install_calibs/CI_master_dark-no_pdu_image.fits'
fname = '/project/projectdirs/desi/users/ameisner/CI/post_install_calibs/CI_master_dark.fits'

hdul_steve = fits.open(fname_steve)
assert(hdul_steve[0].data is None)


assert(len(hdul_steve) == 6)

for extname in ['CIE', 'CIN', 'CIC', 'CIS', 'CIW']:
    old = fits.getdata(fname, extname=extname)
    new = fits.getdata(fname_steve, extname=extname)
    ndiff = np.sum(old != new)

    print(len(np.ravel(old)))
    print(len(np.ravel(new)))
    print(np.median(old), np.median(new))
    assert(ndiff == 0)
