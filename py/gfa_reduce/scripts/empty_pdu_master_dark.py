import astropy.io.fits as fits
import os

fname_master = '/project/projectdirs/desi/users/ameisner/CI/post_install_calibs/CI_master_dark.fits'


hdul_raw = fits.open('/project/projectdirs/desi/spectro/data/20190405/00004111/ci-00004111.fits.fz')

dummy_hdu = hdul_raw[0]
dummy_hdu.header['BITPIX'] = (-32, 'array data type')

hdul_master = fits.open(fname_master)

hdul_out = [dummy_hdu]
for hdu in hdul_master:
    hdul_out.append(fits.ImageHDU(hdu.data, header=hdu.header))

outname = '/project/projectdirs/desi/users/ameisner/CI/post_install_calibs/CI_master_dark-no_pdu_image.fits'
assert(not os.path.exists(outname))

hdul_out = fits.HDUList(hdul_out)
hdul_out.writeto(outname)
