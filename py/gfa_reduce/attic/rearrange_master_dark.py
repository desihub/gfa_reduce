import astropy.io.fits as fits
import os

"""
The CI_master_dark-20190330.fits file was generated using data where
the 'CIE' extension name referred to the sky west camera, the 'CIW' extension 
name referred to the sky east camera, the 'CIN' extension name referred to the 
central camera, and the 'CIC' extension name referred to the sky north camera.

The point of this code is to pair the dark images with extension names
such that the extension name truly does indicate the sky location in all cases
(since the ICS was updated after the 300 s darks were taken to enforce
this matching)
"""

fname_orig = '/project/projectdirs/desi/users/ameisner/CI/post_install_calibs/CI_master_dark-20190330.fits'

# key = old, value = new
mapping = {'CIE' : 'CIW', 
           'CIN' : 'CIC', 
           'CIC' : 'CIN',
           'CIS' : 'CIS', 
           'CIW' : 'CIE'}

hdulist = []
for k,v in mapping.items():
    im, h = fits.getdata(fname_orig, header=True, extname=k)

    h['EXTNAME'] = v

    if len(hdulist) == 0:
        hdu = fits.PrimaryHDU(im, header=h)
    else:
        hdu = fits.ImageHDU(im, header=h)
    hdulist.append(hdu)

outname = '/project/projectdirs/desi/users/ameisner/CI/post_install_calibs/CI_master_dark.fits'
assert(not os.path.exists(outname))
hdulist = fits.HDUList(hdulist)
hdulist.writeto(outname)
