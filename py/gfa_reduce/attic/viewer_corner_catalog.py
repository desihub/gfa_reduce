import ci_reduce.common as common
import ci_reduce.analysis.util as util
import astropy.io.fits as fits
from astropy import wcs
import argparse
import ci_reduce.ci_wcs as ci_wcs
from astropy.table import Table, vstack
import os

def corner_catalog_1ext(telra, teldec, extname):
    # LL -> UL -> UR -> LR
    xpix, ypix = util.ci_corner_pixel_coords()

    wcs = ci_wcs.nominal_tan_wcs(telra, teldec, extname)

    ra, dec = wcs.all_pix2world(xpix, ypix, 0)

    name = [(extname + '_' + corner) for corner in ['LL', 'UL', 'UR', 'LR']]

    tab = Table([ra, dec, name], names=('ra', 'dec', 'name'))
    
    return tab

def create_catalog_1exp(telra, teldec, outname=None):

    tabs = []
    for extname in common.valid_image_extname_list():
        tab = corner_catalog_1ext(telra, teldec, extname)
        tabs.append(tab)

    tab = vstack(tabs)

    if outname is not None:
        assert(not os.path.exists(outname))
        tab.write(outname, format='fits')

    return tab

if __name__ == "__main__":
    descr = 'run full ci_reduce pipeline on a CI exposure'
    parser = argparse.ArgumentParser(description=descr)

# option to either input a desi tileid or else input both telra and teldec
