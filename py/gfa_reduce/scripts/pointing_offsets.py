#!/usr/bin/env python

import argparse
import gfa_reduce.gfa_red as gfa_red
import astropy.io.fits as fits
import numpy as np
import os

basedir = '/exposures/desi'

def pointing_offsets(fm, fname_in):
    if fm is None:
        print('did not obtain a good FieldModel???')

    h = fits.getheader(fname_in, extname='GFA')
    target_ra = h['TARGTRA']
    target_dec = h['TARGTDEC']

    print('cos(Dec) factor = ', np.cos(target_dec/(180.0/np.pi)))

    dra_true_asec = 3600.0*(target_ra - fm.ra)*np.cos(target_dec/(180.0/np.pi))
    ddec_true_asec = 3600.0*(target_dec - fm.dec)
    print('RA POINTING OFFSET : ', '{:.2f}'.format(dra_true_asec), ' asec')
    print('DEC POINTING OFFSET : ', '{:.2f}'.format(ddec_true_asec), ' asec')

def _get_raw_filename(expid, night):
    # assumes file name of the form gfa-????????.fits.fz
    # as opposed to one alternative which I guess would be guide-????????-0000.fits.fz
    # acquisition images

    fname = 'gfa-' + str(expid).zfill(8) + '.fits.fz'
    nightdir = os.path.join(basedir, str(night).zfill(8))
    dir = os.path.join(nightdir, str(expid).zfill(8))
    fname = os.path.join(dir, fname)

    return fname

if __name__ == "__main__":

    # eventually should add a parameter for what angular
    # search radius to use, think default is 1.5 arcmin

    descr = 'compute GFA pointing offsets'

    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('expid', type=int, nargs=1)

    parser.add_argument('night', type=str, nargs=1)

    args = parser.parse_args()

    fname_in = _get_raw_filename(args.expid[0], args.night[0])

    print(fname_in)
    if not os.path.exists(fname_in):
        print('raw gfa-????????.fits.fz file not found??')
        assert(False)

    fm = gfa_red.acquire_field(fname_in)

    pointing_offsets(fm, fname_in)
