#!/usr/bin/env python

import argparse
import os
import gfa_reduce.io as io
import glob
import time
import numpy as np
import copy
import gfa_reduce.common as common
import astropy.io.fits as fits
import gfa_reduce.dark_current as dark_current
from gfa_reduce.analysis.sky import adu_to_surface_brightness

def print_imstats_1exp(imstats, fname_in, verbose=False):

    # this is a small amount of unnecessary I/O but whatever
    h = fits.getheader(fname_in, extname='CI')

    exptime = h['EXPTIME']

    if verbose:
        for row in imstats:
            print(row['camera'], fname_in)
            for c in row.colnames:
                print(' '*5, '{:16}'.format(c), ' : ', row[c])

    cols = ['expid', 'camera', 'median', 'max', 'min', 'sigma', 'med-bias',
            'dark_tot_adu', 'sky_ab']

    _imstats = copy.deepcopy(imstats)

    # save horizntal space on printouts
    _imstats['median'] = np.round(_imstats['median']).astype('int')
    _imstats['min'] = np.round(_imstats['min']).astype('int')
    _imstats['max'] = np.round(_imstats['max']).astype('int')
    _imstats['sigma'] = np.round(_imstats['sig_robust']).astype('int')
    _imstats['expid'] = common.expid_from_filename(fname_in)

    median_minus_bias = np.zeros(len(_imstats))
    total_dark_adu = np.zeros(len(_imstats))
    sky_mag_ab = np.zeros(len(_imstats))
    for i in range(len(_imstats)):
        median_minus_bias[i] = np.round(imstats['median'][i] - common.get_median_bias_adu(_imstats['camera'][i])).astype(int)
        ccdtemp = io.get_temperature_celsius(fname_in, _imstats['camera'][i])
        # should really get rid of the hardcoding of 7.5 Celsius below !!!
        total_dark_adu[i] = exptime*common.get_median_dark_current(_imstats['camera'][i])*dark_current.dark_current_rate(ccdtemp)/dark_current.dark_current_rate(7.5)
        sky_mag_ab[i] = adu_to_surface_brightness(median_minus_bias[i]-total_dark_adu[i], exptime,_imstats['camera'][i])

    _imstats['med-bias'] = median_minus_bias
    _imstats['dark_tot_adu'] = total_dark_adu
    _imstats['sky_ab'] = sky_mag_ab

    print(_imstats[cols])
    print('*sigma column is a robust standard deviation measurement')
    print('**all pixel values quoted are in ADU')

def imstats_1exp(fname_in, verbose=False):

    for i in range(5):
        exp = io.load_exposure(fname_in, realtime=True)
        if exp is not None:
            break
        else:
            time.sleep(2.0)

    imstats = io.gather_pixel_stats(exp)

    colnames = imstats.colnames
    colnames.insert(0, colnames.pop())

    imstats = imstats[colnames]
    print_imstats_1exp(imstats, fname_in, verbose=verbose)


if __name__=="__main__":
    # works on simulated data
    #fname_in = '/project/projectdirs/desi/users/ameisner/CI/ci_data_challenge/sims/dci-01402.fits'

    # try a real frame
    #fname_in = '/project/projectdirs/desi/spectro/data/20190330/00002930/ci-00002930.fits.fz'

    flist = glob.glob('/project/projectdirs/desi/spectro/data/20190330/*/ci*.fits.fz')

    for f in flist:
        imstats_1exp(f)

    #fname = '/project/projectdirs/desi/spectro/data/20190403/00003746/ci-00003746.fits.fz'
    imstats_1exp(fname)
