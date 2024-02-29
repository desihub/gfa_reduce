# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.analysis.recalib_astrom
==================================

Run astrometric calibration given a catalog with the centroids and
an initial guess (SKYRA, SKYDEC) of the field of view center.
This is mainly going to be a wrapper for :func:`~gfa_reduce.analysis.asterisms.pattern_match`.
"""
from .asterisms import pattern_match, gaia_cat_for_exp
import numpy as np
import astropy.io.fits as fits
from multiprocessing import Pool
import time
from astropy.table import Table
from desiutil.log import get_logger


def recalib_astrom(cat, fname_raw, mjd=None, h=None, mp=False,
                   arcmin_max=6.0, gfa_targets=None):
    # cat should be the catalog for an entire exposure
    log = get_logger()
    if cat is None:
        return None

    extnames = np.unique(cat['camera'])

    if h is None:
        try:
            h = fits.getheader(fname_raw, extname='GFA')
        except:
            h = fits.getheader(fname_raw, extname='GUIDER')

    if gfa_targets is None:
        gaia = gaia_cat_for_exp(h['SKYRA'], h['SKYDEC'], mjd=mjd)
        # conserve memory
        gaia = gaia[['ra', 'dec']] # only columns needed for pattern matching
    else:
        gaia = Table()
        gaia['ra'] = gfa_targets['TARGET_RA']
        gaia['dec'] = gfa_targets['TARGET_DEC']

    log.info('astrometry search using %.1f arcminute radius', arcmin_max)
    args = []
    for extname in extnames:
        _cat = cat[(cat['camera'] == extname) & cat['valid_astrom_calibrator']]
        if len(_cat) < 2:
            _cat = cat[cat['camera'] == extname]
        args.append((_cat, h['SKYRA'], h['SKYDEC'], extname, gaia, arcmin_max))

    t0 = time.time()
    if not mp:
        result = []
        for tup in args:
            result.append(pattern_match(*tup))
    else:
        log.info('Running astrometric pattern matching for all guide cameras in parallel...')
        nproc = len(args)
        assert(nproc <= 6)
        p = Pool(nproc)
        result = p.starmap(pattern_match, args)

    det_ids_used = set(np.concatenate([r['det_ids_used'] for r in result]))
    used_astrom_calibrator = np.array([(d in det_ids_used) for d in cat['det_id']])
    used_astrom_calibrator = used_astrom_calibrator.astype('uint8')

    cat['used_astrom_calibrator'] = used_astrom_calibrator

    dt = time.time()-t0

    log.info('time taken to astrometrically recalibrate all cameras: %s seconds', dt)

    for r in result:
        if r is not None:
            log.info('%s CONTRAST = %s', r['extname'], r['contrast'])

    return result
