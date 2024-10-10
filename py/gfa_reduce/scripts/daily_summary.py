# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.scripts.daily_summary
================================

This module is designed to be the equivalent of::

    python -u ${GFA_REDUCE}/py/gfa_reduce/scripts/gfa_recent_nights.py
    python ${GFA_REDUCE}/py/gfa_reduce/scripts/run_concat_ccds_nersc.py

and to be the code behind :command:`desi_gfa_reduce_daily_summary`.
"""
import os
import sys
import time
from argparse import ArgumentParser

import numpy as np

from desiutil.iers import freeze_iers
from desiutil.log import get_logger, DEBUG

from .concat_ccds import nights_list, _get_default_basedir, _append_many_nights
from .gfa_single_night import _gfa_single_night


def gfa_recent_nights(add_nights=None, workers=8):
    """Process all recent nights with GFA data that have not been processed yet.

    This is a wrapper on several other functions.

    Parameters
    ----------
    add_nights : :class:`list`, optional
        If set, explicitly add the nights in `add_nights` to the
        list of nights to process.
    workers : :class:`int`, optional
        Number of workers to create, default 8.
    """
    log = get_logger()
    basedir = os.environ['DESI_SPECTRO_DATA']
    out_basedir = os.environ['DEFAULT_REDUX_DIR']

    log.debug('basedir = %s', basedir)
    log.debug('out_basedir = %s', out_basedir)

    nights = nights_list('20210405', '21000101', basedir=basedir)
    processed_nights = nights_list('20210405', '21000101', basedir=out_basedir, empty=True)

    log.debug('Most recent night with raw GFA data: %s', max(nights))
    log.debug('Most recent night with processed GFA data: %s', max(processed_nights))

    if (max(processed_nights) >= max(nights)) and (add_nights is None or len(add_nights) == 0):
        log.info('No new raw GFA data to process.')
        return

    _nights = np.array(nights)
    _nights = _nights[(_nights > max(processed_nights))]
    if add_nights is not None and len(add_nights) > 0:
        _nights = np.concatenate((_nights, np.array(add_nights)))
    nights = np.unique(_nights).tolist()

    nights.sort()

    t0 = time.time()
    num_processed = 0
    for night in nights:
        _n = _gfa_single_night(night=night, numworkers=workers,
                               out_basedir=out_basedir, guider=True, focus=False,
                               indir=basedir)
        num_processed += _n

    dt = time.time() - t0

    log.info('gfa_recent_nights took %.2f seconds', dt)
    log.info('Number of GFA guide images processed: %d', num_processed)
    return


def gfa_daily_summary(workers=8):
    """Produce summary file of all recent nights.

    This is a wrapper on several other functions.

    Parameters
    ----------
    workers : :class:`int`, optional
        Number of workers to create, default 8.
    """
    log = get_logger()
    basedir = _get_default_basedir(acq=False)
    outdir = os.environ['GFA_SUMMRY_FILE_DIR']
    log.debug('basedir = %s', basedir)
    log.debug('outdir = %s', outdir)
    _append_many_nights(night_min='20210405', basedir=basedir,
                        outdir=outdir, workers=workers)
    return


def _options():
    """Parse command-line options.
    """
    parser = ArgumentParser(description="Reduce recent GFA exposures, possibly from multiple nights, and produce a summary file.",
                            prog=os.path.basename(sys.argv[0]))
    parser.add_argument("-n", '--nights', nargs='*', metavar='NIGHT',
                        help='Process one or more NIGHT in addition to the most recent night(s).')
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="Print extra debug information.")
    parser.add_argument('-w', '--workers', type=int, metavar='N', default=8,
                        help='Number of worker instances, default %(default)s.')
    return parser.parse_args()


def main():
    """Entry-point for command-line scripts.

    Returns
    -------
    :class:`int`
        An integer suitable for passing to :func:`sys.exit`.
    """
    options = _options()
    if options.verbose:
        os.environ['DESI_LOGLEVEL'] = 'DEBUG'
    freeze_iers()
    gfa_recent_nights(add_nights=options.nights, workers=options.workers)
    gfa_daily_summary(workers=options.workers)
    return 0
