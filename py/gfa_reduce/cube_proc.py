# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.cube_proc
====================

Run gfa_reduce pipeline on multiple GFA guide cube frames.

This module is more of a script and may be moved to the attic.
"""
import argparse
from gfa_reduce.gfa_red import _proc
import astropy.io.fits as fits
import os
import numpy as np

if __name__ == "__main__":
    descr = 'run gfa_reduce pipeline on multiple GFA guide cube frames'

    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('fname_in', type=str, nargs=1)

    parser.add_argument('--outdir', default=None, type=str,
                        help='directory to write outputs in')

    parser.add_argument('--indstart', default=0, type=int,
                        help='starting frame index (zero indexed)')

    parser.add_argument('--nproc', default=None, type=int,
                        help='number of frames to process')

    parser.add_argument('--careful_sky', default=False, action='store_true',
        help='use image segmentation when deriving sky quantities')

    parser.add_argument('--no_cataloging', default=False, action='store_true',
                        help='reduce image without cataloging sources')

    parser.add_argument('--no_gaia_xmatch', default=False, action='store_true',
                        help='skip Gaia cross-match')

    parser.add_argument('--no_ps1_xmatch', default=False, action='store_true',
                        help='skip PS1 cross-match')

    parser.add_argument('--write_image_outputs', default=False,
                        action='store_true',
                        help='write image-level reduction outputs')

    parser.add_argument('--realtime', default=False,
                        action='store_true',
                        help='avoid crashing on partially written raw files')

    parser.add_argument('--no_dark_rescaling', default=False,
                        action='store_true',
                        help='skip empirical rescaling of dark current')

    parser.add_argument('--dont_write_invvar', default=False,
                        action='store_true',
                        help="don't write out invvar maps")

    parser.add_argument('--skip_psf_models', default=False,
                        action='store_true',
                        help="skip generating per-camera PSF models")

    parser.add_argument('--compress_reduced_image', default=False,
                        action='store_true',
                        help="compress reduced image output file")

    parser.add_argument('--skip_raw_imstats', default=False,
                        action='store_true',
                        help="skip computing of raw image pixel statistics")

    parser.add_argument('--skip_astrometry', default=False,
                        action='store_true',
                        help='skip astrometric recalibration')

    parser.add_argument('--no_pm_pi_corr', default=False, action='store_true',
        help="do not correct Gaia positions for proper motion or parallax")

    parser.add_argument('--write_psf_cubes', default=False, action='store_true',
        help="write image cubes of sources used to build PSF models")

    parser.add_argument('--write_detmap', default=False, action='store_true',
                        help="write detection map")

    parser.add_argument('--write_full_detlist', default=False, action='store_true',
                        help="write out the initial, full list of detections")

    args = parser.parse_args()

    fname_in = args.fname_in[0]

    assert(os.path.exists(fname_in))

    h = fits.getheader(fname_in, extname='GUIDER')

    nframes = h['FRAMES']
    assert(args.indstart < nframes)

    nproc = args.nproc if args.nproc is not None else nframes

    indend = min(nframes, args.indstart + nproc)

    skip_image_outputs = not args.write_image_outputs
    for i in np.arange(args.indstart, indend):
        # note the hardcoding of certain arguments, especially
        # skipping write-out of image-level outputs (saves disk space...)
        print('WORKING ON FRAME ' + str(i+1) + ' OF ' + str(nframes))
        _proc(fname_in, outdir=args.outdir,
              careful_sky=args.careful_sky, no_cataloging=args.no_cataloging,
              no_gaia_xmatch=args.no_gaia_xmatch,
              no_ps1_xmatch=args.no_ps1_xmatch,
              cube_index=i, skip_image_outputs=skip_image_outputs,
              realtime=args.realtime, no_dark_rescaling=args.no_dark_rescaling,
              dont_write_invvar=args.dont_write_invvar,
              skip_psf_models=args.skip_psf_models,
              compress_reduced_image=args.compress_reduced_image,
              skip_raw_imstats=args.skip_raw_imstats,
              skip_astrometry=args.skip_astrometry,
              no_pm_pi_corr=args.no_pm_pi_corr,
              write_psf_cubes=args.write_psf_cubes,
              write_detmap=args.write_detmap,
              write_full_detlist=args.write_full_detlist)
