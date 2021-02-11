#!/usr/bin/env python

import argparse
import os
import gfa_reduce.io as io
from datetime import datetime
import gfa_reduce.analysis.util as util
import gfa_reduce.common as common
import gfa_reduce.analysis.recalib_astrom as wcs
import gfa_reduce.analysis.dm as dm
import time
import desiutil.iers
desiutil.iers.freeze_iers()

class ProcObj():
    def __init__(self, fname_in, gitrev):
        self.fname_in = fname_in
        self.gitrev = gitrev

def acquire_field(fname_in=None, gfa_targets=None, exp_data=None):
    """Process acquisition image to obtain a global astrometric solution
    Args:
        fname_in    : input acquisition image filename
        gfa_targets : FITS_rec table of relevant astrometric calibrators
        exp_data    : HDUList encapsulating the acquisition image

    Notes:
        exactly one of fname_in and exp_data should be specified
        if gfa_targets is not provided, then Gaia HEALPix chunk files are used

    returns a desimeter fieldmodel object
    """

    assert((fname_in is not None) ^ (exp_data is not None))

    fm = _proc(fname_in=fname_in, no_ps1_xmatch=True, skip_image_outputs=True,
               dont_write_invvar=True, skip_psf_models=True,
               skip_raw_imstats=True,
               dont_write_catalog=True, dont_write_ccds=True,
               return_fieldmodel=True, multiproc=True,
               skip_aper_phot=True, det_sn_thresh=10.0, apply_flatfield=False,
               search_rad_arcmin=1.5, do_sky_mag=False,
               gfa_targets=gfa_targets, exp_data=exp_data,
               minimal_ccds_metadata=True, skip_2d_gaussians=True)

    return fm

def _proc(fname_in=None, outdir=None, careful_sky=False,
          no_cataloging=False,
          no_gaia_xmatch=False, no_ps1_xmatch=False,
          cube_index=None, skip_image_outputs=False,
          realtime=False, no_dark_rescaling=False, 
          dont_write_invvar=False, skip_psf_models=False,
          compress_reduced_image=False, skip_raw_imstats=False,
          skip_astrometry=False, no_pm_pi_corr=False, write_psf_cubes=False,
          write_detmap=False, write_full_detlist=False, max_cbox=31,
          fieldmodel=False, dont_write_catalog=False,
          dont_write_ccds=False, return_fieldmodel=False,
          multiproc=False, skip_aper_phot=False,
          det_sn_thresh=5.0, apply_flatfield=True,
          search_rad_arcmin=6.0, do_sky_mag=True, gfa_targets=None,
          exp_data=None, minimal_ccds_metadata=False,
          skip_2d_gaussians=False, mjdmin=None, mjdmax=None, pmgstars=False):

    print('Starting GFA reduction pipeline at: ' + str(datetime.utcnow()) + 
          ' UTC')

    t0 = time.time()

    try:
        print('Running on host: ' + str(os.environ.get('HOSTNAME')))
    except:
        print('Could not retrieve hostname!')
    
    write_outputs = (outdir is not None)

    # one (but not both) of fname_in and exp_data should be specified
    assert((fname_in is not None) ^ (exp_data is not None))

    if fname_in is not None:
        assert(os.path.exists(fname_in))

    gitrev = io.retrieve_git_rev()

    proc_obj = ProcObj(fname_in, gitrev)
    
    if write_outputs:
        # fail if ANY of expected outputs already exist
        io.check_image_level_outputs_exist(outdir, fname_in, gzip=True,
                                           cube_index=cube_index)

    if (mjdmin is None) or (mjdmax is None):
        mjdrange = None
    else:
        mjdrange = [mjdmin, mjdmax]

    exp = io.load_exposure(fname=fname_in, cube_index=cube_index,
                           realtime=realtime,
                           store_detmap=write_detmap, max_cbox=max_cbox,
                           hdul=exp_data, mjdrange=mjdrange)

    # check for simulated data
    if (exp is None) or util.has_wrong_dimensions(exp):
        # point is to not crash, for sake of real time reductions
        print('EXITING: exposure may be a simulation or contain only focus camera images?!')
        return

    if exp_data is not None:
        fname_in = exp.fname_in # fake file name

    del exp_data

    imstats = io.gather_pixel_stats(exp, skip=skip_raw_imstats)

    # create data quality bitmasks
    # exp.create_all_bitmasks() # revisit this later

    exp_mjd = exp.exp_header['MJD-OBS']
    
    # go from "raw" images to "reduced" images
    exp.calibrate_pixels(do_dark_rescaling=(not no_dark_rescaling),
                         mp=multiproc, do_apply_flatfield=apply_flatfield)

    if do_sky_mag:
        # calculate sky brightness in mag per sq asec
        exp.estimate_all_sky_mags(careful_sky=careful_sky,
                                  flatfield_applied=apply_flatfield)
        exp.estimate_all_sky_sigmas(careful_sky=careful_sky)

    par = common.gfa_misc_params()

    if not no_cataloging:
        catalogs = exp.all_source_catalogs(mp=multiproc,
                                           run_aper_phot=(not skip_aper_phot),
                                           det_sn_thresh=det_sn_thresh,
                                           skip_2dg=skip_2d_gaussians)

        for extname, cat in catalogs.items():
            if cat is not None:
                util.create_det_ids(cat, extname, fname_in,
                                    cube_index=cube_index)

        # reformat the output catalogs into a single merged astropy Table
        catalog = io.combine_per_camera_catalogs(catalogs)

        if not skip_astrometry:
            # run astrometric recalibration
            print('Attempting astrometric recalibration relative to Gaia DR2')

            if gfa_targets is not None:
                util.pm_pi_corr_fiberassign(gfa_targets, exp_mjd)

            astr = wcs.recalib_astrom(catalog, fname_in,
                                      mjd=(None if no_pm_pi_corr else exp_mjd),
                                      h=exp.exp_header, mp=multiproc,
                                      arcmin_max=search_rad_arcmin,
                                      gfa_targets=gfa_targets)
            exp.update_wcs(astr)
            exp.recompute_catalog_radec(catalog)

        if not skip_psf_models:
            exp.compute_psfs(catalog)
            
        if (not no_ps1_xmatch) and (par['ps1_env_var'] in os.environ):
            # probably should look into dec < -30 handling more at some point
            print('Attempting to perform PS1 cross-matching...')
            ps1 = io.get_ps1_matches(catalog)
        else:
            ps1 = None
            
        if (not no_gaia_xmatch) and (par['gaia_env_var'] in os.environ):
            print('Attempting to identify Gaia cross-matches')
            catalog = io.append_gaia_crossmatches(catalog,
                                                  mjd=(None if no_pm_pi_corr else exp_mjd), gfa_targets=gfa_targets)
    else:
        catalog = None
        ps1 = None

    # make this work correctly in the case that --no_cataloging is set
    ccds = io.assemble_ccds_table(imstats, catalog, exp, outdir, proc_obj,
                                  cube_index=cube_index, ps1=ps1,
                                  det_sn_thresh=det_sn_thresh,
                                  sky_mags=do_sky_mag,
                                  minimal=minimal_ccds_metadata,
                                  mjdrange=mjdrange)
    
    if write_outputs:

        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        if not skip_image_outputs:
            print('Attempting to write image-level outputs to directory : ' + 
                  outdir)
            # could add command line arg for turning off gzip compression
            io.write_image_level_outputs(exp, outdir, proc_obj, gzip=True,
                                         cube_index=cube_index,
                                         dont_write_invvar=dont_write_invvar,
                                         compress_reduced_image=compress_reduced_image,
                                         write_detmap=write_detmap, mjdrange=mjdrange)

        if not dont_write_ccds:
            io.write_ccds_table(ccds, outdir, proc_obj, cube_index=cube_index)

        if not no_cataloging:
            if not dont_write_catalog:
                io.write_exposure_source_catalog(catalog, outdir, proc_obj, exp,
                                                 cube_index=cube_index)
            io.write_ps1_matches(ps1, outdir, fname_in,
                                 cube_index=cube_index)
            if not skip_psf_models:
                io.write_psfs(exp, outdir, fname_in,
                                    cube_index=cube_index)
                if write_psf_cubes:
                    io.write_psfs(exp, outdir, fname_in,
                                        cube_index=cube_index, cubes=True)

            if write_full_detlist:
                io.write_full_detlists(exp, outdir, fname_in, cube_index=cube_index)

        if pmgstars:
            io.write_pmgstars(exp, outdir, fname_in, cube_index)

    # desimeter fieldmodel if applicable
    if (not skip_astrometry) and (fieldmodel or return_fieldmodel):
        # should probably log the timing of this step
        check_header_cards = (cube_index is not None) or (fname_in.find('-0000.fits.fz') != -1)
        fm = dm.fit_dm_fieldmodel(exp.exp_header, ccds, catalog,
                                  check_header_cards=check_header_cards)
        if not return_fieldmodel:
            io.write_dm_fieldmodel(fm, outdir, fname_in, cube_index=cube_index)
    
    print('Successfully finished reducing ' + fname_in)

    dt = time.time() - t0
    print('GFA reduction pipeline took ' + '{:.2f}'.format(dt) + ' seconds')
    print('GFA reduction pipeline completed at: ' + str(datetime.utcnow()) + 
          ' UTC')

    # for field acquisition mode
    if return_fieldmodel:
        return fm

if __name__ == "__main__":
    descr = 'run the gfa_reduce pipeline on a GFA exposure'
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('fname_in', type=str, nargs=1)

    parser.add_argument('--outdir', default=None, type=str,
                        help='directory to write outputs in')

    parser.add_argument('--careful_sky', default=False, action='store_true',
        help='use image segmentation when deriving sky quantities')

    parser.add_argument('--no_cataloging', default=False, action='store_true', 
                        help='reduce image without cataloging sources')

    parser.add_argument('--no_gaia_xmatch', default=False, action='store_true',
                        help='skip Gaia cross-match')

    parser.add_argument('--no_ps1_xmatch', default=False, action='store_true',
                        help='skip PS1 cross-match')

    parser.add_argument('--cube_index', default=None, type=int,
                        help='guide cube index')

    parser.add_argument('--skip_image_outputs', default=False,
                        action='store_true',
                        help='skip writing of full-frame image outputs')

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

    parser.add_argument('--max_cbox', default=31, type=int,
                        help="maximum centroiding box size (pixels)")

    parser.add_argument('--fieldmodel', default=False, action='store_true',
                        help="fit and write desimeter field model")

    parser.add_argument('--dont_write_catalog', default=False, action='store_true',
                        help="don't write source catalog")

    parser.add_argument('--dont_write_ccds', default=False, action='store_true',
                        help="don't write CCDs table")

    parser.add_argument('--multiproc', default=False, action='store_true',
                        help="use multiprocessing to decrease wall time")

    parser.add_argument('--skip_aper_phot', default=False, action='store_true',
                        help="don't perform aperture photometry")

    parser.add_argument('--det_sn_thresh', default=5.0, type=float,
                        help="source detection significance threshold")

    parser.add_argument('--skip_flatfield', default=False, action='store_true',
                        help="skip flatfielding during pixel-level detrending")

    parser.add_argument('--search_rad_arcmin', default=6.0, type=float,
                        help="astrometric pattern match search radius (arcmin)")

    parser.add_argument('--skip_sky_mag', default=False, action='store_true',
                        help="skip sky surface brightness estimation")

    parser.add_argument('--mjdmin', default=None, type=float,
                        help='minimum MJD for guide cube coaddition')

    parser.add_argument('--mjdmax', default=None, type=float,
                        help='maximum MJD for guide cube coaddition')

    parser.add_argument('--pmgstars', default=False, action='store_true',
                        help='do PMGSTARS forced photometry')

    args = parser.parse_args()

    fname_in = args.fname_in[0]

    _proc(fname_in=fname_in, outdir=args.outdir, careful_sky=args.careful_sky,
          no_cataloging=args.no_cataloging, no_gaia_xmatch=args.no_gaia_xmatch,
          no_ps1_xmatch=args.no_ps1_xmatch, cube_index=args.cube_index,
          skip_image_outputs=args.skip_image_outputs, realtime=args.realtime,
          no_dark_rescaling=args.no_dark_rescaling, 
          dont_write_invvar=args.dont_write_invvar,
          skip_psf_models=args.skip_psf_models,
          compress_reduced_image=args.compress_reduced_image,
          skip_raw_imstats=args.skip_raw_imstats,
          skip_astrometry=args.skip_astrometry,
          no_pm_pi_corr=args.no_pm_pi_corr,
          write_psf_cubes=args.write_psf_cubes,
          write_detmap=args.write_detmap,
          write_full_detlist=args.write_full_detlist,
          max_cbox=args.max_cbox, fieldmodel=args.fieldmodel,
          dont_write_catalog=args.dont_write_catalog,
          dont_write_ccds=args.dont_write_ccds,
          multiproc=args.multiproc, skip_aper_phot=args.skip_aper_phot,
          det_sn_thresh=args.det_sn_thresh,
          apply_flatfield=(not args.skip_flatfield),
          search_rad_arcmin=args.search_rad_arcmin,
          do_sky_mag=(not args.skip_sky_mag), mjdmin=args.mjdmin,
          mjdmax=args.mjdmax, pmgstars=args.pmgstars)
