from gfa_reduce.image import GFA_image
from gfa_reduce.exposure import GFA_exposure
import gfa_reduce.common as common
import gfa_reduce.xmatch.gaia as gaia
import astropy.io.fits as fits
from astropy.table import Table, vstack, hstack
import os
import gfa_reduce.analysis.basic_image_stats as bis
import gfa_reduce.analysis.basic_catalog_stats as bcs
import gfa_reduce.analysis.util as util
import numpy as np
import time
from gfa_reduce.gfa_wcs import ccd_center_radec
import json
import astropy

# in the context of this file, "image" and "exposure" generally refer to 
# GFA_image and GFA_exposure objects

def loading_image_extension_message(extname):
    print('Attempting to load image extension : ' + extname)

def load_image_from_hdu(hdu, verbose=True, cube_index=None, store_detmap=False,
                        exp_header=None, coadd_index_range=None):
    loading_image_extension_message(hdu.header['EXTNAME'])

    header = hdu.header

    # hack for PlateMaker acquisition image file format
    if 'SKYRA' not in header:
        header['SKYRA'] = exp_header['SKYRA']
        header['SKYDEC'] = exp_header['SKYDEC']
        header['EXPTIME'] = exp_header['EXPTIME']
        header['REQTIME'] = exp_header['REQTIME']
        header['MJD-OBS'] = exp_header['MJD-OBS']

    if 'EXPID' not in header:
        header['EXPID'] = exp_header['EXPID']

    # found some cases on 20201128 and 20201129 where SKYRA, SKYDEC
    # are correct in the exposure-level header but incorrect in the
    # per-image headers, for example EXPID = 65021
    if ('SKYRA' in exp_header) and ('SKYDEC' in exp_header):
        header['SKYRA'] = exp_header['SKYRA'] if exp_header['SKYRA'] is not None else exp_header['REQRA']
        header['SKYDEC'] = exp_header['SKYDEC'] if exp_header['SKYDEC'] is not None else exp_header['REQDEC']

    return GFA_image(hdu.data, header, cube_index=cube_index,
                     store_detmap=store_detmap,
                     coadd_index_range=coadd_index_range)

def load_image_from_filename(fname, extname):
    assert(os.path.exists(fname))

    loading_image_extension_message(extname)
    assert(common.is_valid_image_extname(extname))

    data, header = fits.getdata(fname, extname=extname, header=True)
    return GFA_image(data, header)

def _atomic_write(data, outname):
    # data should be either an astropy Table or an hdulist

    if outname[-3:] == '.gz':
        outname_tmp = outname.replace('.gz', '.tmp.gz')
    else:
        outname_tmp = outname + '.tmp'
    
    if isinstance(data, Table):
        data.write(outname_tmp, format='fits')
    else:
        data.writeto(outname_tmp)

    os.rename(outname_tmp, outname)
    
def realtime_raw_read(fname, delay=2.0, max_attempts=5):
    """
    attempt to avoid getting burned by partially written files when
    trying to analyze data in real time

    delay is in seconds
    """

    # something has gone badly wrong if the filename doesn't even exist
    # that's not the scenario I'm trying to address here
    assert(os.path.exists(fname))

    hdul = None
    for i in range(max_attempts):
        try:
            hdul = fits.open(fname, lazy_load_hdus=False)
            hdul.verify(option='exception')
            for hdu in hdul:
                _, __ = hdu.data, hdu.header
                ___ = hdu.data.shape
        except:
            print('encountered problem reading ' + fname)
            time.sleep(delay)
        if hdul is not None:
            break

    # die if unable to read file after max_attempts attempts
    assert(hdul is not None)

    return hdul

def load_exposure(fname=None, verbose=True, realtime=False, cube_index=None,
                  store_detmap=False, max_cbox=31, hdul=None, mjdrange=None):

    # mjdrange should be 2 element list : [mjdmin, mjdmax]

    # exactly one of fname, hdul should be specified
    assert((fname is not None) ^ (hdul is not None))

    from_file = fname is not None

    if fname is not None:
        assert(os.path.exists(fname))

        print('Attempting to load exposure : ' + fname)

        if not realtime:
            hdul = fits.open(fname)
        else:
            hdul = realtime_raw_read(fname)
    else:
        # create a fake file name
        fname = 'gfa-' + str(hdul[0].header['EXPID']) + '.3.fits'

    exp_header = None

    par = common.gfa_misc_params()

    is_image_hdu = np.zeros(len(hdul), dtype=bool)
    for i, hdu in enumerate(hdul):
        # real data has another dummy extension added with no EXTNAME
        keywords = [c[0] for c in hdu.header.cards]
        if not ('EXTNAME' in keywords):
            continue
        if hdu.header['EXTNAME'] not in common.valid_extname_list():
            continue
        if (hdu.header['EXTNAME']).strip() in [par['gfa_exp_extname'], par['guider_exp_extname']]:
            exp_header = hdu.header
            continue
        is_image_hdu[i] = True

    # hacks for PlateMaker acquisition image file format
    if exp_header is None:
        exp_header = hdul['PRIMARY'].header
        
    if np.sum(is_image_hdu) == 0:
        print('exposure may contain only focus cameras?')
        return None
    
    w_im = np.where(is_image_hdu)[0]

    is_cube = (len(hdul[w_im[0]].data.shape) == 3)

    assert((is_cube and (cube_index is None)) == False)
    assert(((not is_cube) and (cube_index is not None)) == False)

    bintables = None
    if is_cube:
        bintables = {}
        coadd_ind_ranges = []
        for ind in w_im:
            extname_im = hdul[ind].header['EXTNAME'].strip()
            extname_tab = extname_im + 'T'
            # this will crash if the binary table extension is missing...
            bintables[extname_im] = hdul[extname_tab].data
            coadd_ind_ranges.append(util.coadd_cube_index_range(bintables[extname_im], cube_index, mjdrange))
    else:
        coadd_ind_ranges = [None]*len(w_im)
            
    try:
        imlist = [load_image_from_hdu(hdul[ind], verbose=verbose, cube_index=cube_index, store_detmap=store_detmap, exp_header=exp_header, coadd_index_range=coadd_ind_ranges[i]) for i, ind in enumerate(w_im)]
    except Exception as e:
        print('failed to load exposure at image list creation stage')
        print(e)
        return None

    exp = GFA_exposure(imlist, exp_header=exp_header, bintables=bintables,
                       max_cbox=max_cbox)

    exp.set_bintable_rows()

    # fake file name
    if not from_file:
        exp.assign_input_filename(fname)

    if cube_index != None:
        util._patch_guider_mjd_obs(exp)

    print('Successfully loaded exposure : ' + fname)
    print('Exposure has ' + str(exp.num_images_populated()) + 
          ' image extensions populated')
    print('Populated image extension names are : ' + 
          str(exp.populated_extnames()))

    return exp

def reduced_image_fname(outdir, fname_in, flavor, gzip=True,
                        cube_index=None, outdir_not_needed=False):

    if not outdir_not_needed:
        assert(os.path.exists(outdir))

    outname = os.path.join(outdir, os.path.basename(fname_in))

    # get rid of any ".fz" or ".gz" present in input filename
    outname = outname.replace('.fz', '')
    outname = outname.replace('.gz', '')

    assert(outname[-5:] == '.fits')

    outname = outname.replace('.fits', 
        common.reduced_image_filename_label(flavor) + '.fits')

    if gzip:
        outname += '.gz'

    if cube_index is not None:
        outname = outname.replace('.fits',
                                  '-' + str(cube_index).zfill(5) + '.fits')

    assert(not os.path.exists(outname))

    return outname

def check_image_level_outputs_exist(outdir, fname_in, gzip=True,
                                    cube_index=None):
    par = common.gfa_misc_params()

    for flavor in par['reduced_image_flavors']:
        _ = reduced_image_fname(outdir, fname_in, flavor, gzip=gzip,
                                cube_index=cube_index, outdir_not_needed=True)

def retrieve_git_rev(fname=None):

    if fname is None:
        fname = __file__

    code_dir = os.path.dirname(os.path.realpath(fname))
    cwd = os.getcwd()
    do_chdir = (cwd[0:len(code_dir)] != code_dir)
    if do_chdir:
        os.chdir(code_dir)
    gitrev = os.popen("git rev-parse --short HEAD").read().replace('\n','')
    if do_chdir:
        os.chdir(cwd)
    print('"git rev" version info:', gitrev)

    return gitrev

def write_image_level_outputs(exp, outdir, proc_obj, gzip=True,
                              cube_index=None, dont_write_invvar=False,
                              compress_reduced_image=False,
                              write_detmap=False, mjdrange=None):
    # exp is a GFA_exposure object
    # outdir is the output directory (string)

    par = common.gfa_misc_params()

    flavors_list = par['reduced_image_flavors']

    if not write_detmap:
        flavors_list.remove('DETMAP')

    if dont_write_invvar:
        flavors_list.remove('INVVAR')

    for flavor in par['reduced_image_flavors']:
        _gzip = (gzip if (flavor != 'REDUCED') else compress_reduced_image)
        outname = reduced_image_fname(outdir, proc_obj.fname_in, flavor,
                                      gzip=_gzip, cube_index=cube_index)

        hdulist = exp.to_hdulist(flavor=flavor)

        for hdu in hdulist:
            hdu.header['GITREV'] = proc_obj.gitrev
            if mjdrange is not None:
                hdu.header['REQMJDLO'] = mjdrange[0]
                hdu.header['REQMJDHI'] = mjdrange[1]

        print('Attempting to write ' + flavor + ' image output to ' + 
              outname)

        _atomic_write(hdulist, outname)
        
        print('Successfully wrote ' + flavor + ' image output to ' + 
              outname)

def strip_none_columns(table):
    # can't write an astropy table to FITS if it has columns with None
    # values

    for c in table.colnames:
        if table[c].dtype.str == '|O':
            table.remove_column(c)

    return table

def combine_per_camera_catalogs(catalogs):
    # catalogs is the output of GFA_exposure's all_source_catalogs() method
    # which is a dictionary of astropy QTable's, with the keys
    # being the GFA camera extension names

    # want to add a column to each table giving the GFA camera name, then
    # append the all into one per-exposure table

    assert(type(catalogs).__name__ == 'dict')

    composite_list = []
    for extname, tab in catalogs.items():
        if tab is not None:
            tab['camera'] = extname
            tab['petal_loc'] = np.array([common.gfa_extname_to_gfa_number(extname) for extname in tab['camera']], dtype='uint8')
            composite_list.append(tab)

    # handle case of no sources in any image
    if len(composite_list) == 0:
        return None
    
    composite = vstack(composite_list)
    composite = strip_none_columns(composite)

    composite['extname'] = composite['camera']
    return composite

def write_exposure_source_catalog(catalog, outdir, proc_obj, exp, 
                                  cube_index=None):

    # handle case of exposure with no retained sources
    if catalog is None:
        return
    
    assert(os.path.exists(outdir))

    outname = os.path.join(outdir, os.path.basename(proc_obj.fname_in))

    # get rid of any ".fz" or ".gz" present in input filename
    outname = outname.replace('.fz', '')
    outname = outname.replace('.gz', '')

    assert(outname[-5:] == '.fits')

    outname = outname.replace('.fits', '_catalog.fits')

    if cube_index is not None:
        outname = outname.replace('.fits',
                                  '-' + str(cube_index).zfill(5) + '.fits')
    
    assert(not os.path.exists(outname))

    catalog['fname_in'] = proc_obj.fname_in
    catalog['gitrev'] = proc_obj.gitrev
    expid = util.expid_from_raw_filename(proc_obj.fname_in)
    catalog['expid'] = expid
    catalog['cube_index'] = np.nan if cube_index is None else float(cube_index)

    mjd = np.zeros(len(catalog), dtype=float)
    for camera in np.unique(catalog['camera']):
        mjd[catalog['camera'] == camera] = exp.images[camera].try_retrieve_meta_keyword('MJD-OBS', placeholder=0.0)

    catalog['mjd'] = mjd

    catalog['valid_astrom_calibrator'] = catalog['valid_astrom_calibrator'].astype('uint8')
    
    hdul = []
    hdul.append(fits.PrimaryHDU())
    hdul.append(fits.BinTableHDU(data=catalog, name='CATALOG'))

    for image in exp.images.values():
        if image is None:
            continue
        hdul.append(fits.ImageHDU(data=None, header=image.header, name=image.extname))

    # dummy extension with copy of exposure-level raw data header
    hdul.append(fits.ImageHDU(data=None, header=exp.exp_header, name=exp.exp_header['EXTNAME']))
    
    hdul = fits.HDUList(hdul)
        
    print('Attempting to write source catalog to ' + outname)

    _atomic_write(hdul, outname)

# since this function doesn't do the writing of the PS1 cross-matches
# it may belong somewhere else rather than in 'io'
def get_ps1_matches(catalog):
    # handle case of exposure with no retained sources
    if catalog is None:
        return None

    # intentionally not sending any MJD info when doing ps1 cross-match
    ps1 = gaia.gaia_xmatch(catalog['ra'], catalog['dec'], ps1=True)

    if ps1 is None:
        print('No PS1 matches available -- Dec may be too low??')
        return None
    
    ps1.rename_column('ra', 'ra_ps1')
    ps1.rename_column('dec', 'dec_ps1')

    ps1_matches = hstack([catalog, ps1])

    return ps1_matches
    
def write_ps1_matches(ps1_matches, outdir, fname_in, cube_index=None):

    if ps1_matches is None:
        print('no PS1 cross-matches to write...')
        return
    
    assert(os.path.exists(outdir))

    outname = os.path.join(outdir, os.path.basename(fname_in))

    # get rid of any ".fz" or ".gz" present in input filename
    outname = outname.replace('.fz', '')
    outname = outname.replace('.gz', '')

    assert(outname[-5:] == '.fits')

    outname = outname.replace('.fits', '_ps1.fits')

    if cube_index is not None:
        outname = outname.replace('.fits',
                                  '-' + str(cube_index).zfill(5) + '.fits')
    
    assert(not os.path.exists(outname))

    _atomic_write(ps1_matches, outname)
    
def gather_gaia_crossmatches(catalog, mjd=None, gfa_targets=None):

    # handle case of exposure with no retained sources
    if catalog is None:
        return None
    
    gaia_matches = gaia.gaia_xmatch(catalog['ra'], catalog['dec'], mjd=mjd,
                                    gfa_targets=gfa_targets)

    # avoid downstream conflict with 'ra', 'dec' columns that refer
    # to the world coordinates of the GFA detections
    gaia_matches.rename_column('ra', 'ra_gaia')
    gaia_matches.rename_column('dec', 'dec_gaia')

    if mjd is not None:
        gaia_matches['mjd_for_gaia'] = mjd

    catalog['gaia_motion_corr'] = int(mjd is not None)
        
    return gaia_matches

def append_gaia_crossmatches(catalog, mjd=None, gfa_targets=None):
    gaia_matches = gather_gaia_crossmatches(catalog, mjd=mjd,
                                            gfa_targets=gfa_targets)

    # handle case of exposure with no retained sources
    if gaia_matches is None:
        return None
    
    # I believe that there will always be a Gaia match for each
    # detected source, but will need to see if that assumption breaks
    # at any point

    catalog = hstack([catalog, gaia_matches])
    
    return catalog

def gather_pixel_stats(exp, skip=False):

    if not skip:
        print('Attempting to compute basic statistics of raw pixel data')
        t = None
        for extname, im in exp.images.items():
            if im is None:
                continue

            print('Computing pixel statistics for ' + extname)
            t_im = bis.compute_all_stats(im.image, extname=extname)
            if t is None:
                t = t_im
            else:
                t = vstack([t, t_im])
    else:
        t = Table()
        t['camera'] = list(exp.images.keys())
            
                
    return t
    
def high_level_ccds_metrics(tab, catalog, exp):
    
    nrows = len(tab)

    fwhm_major_pix = np.zeros(nrows)
    fwhm_minor_pix = np.zeros(nrows)
    fwhm_pix = np.zeros(nrows)
    fwhm_asec = np.zeros(nrows)
    n_sources = np.zeros(nrows, dtype=int)
    n_sources_for_shape = np.zeros(nrows, dtype=int)
    # boolean mask listing whether each source in 'catalog' was
    # used in computing its camera's overall FWHM_ASEC value
    if catalog is not None:
        used_for_fwhm_meas = np.zeros(len(catalog), dtype=bool)

    for i, row in enumerate(tab):
        if (catalog is None) or (np.sum(catalog['camera'] == row['camera']) == 0):
            fwhm_major_pix[i] = np.nan
            fwhm_minor_pix[i] = np.nan
            fwhm_pix[i] = np.nan
            fwhm_asec[i] = np.nan
            continue

        bad_amps = exp.images[row['camera']].overscan.bad_amps_list()
        fwhm_stats = bcs.overall_image_fwhm(catalog[catalog['camera'] == row['camera']], bad_amps=bad_amps)
        fwhm_major_pix[i] = fwhm_stats[0]
        fwhm_minor_pix[i] = fwhm_stats[1]
        fwhm_pix[i] = fwhm_stats[2]
        fwhm_asec[i] = fwhm_stats[3]
        n_sources[i] = int(np.sum(catalog['camera'] == row['camera']))
        n_sources_for_shape[i] = fwhm_stats[4]
        used_for_fwhm_meas[catalog['camera'] == row['camera']] = fwhm_stats[5]

    tab['fwhm_major_pix'] = fwhm_major_pix
    tab['fwhm_minor_pix'] = fwhm_minor_pix
    tab['fwhm_pix'] = fwhm_pix
    tab['fwhm_asec'] = fwhm_asec
    tab['n_sources'] = n_sources
    tab['n_sources_for_shape'] = n_sources_for_shape

    if catalog is not None:
        catalog['used_for_fwhm_meas'] = used_for_fwhm_meas.astype('uint8')

def prescan_overscan_ccds_table(tab, exp):
    # add information about bad pixels in overscan/prescan to
    # CCDs table, should be useful in identifying reasons for problematic
    # reductions...

    tab['npix_bad_total'] = [exp.images[extname].overscan.n_badpix_all for extname in tab['extname']]

    ampnames = common.valid_amps_list()
    npix_bad_per_amp = np.zeros((len(tab), len(ampnames)), dtype=int)
    overscan_medians = np.zeros((len(tab), len(ampnames)), dtype='float32')
    prescan_medians = np.zeros((len(tab), len(ampnames)), dtype='float32')
    
    for i, t in enumerate(tab):
        npix_bad_per_amp[i, :] = np.array([exp.images[t['extname']].overscan.n_badpix[amp] for amp in ampnames])
        overscan_medians[i, :] = np.array([exp.images[t['extname']].overscan.overscan_medians[amp] for amp in ampnames])
        prescan_medians[i, :] = np.array([exp.images[t['extname']].overscan.prescan_medians[amp] for amp in ampnames])

    tab['npix_bad_per_amp'] = npix_bad_per_amp
    tab['overscan_medians_adu'] = overscan_medians
    tab['prescan_medians_adu'] = prescan_medians

def astrom_ccds_table(tab, exp):
    # package WCS solutions into CCDs table

    nrows = len(tab)
    crvals = np.zeros((nrows, 2), dtype=float) # double
    naxis = np.zeros((nrows, 2), dtype=int)
    naxis[:, 0] = 2048 # prescan/overscan removed
    naxis[:, 1] = 1032 # prescan/overscan removed
    cds = np.zeros((nrows, 2, 2), dtype=float) # double
    cdelts = np.ones((nrows, 2), dtype=float) # double
    crpixs = np.zeros((nrows, 2), dtype=float) # double
    ctypes = np.zeros((nrows, 2), dtype='U8') # double
    ctypes[:, 0] = 'RA---TAN'
    ctypes[:, 1] = 'DEC--TAN'
    longpoles = np.zeros(nrows, dtype=float) + 180 # double
    latpoles = np.zeros(nrows, dtype=float) + 90 # double
    pv2s = np.zeros((nrows, 2), dtype=float) # double
    
    
    tab['NAXIS'] = naxis
    
    for i, extname in enumerate(tab['extname']):
        crvals[i, :] = exp.images[extname].wcs.wcs.crval
        cds[i, :, :] = np.transpose(exp.images[extname].wcs.wcs.cd)
        crpixs[i, :] = exp.images[extname].wcs.wcs.crpix

    tab['cd'] = cds
    tab['cdelt'] = cdelts
    tab['crpix'] = crpixs
    tab['crval'] = crvals
    tab['ctype'] = ctypes
    tab['longpole'] = longpoles
    tab['latpole'] = latpoles
    tab['pv2'] = pv2s


def radprof_ccds_table(tab, exp):
    # package PSF radial profiles into CCDs table
    nrows = len(tab)

    # number of radial profile radius bins
    nrad = 26 # special number -- eventually handle this better

    all_radii = np.zeros((nrows, nrad), dtype='float32')
    all_profiles = np.zeros((nrows, nrad), dtype='float32')

    for i, t in enumerate(tab):
        psf = exp.images[t['extname']].psf

        if psf is None:
            continue

        all_radii[i, :] = psf.profile_radius_pix
        all_profiles[i, :] = psf.radial_profile

    tab['profile_radius_pix'] = all_radii
    tab['psf_radial_profile'] = all_profiles

def dark_current_ccds_table(tab, exp):
    fname_master_dark = []
    nrows = len(tab)
    do_fit_dark_scaling = np.zeros(nrows, dtype='uint8')
    origtime = np.zeros(nrows)
    master_dark_gccdtemp = np.zeros(nrows)
    dark_temp_scaling_factor = np.zeros(nrows)
    total_dark_scaling_factor = np.zeros(nrows)
    rescale_factors_per_amp = np.zeros((nrows, 4), dtype=float)
    dark_rescale_factor_bestfit = np.zeros(nrows)
    dark_rescale_factor_adopted = np.zeros(nrows)
    apply_rescale_factor = np.zeros(nrows, dtype='uint8')
    dark_rescale_ncalls = np.zeros((nrows, 4), dtype=int) # this is per-amp
    dark_rescale_converged = np.zeros((nrows, 4), dtype='uint8')
    
    for i, t in enumerate(tab):
        dc = exp.dark_current_objs[t['extname']]
        fname_master_dark.append(dc.fname_master_dark)
        do_fit_dark_scaling[i] = dc.do_fit_dark_scaling
        origtime[i] = dc.header['ORIGTIME']
        master_dark_gccdtemp[i] = dc.header['GCCDTEMP']
        dark_temp_scaling_factor[i] = dc.temp_scaling_factor
        total_dark_scaling_factor[i] = dc.total_dark_scaling
        rescale_factors_per_amp[i, :] = dc.rescale_factors
        dark_rescale_factor_bestfit[i] = dc.dark_rescale_factor_bestfit
        apply_rescale_factor[i] = dc.apply_rescale_fac
        dark_rescale_factor_adopted[i] = dc.dark_rescale_factor_adopted
        dark_rescale_ncalls[i, :] = dc.dark_rescale_ncalls
        dark_rescale_converged[i, :] = dc.dark_rescale_converged
        
    tab['fname_master_dark'] = fname_master_dark
    tab['do_fit_dark_scaling'] = do_fit_dark_scaling
    tab['master_dark_exptime'] = origtime
    tab['master_dark_gccdtemp'] = master_dark_gccdtemp
    tab['dark_temp_scaling_factor'] = dark_temp_scaling_factor
    tab['total_dark_scaling_factor'] = total_dark_scaling_factor
    tab['dark_rescale_factors_per_amp'] = rescale_factors_per_amp
    tab['dark_rescale_factor_bestfit'] = dark_rescale_factor_bestfit
    tab['dark_rescale_factor_adopted'] = dark_rescale_factor_adopted
    tab['apply_dark_rescale_factor'] = apply_rescale_factor
    tab['dark_rescale_ncalls'] = dark_rescale_ncalls
    tab['dark_rescale_converged'] = dark_rescale_converged
    
def assemble_ccds_table(tab, catalog, exp, outdir, proc_obj, cube_index=None,
                        ps1=None, det_sn_thresh=5.0, sky_mags=True,
                        minimal=False, mjdrange=None):

    nrows = len(tab)

    tab['extname'] = tab['camera']

    tab['contrast'] = [exp.images[extname].header['CONTRAST'] for extname in tab['camera']]

    if minimal:
        return tab

    if sky_mags:
        tab['sky_mag_ab'] = [exp.images[extname].sky_mag for extname in tab['camera']]
        tab['sky_mag_ab_subregion'] = [exp.images[extname].sky_mag_upper for extname in tab['camera']]

    amps = common.valid_amps_list()

    if sky_mags:
        sky_mag_ab_per_amp = np.zeros((nrows, len(amps)), dtype='float32') # 4 amps
        for i, extname in enumerate(tab['camera']):
            for j in range(len(amps)):
                _mag = np.nan if exp.images[extname].sky_mag_per_amp is None else exp.images[extname].sky_mag_per_amp[j]
                sky_mag_ab_per_amp[i, j] = _mag

        tab['sky_mag_ab_per_amp'] = sky_mag_ab_per_amp
    
    tab['petal_loc'] = np.array([common.gfa_extname_to_gfa_number(extname) for extname in tab['camera']], dtype='uint8')

    tab['expid'] = [exp.images[extname].header['EXPID'] for extname in tab['camera']]

    # should work except if early versions of guide cubes lacked
    # MJD information...
    tab['mjd'] = [exp.images[extname].try_retrieve_meta_keyword('MJD-OBS', placeholder=0.0) for extname in tab['camera']]

    eph = util.load_lst()
    tab['lst_deg'] = [util.interp_ephemeris(t['mjd'], eph=eph) for t in tab]
    tab['moon_illumination'] = [util.interp_ephemeris(t['mjd'], eph=eph, mphase=True) for t in tab]

    tab['program'] = [str(exp.images[extname].try_retrieve_meta_keyword('PROGRAM', placeholder='')) for extname in tab['camera']]

    tab['skyra'] = [exp.images[extname].try_retrieve_meta_keyword('SKYRA', placeholder=np.nan) for extname in tab['camera']]

    tab['skydec'] = [exp.images[extname].try_retrieve_meta_keyword('SKYDEC', placeholder=np.nan) for extname in tab['camera']]
    
    # zenith distance using approximate center of the field given by SKYRA, SKYDEC

    tab['zenith_dist_deg'] = [util._zenith_distance(t['skyra'], t['skydec'], t['lst_deg']) for t in tab]
    
    tab['domshutl'] = np.array(nrows*[exp.try_retrieve_header_card('DOMSHUTL', placeholder='')], dtype='U8')
    tab['domshutu'] = np.array(nrows*[exp.try_retrieve_header_card('DOMSHUTU', placeholder='')], dtype='U8')
    tab['pmcover'] = exp.try_retrieve_header_card('PMCOVER', placeholder='')
    tab['moonra'] = exp.try_retrieve_header_card('MOONRA', placeholder=np.nan)
    tab['moondec'] = exp.try_retrieve_header_card('MOONDEC', placeholder=np.nan)

    tab['moon_zd_deg'] = util._zenith_distance(tab['moonra'][0],
                                               tab['moondec'][0],
                                               tab['lst_deg'][0])
    
    tab['t_c_for_dark'] = [exp.images[extname].t_c_for_dark for extname in tab['camera']]
    tab['t_c_for_dark_is_guess'] = [int(exp.images[extname].t_c_for_dark_is_guess) for extname in tab['camera']]
    tab['time_s_for_dark'] = [exp.images[extname].time_s_for_dark for extname in tab['camera']]

    tab['night'] = exp.try_retrieve_header_card('NIGHT', placeholder='')

    tab['focus'] = exp.try_retrieve_header_card('FOCUS', placeholder='')

    tab['exptime'] = [exp.images[extname].try_retrieve_meta_keyword('EXPTIME', placeholder=np.nan) for extname in tab['camera']]
    
    tab['cube_index'] = np.nan if cube_index is None else int(cube_index)

    if cube_index == -1:
        tab['coadd_index_start'] = [exp.images[extname].coadd_index_range[0] for extname in tab['camera']]
        tab['coadd_index_end'] = [exp.images[extname].coadd_index_range[1] for extname in tab['camera']]
        tab['coadd_mjdobs_min'] = [exp.bintables[extname]['MJD-OBS'][ind] for extname, ind in zip(tab['camera'], tab['coadd_index_start'])]
        tab['coadd_mjdobs_max'] =  [exp.bintables[extname]['MJD-OBS'][ind] for extname, ind in zip(tab['camera'], tab['coadd_index_end'])]
    
    tab['racen'] = np.zeros(len(tab), dtype=float)
    tab['deccen'] = np.zeros(len(tab), dtype=float)

    tab['fname_raw'] = proc_obj.fname_in
    tab['gitrev']  = proc_obj.gitrev

    tab['fiber_fracflux'] = [(exp.images[extname].psf.fiber_fracflux if exp.images[extname].psf is not None else np.nan) for extname in tab['camera']]

    tab['fiber_fracflux_elg'] = [(exp.images[extname].psf.fiber_fracflux_elg if exp.images[extname].psf is not None else np.nan) for extname in tab['camera']]

    tab['n_sources_for_psf'] = [(exp.images[extname].psf.nstars if exp.images[extname].psf is not None else 0) for extname in tab['camera']]
    
    # this pertains to aperture _3 which is 1.5 asec radius
    tab['aper_corr_fac'] = [(exp.images[extname].psf.aper_corr_fac if exp.images[extname].psf is not None else np.nan) for extname in tab['camera']]

    tab['xcentroid_psf'] = [(exp.images[extname].psf.xcen_flux_weighted if exp.images[extname].psf is not None else np.nan) for extname in tab['camera']]
    tab['ycentroid_psf'] = [(exp.images[extname].psf.ycen_flux_weighted if exp.images[extname].psf is not None else np.nan) for extname in tab['camera']]

    tab['psf_fwhm_pix'] =  [(exp.images[extname].psf.moffat_fwhm_pix if exp.images[extname].psf is not None else np.nan) for extname in tab['camera']]

    tab['psf_fwhm_asec'] = [(exp.images[extname].psf.moffat_fwhm_asec if exp.images[extname].psf is not None else np.nan) for extname in tab['camera']]

    tab['psf_centroid_cbox'] = [(float(exp.images[extname].psf.cbox) if exp.images[extname].psf is not None else np.nan) for extname in tab['camera']]

    tab['psf_centroid_failed'] =  [(exp.images[extname].psf.psf_centroiding_failed if exp.images[extname].psf is not None else 0) for extname in tab['camera']]

    tab['radprof_fwhm_asec'] =  [(exp.images[extname].psf.radprof_fwhm_asec if exp.images[extname].psf is not None else np.nan) for extname in tab['camera']]

    # is 0 the best placeholder value here when no PSF exists?
    tab['psf_centroiding_flag'] =  [(exp.images[extname].psf.psf_centroiding_flag if exp.images[extname].psf is not None else 0) for extname in tab['camera']]

    tab['psf_asymmetry_ratio'] = [(exp.images[extname].psf.psf_asymmetry_ratio if exp.images[extname].psf is not None else np.float32(np.nan)) for extname in tab['camera']]

    tab['psf_asymmetry_numerator'] = [(exp.images[extname].psf.psf_asymmetry_numerator if exp.images[extname].psf is not None else np.float32(np.nan)) for extname in tab['camera']]

    tab['psf_asymmetry_denominator'] =  [(exp.images[extname].psf.psf_asymmetry_denominator if exp.images[extname].psf is not None else np.float32(np.nan)) for extname in tab['camera']]

    tab['psf_total_flux'] =  [(exp.images[extname].psf.psf_total_flux if exp.images[extname].psf is not None else np.float32(np.nan)) for extname in tab['camera']]

    radprof_ccds_table(tab, exp)

    for i, extname in enumerate(tab['camera']):
        racen, deccen = ccd_center_radec(exp.images[extname].wcs)
        tab['racen'][i] = racen
        tab['deccen'][i] = deccen

    tab['mountha_header'] = exp.try_retrieve_header_card('MOUNTHA', placeholder=np.nan)
    tab['mountdec_header'] = exp.try_retrieve_header_card('MOUNTDEC', placeholder=np.nan)

    tab['ha_deg'] = [util._get_ha(t['skyra'], t['lst_deg'], t['mountdec_header']) for t in tab]

    tab['ha_deg_per_gfa'] = [util._get_ha(t['racen'], t['lst_deg'], t['mountdec_header']) for t in tab]
    
    tab['moon_sep_deg'] = util.moon_separation(tab['moonra'], tab['moondec'],
                                               tab['racen'], tab['deccen'])

    # per-camera zenith distance -- will only be accurate to the extent that
    # each camera's WCS recalibration succeeded
    tab['zd_deg_per_gfa'] = [util._zenith_distance(t['racen'], t['deccen'], t['lst_deg']) for t in tab]

    tab['header_airmass'] = exp.try_retrieve_header_card('AIRMASS', placeholder=np.nan)

    # this should make the airmass properly evolve
    # with time for guide cube outputs
    tab['airmass'] = 1.0/np.cos(tab['zenith_dist_deg']/(180.0/np.pi))

    # this per-camera version of airmass should also evolve
    # properly with time for guide cube outputs
    tab['airmass_per_gfa'] = 1.0/np.cos(tab['zd_deg_per_gfa']/(180.0/np.pi))
    
    tab['zp_adu_per_s'] = [exp.images[extname].compute_zeropoint(ps1) for extname in tab['camera']]

    tab['transparency'] = [util.transparency_from_zeropoint(tab[i]['zp_adu_per_s'], tab[i]['airmass_per_gfa'], tab[i]['camera']) for i in range(len(tab))]

    par = common.gfa_misc_params()
    tab['kterm'] = np.float32(par['kterm'])

    tab['det_sn_thresh'] = det_sn_thresh
    
    prescan_overscan_ccds_table(tab, exp)
    high_level_ccds_metrics(tab, catalog, exp)
    astrom_ccds_table(tab, exp)
    dark_current_ccds_table(tab, exp)

    if mjdrange is not None:
        tab['req_mjd_min'] = mjdrange[0]
        tab['req_mjd_max'] = mjdrange[1]

    return tab

def write_ccds_table(tab, outdir, proc_obj, cube_index=None):

    assert(os.path.exists(outdir))

    outname = os.path.join(outdir, os.path.basename(proc_obj.fname_in))

    # get rid of any ".fz" or ".gz" present in input filename
    outname = outname.replace('.fz', '')
    outname = outname.replace('.gz', '')

    assert(outname[-5:] == '.fits')

    outname = outname.replace('.fits', '_ccds.fits')

    if cube_index is not None:
        outname = outname.replace('.fits',
                                  '-' + str(cube_index).zfill(5) + '.fits')

    assert(not os.path.exists(outname))
    
    print('Attempting to write CCDs table to ' + outname)

    _atomic_write(tab, outname)

def write_psfs(exp, outdir, fname_in, cube_index=None, cubes=False):

    assert(os.path.exists(outdir))

    outname = os.path.join(outdir, os.path.basename(fname_in))

    # get rid of any ".fz" or ".gz" present in input filename
    outname = outname.replace('.fz', '')
    outname = outname.replace('.gz', '')

    assert(outname[-5:] == '.fits')

    suffix = '_psfcubes.fits' if cubes else '_psfs.fits'

    outname = outname.replace('.fits', suffix)

    if cube_index is not None:
        outname = outname.replace('.fits',
                                  '-' + str(cube_index).zfill(5) + '.fits')
    
    assert(not os.path.exists(outname))

    hdul = []
    for image in exp.images.values():
        if image is None:
            continue
        if image.psf is not None:
            is_primary = (len(hdul) == 0)

            if not cubes:
                hdu = image.psf.to_hdu(primary=is_primary)
            else:
                hdu = image.psf.cube_to_hdu(primary=is_primary)

            hdul.append(hdu)

    if len(hdul) > 0:
        hdul = fits.HDUList(hdul)
        _atomic_write(hdul, outname)

def write_full_detlists(exp, outdir, fname_in, cube_index=None):

    assert(os.path.exists(outdir))

    outname = os.path.join(outdir, os.path.basename(fname_in))

    # get rid of any ".fz" or ".gz" present in input filename
    outname = outname.replace('.fz', '')
    outname = outname.replace('.gz', '')

    assert(outname[-5:] == '.fits')

    outname = outname.replace('.fits', '_detlist.fits')

    if cube_index is not None:
        outname = outname.replace('.fits',
                                  '-' + str(cube_index).zfill(5) + '.fits')

    tables = [image.full_detlist for image in exp.images.values()]

    tab = vstack(tables)

    _atomic_write(tab, outname)

def write_dm_fieldmodel(fm, outdir, fname_in, cube_index=None):
    # fm should be a desimeter fieldmodel object

    if fm is None:
        print('No desimeter field model JSON file will be written')
        return

    assert(os.path.exists(outdir))

    outname = os.path.join(outdir, os.path.basename(fname_in))

    # get rid of any ".fz" or ".gz" present in input filename
    outname = outname.replace('.fz', '')
    outname = outname.replace('.gz', '')

    assert(outname[-5:] == '.fits')

    outname = outname.replace('.fits', '_fieldmodel.json')

    if cube_index is not None:
        outname = outname.replace('.json',
                                  '-' + str(cube_index).zfill(5) + '.json')

    outname_tmp = outname + '.tmp'
    with open(outname_tmp, 'w') as file:
        _json = json.loads(fm.tojson())
        _json['desimeter_gitrev'] = fm.desimeter_gitrev
        _json['n_cameras'] = fm.n_cameras
        file.write(json.dumps(_json))

    os.rename(outname_tmp, outname)
