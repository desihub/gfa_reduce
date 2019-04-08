from ci_reduce.image import CI_image
from ci_reduce.exposure import CI_exposure
import ci_reduce.common as common
import ci_reduce.xmatch.gaia as gaia
import astropy.io.fits as fits
from astropy.table import vstack, hstack
import os
import ci_reduce.analysis.basic_image_stats as bis
import ci_reduce.analysis.basic_catalog_stats as bcs
import numpy as np
import time

# in the context of this file, "image" and "exposure" generally refer to 
# CI_image and CI_exposure objects

def loading_image_extension_message(extname):
    print('Attempting to load image extension : ' + extname)

def load_image_from_hdu(hdu, verbose=True):
    loading_image_extension_message(hdu.header['EXTNAME'])

    if verbose:
        print(repr(hdu.header))

    return CI_image(hdu.data, hdu.header)

def load_image_from_filename(fname, extname):
    assert(os.path.exists(fname))

    loading_image_extension_message(extname)
    assert(common.is_valid_image_extname(extname))

    data, header = fits.getdata(fname, extname=extname, header=True)
    return CI_image(data, header)

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
        except:
            print('encountered problem reading ' + fname)
            time.sleep(delay)
        if hdul is not None:
            break

    # die if unable to read file after max_attempts attempts
    assert(hdul is not None)

    return hdul

def load_exposure(fname, verbose=True, realtime=False):
    assert(os.path.exists(fname))

    print('Attempting to load exposure : ' + fname)

    par = common.ci_misc_params()

    if not realtime:
        hdul = fits.open(fname)
    else:
        hdul = realtime_raw_read(fname)

    dummy_fz_header = None

    is_image_hdu = np.zeros(len(hdul), dtype=bool)
    for i, hdu in enumerate(hdul):
        # real data has another dummy extension added with no EXTNAME
        keywords = [c[0] for c in hdu.header.cards]
        if not ('EXTNAME' in keywords):
            continue
        if (hdu.header['EXTNAME']).strip() == par['fz_dummy_extname']:
            dummy_fz_header = hdu.header
            continue
        is_image_hdu[i] = True

    w_im = np.where(is_image_hdu)[0]

    imlist = [load_image_from_hdu(hdul[ind], verbose=verbose) for ind in w_im]

    exp = CI_exposure(imlist, dummy_fz_header=dummy_fz_header)

    print('Successfully loaded exposure : ' + fname)
    print('Exposure has ' + str(exp.num_images_populated()) + 
          ' image extensions populated')
    print('Populated image extension names are : ' + 
          str(exp.populated_extnames()))

    return exp

def reduced_image_fname(outdir, fname_in, flavor, gzip=True):
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

    assert(not os.path.exists(outname))

    return outname

def check_image_level_outputs_exist(outdir, fname_in, gzip=True):
    par = common.ci_misc_params()

    for flavor in par['reduced_image_flavors']:
        _ = reduced_image_fname(outdir, fname_in, flavor, gzip=gzip)

def retrieve_git_rev():
    code_dir = os.path.dirname(os.path.realpath(__file__))
    cwd = os.getcwd()
    do_chdir = (cwd[0:len(code_dir)] != code_dir)
    if do_chdir:
        os.chdir(code_dir)
    gitrev = os.popen("git rev-parse --short HEAD").read().replace('\n','')
    if do_chdir:
        os.chdir(cwd)
    print('"git rev" version info:', gitrev)

    return gitrev

def write_image_level_outputs(exp, outdir, fname_in, gzip=True):
    # exp is a CI_exposure object
    # outdir is the output directory (string)
    # fname_in is the input filename (string)

    par = common.ci_misc_params()

    for flavor in par['reduced_image_flavors']:
        outname = reduced_image_fname(outdir, fname_in, flavor, gzip=gzip)

        hdulist = exp.to_hdulist(flavor=flavor)

        print('Attempting to write ' + flavor + ' image output to ' + 
              outname)

        hdulist.writeto(outname)

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
    # catalogs is the output of CI_exposure's all_source_catalogs() method
    # which is a dictionary of astropy QTable's, with the keys
    # being the CI camera extension names

    # want to add a column to each table giving the CI camera name, then
    # append the all into one per-exposure table

    assert(type(catalogs).__name__ == 'dict')

    for extname, tab in catalogs.items():
        tab['camera'] = extname
        tab['ci_number'] = [common.ci_extname_to_ci_number(extname) for extname in tab['camera']]

    composite = vstack([tab for tab in catalogs.values()])
    composite = strip_none_columns(composite)

    return composite

def write_exposure_source_catalog(catalog, outdir, fname_in):

    assert(os.path.exists(outdir))

    outname = os.path.join(outdir, os.path.basename(fname_in))

    # get rid of any ".fz" or ".gz" present in input filename
    outname = outname.replace('.fz', '')
    outname = outname.replace('.gz', '')

    assert(outname[-5:] == '.fits')

    outname = outname.replace('.fits', '_catalog.fits')

    assert(not os.path.exists(outname))

    print('Attempting to write source catalog to ' + outname)
    catalog.write(outname, format='fits')

def gather_gaia_crossmatches(catalog):
    gaia_matches = gaia.gaia_xmatch(catalog['ra'], catalog['dec'])

    # avoid downstream conflict with 'ra', 'dec' columns that refer
    # to the world coordinates of the CI detections
    gaia_matches.rename_column('ra', 'ra_gaia')
    gaia_matches.rename_column('dec', 'dec_gaia')

    return gaia_matches

def append_gaia_crossmatches(catalog):
    gaia_matches = gather_gaia_crossmatches(catalog)

    # I believe that there will always be a Gaia match for each
    # detected source, but will need to see if that assumption breaks
    # at any point

    catalog = hstack([catalog, gaia_matches])
    
    return catalog

def gather_pixel_stats(exp):

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

    return t

def high_level_ccds_metrics(tab, catalog):

    nrows = len(tab)

    fwhm_major_pix = np.zeros(nrows)
    fwhm_minor_pix = np.zeros(nrows)
    fwhm_pix = np.zeros(nrows)
    fwhm_asec = np.zeros(nrows)
    n_sources = np.zeros(nrows, dtype=int)
    n_sources_for_shape = np.zeros(nrows, dtype=int)

    for i, row in enumerate(tab):
        fwhm_stats = bcs.overall_image_fwhm(catalog[catalog['camera'] == row['camera']])
        fwhm_major_pix[i] = fwhm_stats[0]
        fwhm_minor_pix[i] = fwhm_stats[1]
        fwhm_pix[i] = fwhm_stats[2]
        fwhm_asec[i] = fwhm_stats[3]
        n_sources[i] = int(np.sum(catalog['camera'] == row['camera']))
        n_sources_for_shape[i] = fwhm_stats[4]

    tab['fwhm_major_pix'] = fwhm_major_pix
    tab['fwhm_minor_pix'] = fwhm_minor_pix
    tab['fwhm_pix'] = fwhm_pix
    tab['fwhm_asec'] = fwhm_asec
    tab['n_sources'] = n_sources
    tab['n_sources_for_shape'] = n_sources_for_shape

def write_ccds_table(tab, catalog, exp, outdir, fname_in):

    assert(os.path.exists(outdir))

    outname = os.path.join(outdir, os.path.basename(fname_in))

    # get rid of any ".fz" or ".gz" present in input filename
    outname = outname.replace('.fz', '')
    outname = outname.replace('.gz', '')

    assert(outname[-5:] == '.fits')

    outname = outname.replace('.fits', '_ccds.fits')

    assert(not os.path.exists(outname))

    tab['sky_mag_ab'] = [im.sky_mag for im in exp.images.values()]
    tab['ci_number'] = [common.ci_extname_to_ci_number(extname) for extname in tab['camera']]

    high_level_ccds_metrics(tab, catalog)

    print('Attempting to write CCDs table to ' + outname)
    tab.write(outname, format='fits')
