import gfa_reduce.xmatch.gaia as gaia_xmatch
import astropy.io.fits as fits
import os
import healpy
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
from astropy import wcs
from .amm_2dhist import amm_2dhist
from scipy.ndimage import gaussian_filter
from .center_contrast import center_contrast
import gfa_reduce.common as common
import copy
from gfa_reduce.gfa_wcs import nominal_tan_wcs
from desiutil.log import get_logger


def downselected_star_sample(cat, n_desi_max):
    assert(len(np.unique(cat['camera'])) == 1)

    if 'DQ_FLAGS' in cat.colnames:
        _cat = cat[cat['DQ_FLAGS'] == 0]
    else:
        _cat = cat[cat['dq_flags'] == 0]

    if 'APER_SUM_BKGSUB_3' in cat.colnames:
        brightness = _cat['APER_SUM_BKGSUB_3']
    elif 'aper_sum_bkgsub_3' in cat.colnames:
        brightness = _cat['aper_sum_bkgsub_3']
    elif 'DETMAP_PEAK' in cat.colnames:
        brightness = _cat['DETMAP_PEAK']
    else:
        brightness = _cat['detmap_peak']

    sind = np.argsort(-1.0*brightness)

    n = len(sind)

    result = _cat[sind[0:(min(n_desi_max, n))]]

    return result


def gaia_cat_for_exp(racen, deccen, mjd=None):

    nside = 32
    ipix = healpy.pixelfunc.get_all_neighbours(nside, racen, phi=deccen,
                                               lonlat=True)

    # remove dummy -1 values from ipix

    ipix = ipix[ipix >= 0]

    # probably also want to include the neighbors to the neighbors
    ipix_all = []
    for _ipix in ipix:
        ipix_all.append(healpy.pixelfunc.get_all_neighbours(nside,
                                                            _ipix))

    ipix = np.unique(ipix_all)

    ipix = ipix[ipix >= 0]

    ra_pixcenters, dec_pixcenters = healpy.pixelfunc.pix2ang(nside, ipix,
                                                              lonlat=True)

    gaia = gaia_xmatch.read_gaia_cat(ra_pixcenters, dec_pixcenters, mjd=mjd)
    return gaia


def pattern_match(catalog, skyra, skydec, extname, gaia, arcmin_max):

    # cat needs to have fields xcentroid and ycentroid
    # skyra, skydec are the initial guesses of the
    # actual center of the field of view; these need to be within
    # arcmin_max of the true FOV center for this routine to succeed

    # output needs to include, at a minimum:
    #     xshift_best
    #     yshift_best
    #     contrast
    #     extname
    # probably also want
    #     expid retrieved from the catalog table
    log = get_logger()
    cat = copy.deepcopy(catalog)
    # should just entirely rename racen, deccen to skyr, skydec ...
    racen = skyra
    deccen = skydec

    # this apparently happened for the CI in some cases...
    if isinstance(racen, str) or isinstance(deccen, str):
        return None

    if gaia is None:
        gaia = gaia_cat_for_exp(racen, deccen)

    g = SkyCoord(gaia['ra']*u.deg, gaia['dec']*u.deg)
    c = SkyCoord(racen*u.deg, deccen*u.deg)

    dangle = g.separation(c)

    gaia = gaia[dangle.degree < 2]
    g = SkyCoord(gaia['ra']*u.deg, gaia['dec']*u.deg)

    cat = cat[cat['camera'] == extname]

    assert(len(cat) > 0)

    n_desi_max = 150

    if len(cat) > n_desi_max:
        cat = downselected_star_sample(cat, n_desi_max)

    par = common.gfa_misc_params()

    astrom = nominal_tan_wcs(racen, deccen, extname)

    x_gaia_guess, y_gaia_guess = astrom.wcs_world2pix(gaia['ra'], gaia['dec'],
                                                      0)
    dx_all = np.array([], dtype='float64')
    dy_all = np.array([], dtype='float64')

    cat = cat[np.argsort(cat['dec'])] # not really necessary

    _g = None
    for i in range(len(cat)):
        log.debug('%d of %d', i+1, len(cat))
        c = SkyCoord(cat[i]['ra']*u.deg, cat[i]['dec']*u.deg)
        dangle = c.separation(g if _g is None else _g)
        w = (np.where(dangle.arcminute < arcmin_max))[0]
        if len(w) == 0:
            continue

        assert(len(dangle) == len(x_gaia_guess))
        assert(len(dangle) == len(y_gaia_guess))

        dx = cat[i]['xcentroid'] - x_gaia_guess[w]
        dy = cat[i]['ycentroid'] - y_gaia_guess[w]
        dx_all = np.concatenate((dx_all, dx))
        dy_all = np.concatenate((dy_all, dy))

        if _g is None:
            padfac = 2.0
            keep = (dangle.arcminute < (8.1 + padfac*arcmin_max))
            x_gaia_guess = x_gaia_guess[keep]
            y_gaia_guess = y_gaia_guess[keep]
            _g = g[keep]

    assert(len(dx_all) == len(dy_all))

    axlim = max(np.round(arcmin_max*60.0/0.214), 1000.0)

    dx = 1.0
    dy = 1.0
    nx = ny = 2*axlim
    counts, x_edges_left, y_edges_left = amm_2dhist(-1.0*axlim, -1.0*axlim,
                                                    nx, ny, dx, dy,
                                                    dx_all, dy_all)

    counts = counts.astype(float)

    # it's possible that this fwhm_pix value could use more tailoring
    fwhm_pix = 4.7
    sigma_pix = fwhm_pix/(2*np.sqrt(2*np.log(2)))
    smth = gaussian_filter(counts, sigma_pix, mode='constant')

    counts_shape = counts.shape
    sidelen = counts_shape[0]

    indmax = np.argmax(smth)

    indx = (indmax // sidelen).astype(int)
    indy = (indmax % sidelen).astype(int)

    xshift_best = x_edges_left[indx] + 0.5*dx
    yshift_best = y_edges_left[indy] + 0.5*dy

    # assumes dx = dy = 1 !!
    ycen_grid, xcen_grid = np.meshgrid(x_edges_left[0:(len(x_edges_left)-1)], y_edges_left[0:(len(y_edges_left)-1)])

    d = np.sqrt(np.power(xcen_grid - xshift_best, 2) + np.power(ycen_grid - yshift_best, 2))

    r_max = 4*4.7 # roughly 4 asec radius

    sind = np.argsort(np.ravel(smth))

    val_90 = np.ravel(smth)[sind[int(np.round(0.925*len(sind)))]]

    high = [(d < r_max) & ((smth == np.max(smth)) | (smth > val_90))]

    assert(np.sum(high) > 0)

    wt = np.sum(high*smth)
    xshift_best = np.sum(high*xcen_grid*smth)/wt
    yshift_best = np.sum(high*ycen_grid*smth)/wt

    contrast = center_contrast(smth)

    result = {'xshift_best': xshift_best,
              'yshift_best': yshift_best,
              'contrast': contrast,
              'extname': extname,
              'astr_guess': astrom,
              'det_ids_used': np.array(cat['det_id'])}

    return result
