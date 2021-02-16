import gfa_reduce.analysis.util as util
from astropy.stats import mad_std
from scipy.ndimage.filters import gaussian_filter
from scipy import ndimage
import numpy as np
from scipy.ndimage.measurements import label, find_objects
from astropy.table import Table
from photutils.centroids.core import fit_2dgaussian
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import CircularAperture, CircularAnnulus, EllipticalAperture
import gfa_reduce.common as common
import gfa_reduce.analysis.djs_maskinterp as djs_maskinterp
from gfa_reduce.analysis.djs_photcen import _loop_djs_photcen
import photutils
import copy
import os
import astropy.io.fits as fits

def slices_to_table(slices, detsn, extname):
    nslc = len(slices)
    xcen = np.zeros(nslc)
    ycen = np.zeros(nslc)
    xmin = np.zeros(nslc, dtype='int')
    xmax = np.zeros(nslc, dtype='int')
    ymin = np.zeros(nslc, dtype='int')
    ymax = np.zeros(nslc, dtype='int')

    detmap_peak = np.zeros(nslc)

    for i, slc in enumerate(slices):
        ycen[i] = (float(slc[0].start) + float(slc[0].stop))/2.0
        xcen[i] = (float(slc[1].start) + float(slc[1].stop))/2.0
        xmin[i] = slc[1].start
        xmax[i] = slc[1].stop
        ymin[i] = slc[0].start
        ymax[i] = slc[0].stop
        detmap_peak[i] = np.max(detsn[slc])

    detmap_sidelen_x = xmax - xmin # no + 1 here due to slice convention
    detmap_sidelen_y = ymax - ymin # no + 1 here due to slice convention

    tab = Table([xcen, ycen, xmin, xmax, ymin, ymax, detmap_sidelen_x, detmap_sidelen_y, detmap_peak], names=('xcen_init', 'ycen_init', 'detmap_xmin', 'detmap_xmax', 'detmap_ymin', 'detmap_ymax', 'detmap_sidelen_x', 'detmap_sidelen_y', 'detmap_peak'))

    return tab

def detmap_centroids(tab, detmap, max_cbox=31):

    # tab should be a table of detections like that output by the
    # slices_to_table function
    
    # construct lists of starting guess (x, y) pixel coordinates
    # and cbox sizes

    # cbox ceiling at 25 pixels?
    # cbox floor at 7 pixels?

    cbox = (tab['detmap_sidelen_x'] + tab['detmap_sidelen_y'])/ 2.0
    cbox_upper = max_cbox
    cbox_lower = 9
    cbox = np.maximum(cbox, cbox_lower)
    cbox = np.minimum(cbox, cbox_upper)

    assert(np.sum(cbox > cbox_upper) == 0)
    assert(np.sum(cbox < cbox_lower) == 0)

    results = _loop_djs_photcen(tab['xcen_init'], tab['ycen_init'], detmap,
                                cbox=cbox)

    tab['detmap_cbox'] = cbox
    tab['xcen_detmap_fw'] = results['x_djs']
    tab['ycen_detmap_fw'] = results['y_djs']

    return tab

def _get_area_from_ap(ap):
    # this is to try and work around the photutils API change
    # between versions 0.6 and 0.7
    if photutils.__version__.find('0.7') != -1:
        area = ap.area # 0.7
    else:
        area = ap.area() # 0.6

    return area


def detection_map(im, fwhm):
    psf_sigma = fwhm/2.355

    smth = gaussian_filter(ndimage.median_filter(im, 3), psf_sigma)

    smth[0, :] = 0.0
    smth[1031, :] = 0.0
    smth[:, 0] = 0.0
    smth[:, 2047] = 0.0
    
    sig_smth = mad_std(smth)

    detsn = (smth-np.median(smth))/sig_smth
    
    return detsn

def detect_sources(detsn, thresh):
    peaks = (detsn > thresh)
    blobs, _ = label(peaks)
    slices = find_objects(blobs)
    return slices

def aper_phot_unc_map(ivar):

    regfac = 0.01
    return np.power(ivar + (ivar == 0)*regfac*np.mean(ivar), -0.5)

def aper_rad_pix(extname):
    # list of aperture radii in pixels
    # list in asec comes from http://legacysurvey.org/dr7/catalogs

    # eventually may want to store this somewhere more useful rather
    # than hardcoding right here
    rad_asec = [0.5, 0.75, 1.0, 1.5, 2.0, 3.5, 5.0, 7.0]

    asec_per_pix = util.nominal_pixel_sidelen_arith()

    rad_pix = [r/asec_per_pix for r in rad_asec]

    return rad_pix

def do_aper_phot(data, catalog, extname, ivar_adu, sig_adu=None):
    # catalog should be the catalog with refined centroids 
    # for **one GFA camera**

    if sig_adu is None:
        sig_adu = aper_phot_unc_map(ivar_adu)

    print('Attempting to do aperture photometry')
    positions = list(zip(catalog['xcentroid'], catalog['ycentroid']))

    radii = aper_rad_pix(extname)

    apertures = [CircularAperture(positions, r=r) for r in radii]
    annulus_apertures = CircularAnnulus(positions, r_in=60.0, r_out=65.0)
    annulus_masks = annulus_apertures.to_mask(method='center')

    par = common.gfa_misc_params()

    b_over_a = par['nominal_mer_cd']/par['nominal_sag_cd']

    apertures_ell = [EllipticalAperture(positions, a, a*b_over_a, theta=np.pi/2) for a in radii]

    # 107 um fiber diam, 9 um on a side for a pixel
    # fiber diam from Table 4.1 of https://arxiv.org/abs/1611.00037
    rad_fiber_pix_sag = (107.0/9.0)/2.0
    deg_to_normal = 5.43 # [desi-commiss 522]

    rad_fiber_pix_mer = rad_fiber_pix_sag*np.sin(deg_to_normal/(180.0/np.pi))

    aper_fib = EllipticalAperture(positions, rad_fiber_pix_sag, rad_fiber_pix_mer, theta=np.pi/2)

    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        # this sigma_clipped_stats call is actually the slow part !!
        _, median_sigclip, std_bg = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)

    bkg_median = np.array(bkg_median)
    phot = aperture_photometry(data, apertures, 
                               error=sig_adu)

    for i, aperture in enumerate(apertures):
        aper_bkg_tot = bkg_median*_get_area_from_ap(aperture)
        catalog['aper_sum_bkgsub_' + str(i)] = phot['aperture_sum_' + str(i)] - aper_bkg_tot

        catalog['aper_bkg_' + str(i)] = aper_bkg_tot
        catalog['aperture_sum_err_' + str(i)] = phot['aperture_sum_err_' + str(i)]

    ###
    del phot
    phot = aperture_photometry(data, apertures_ell, 
                               error=sig_adu)
    for i, aperture in enumerate(apertures_ell):
        aper_bkg_tot = bkg_median*_get_area_from_ap(aperture)
        catalog['aper_ell_sum_bkgsub_' + str(i)] = phot['aperture_sum_' + str(i)] - aper_bkg_tot

        catalog['aper_ell_bkg_' + str(i)] = aper_bkg_tot
        catalog['aperture_ell_sum_err_' + str(i)] = phot['aperture_sum_err_' + str(i)]
    ###

    ###
    del phot
    phot = aperture_photometry(data, aper_fib, 
                               error=sig_adu)

    aper_bkg_tot = bkg_median*_get_area_from_ap(aper_fib)
    catalog['aper_sum_bkgsub_fib'] = phot['aperture_sum'] - aper_bkg_tot

    catalog['aper_bkg_fib'] = aper_bkg_tot
    catalog['aperture_sum_err_fib'] = phot['aperture_sum_err']

    ####

    # is .area() result a vector or scalar ??
    catalog['sky_annulus_area_pix'] = _get_area_from_ap(annulus_apertures)
    catalog['sky_annulus_median'] = bkg_median

def get_nominal_fwhm_pix(extname):
    # roughly 1.25 asec divided by GFA pixel scale
    # 1.25 is same nominal FWHM I used for source detection with CI reductions
    nominal_fwhm_pix = 6.1

    return nominal_fwhm_pix

def refine_centroids(tab, image, bitmask, ivar_adu, sig_adu=None,
                     skip_2dg=False):
    # input table tab gets augmented with additional columns

    print('Attempting to refine initial centroids')

    # this is for field acquisition mode
    if skip_2dg:
        print('Skipping 2D Gaussian source fitting')
        tab['xcentroid'] = tab['xcen_detmap_fw']
        tab['ycentroid'] = tab['ycen_detmap_fw']
        tab['sig_major_pix'] = 5.0 # HACK !!!
        tab['sig_minor_pix'] = 5.0 # HACK !!!
        return

    if sig_adu is None:
        aper_phot_unc_map(ivar_adu)

    # could scale this based on the typical size of the slices
    boxsize = 11

    half = int(np.floor(boxsize/2))

    med = np.median(image)

    nobj = len(tab)
    xcentroid = np.zeros(nobj)
    ycentroid = np.zeros(nobj)
    sig_major_pix = np.zeros(nobj)
    sig_minor_pix = np.zeros(nobj)
    ellipticity = np.zeros(nobj)
    # deg CC from +X, in the range -90 deg to +90 deg
    pos_angle = np.zeros(nobj)
    dxcentroid = np.zeros(nobj)
    dycentroid = np.zeros(nobj)

    no_centroid_refinement = np.zeros(nobj, dtype=bool)
    centroid_refinement_fail = np.zeros(nobj, dtype=bool)

    for i in range(nobj):
        ix_guess = int(round(tab[i]['xcen_detmap_fw']))
        iy_guess = int(round(tab[i]['ycen_detmap_fw']))
        min_edge_dist = util.min_edge_dist_pix(ix_guess, iy_guess)

        # if the centering box extends outside of image, don't
        # try to do any refinement of the centroid
        if min_edge_dist < half:
            xcentroid[i] = tab[i]['xcen_detmap_fw']
            ycentroid[i] = tab[i]['ycen_detmap_fw']
            dxcentroid[i] = np.nan
            dycentroid[i] = np.nan
            sig_major_pix[i] = -1 # dummy value
            sig_minor_pix[i] = -1 # dummy value
            ellipticity[i] = -1
            no_centroid_refinement[i] = True
            pos_angle[i] = np.nan
            continue

        cutout = image[(iy_guess-half):(iy_guess+half+1), (ix_guess-half):(ix_guess+half+1)] - med

        gfit = fit_2dgaussian(cutout)

        _xcentroid = gfit.x_mean.value
        _ycentroid = gfit.y_mean.value
        
        ph, pw = cutout.shape
        px, py = np.meshgrid(np.arange(pw), np.arange(ph))

        var_cutout = np.power(sig_adu[(iy_guess-half):(iy_guess+half+1), (ix_guess-half):(ix_guess+half+1)], 2)

        var_xcen = np.sum(var_cutout*np.power(px - _xcentroid, 2))/(np.sum(cutout)**2)*(2.5**2)
        var_ycen = np.sum(var_cutout*np.power(py - _ycentroid, 2))/(np.sum(cutout)**2)*(2.5**2)

        sig_xcen = np.sqrt(var_xcen)
        sig_ycen = np.sqrt(var_ycen)

        if (not np.isfinite(_xcentroid)) or (not np.isfinite(_ycentroid)):
            xcentroid[i] = tab[i]['xcen_detmap_fw']
            ycentroid[i] = tab[i]['ycen_detmap_fw']
            dxcentroid[i] = np.nan
            dycentroid[i] = np.nan
            sig_major_pix[i] = -1
            sig_minor_pix[i] = -1
            ellipticity[i] = -1
            pos_angle[i] = np.nan
            centroid_refinement_fail[i] = True
        else:
            xcentroid[i] = _xcentroid + ix_guess - half
            ycentroid[i] = _ycentroid + iy_guess - half
            # the x_stddev and y_stddev usages here may look confusing,
            # but that's actually how they're defined ...
            sig_major_pix[i] = gfit.x_stddev.value
            sig_minor_pix[i] = gfit.y_stddev.value
            dxcentroid[i] = sig_xcen
            dycentroid[i] = sig_ycen
            if (gfit.x_stddev.value <= 0) or (gfit.y_stddev.value <= 0):
                ellipticity[i] = -1
            else:
                ellipticity[i] = 1.0 - gfit.y_stddev.value/gfit.x_stddev.value
            # gfit pos angle is in radians CC from +X
            # use atan to confine output pos_angle value to [-90, 90] deg
            pos_angle[i] = (180.0/np.pi)*np.arctan(np.tan(gfit.theta.value))

    tab['no_centroid_refinement'] = no_centroid_refinement.astype(int)
    tab['centroid_refinement_fail'] = centroid_refinement_fail.astype(int)
    tab['xcentroid'] = tab['xcen_detmap_fw'] # HACK !!!!!
    tab['ycentroid'] = tab['ycen_detmap_fw'] # HACK !!!!!
    tab['sig_major_pix'] = sig_major_pix
    tab['sig_minor_pix'] = sig_minor_pix
    tab['ellipticity'] = ellipticity
    tab['pos_angle'] = pos_angle
    tab['dxcentroid'] = dxcentroid
    tab['dycentroid'] = dycentroid

def add_metadata_columns(tab, bitmask):
    # input table tab gets modified

    tab['min_edge_dist_pix'] = [util.min_edge_dist_pix(c[0], c[1]) for c in zip(tab['xcentroid'], tab['ycentroid'])]

    tab['dq_flags'] = util.get_dq_flags(tab, bitmask)
    tab['valid_astrom_calibrator'] = ((tab['min_edge_dist_pix'] > 5) & (tab['sig_major_pix'] > 1.0)) # 5 could use more tuning

def get_source_list(image, bitmask, extname, ivar_adu, max_cbox=31,
                    run_aper_phot=True, thresh=5, skip_2dg=False):

    print('Attempting to catalog sources in ' + extname + ' image')

    assert((thresh >= 4) and (thresh <= 100))

    par = common.mask_bit_dict()
    
    image = djs_maskinterp.average_bilinear(image, (np.bitwise_and(bitmask, 1) != 0))

    nominal_fwhm_pix = get_nominal_fwhm_pix(extname)

    detsn = detection_map(image, nominal_fwhm_pix)

    slices = detect_sources(detsn, thresh)

    all_detections = slices_to_table(slices, detsn, extname)

    all_detections = detmap_centroids(all_detections, detsn, max_cbox=max_cbox)

    if len(all_detections) == 0:
        return None, detsn, all_detections, image

    tab = copy.deepcopy(all_detections)

    # only compute this 'sigma map' once to avoid wasted processing time
    sig_adu_map = aper_phot_unc_map(ivar_adu)
    refine_centroids(tab, image, bitmask, ivar_adu, sig_adu=sig_adu_map,
                     skip_2dg=skip_2dg)

    # attempt to remove hot pixels, think this is safe since i end up
    # rejecting (sig_major_pix <= 1) sources when computing
    # overall FWHM and recalibrating the astrometry
    ### tab = tab[tab['sig_major_pix'] > 1.0] # HACK !!!!!
    
    add_metadata_columns(tab, bitmask)

    tab = tab[(tab['min_edge_dist_pix'] > -5)]

    if len(tab) == 0:
        return None, detsn, all_detections, image

    tab['DETMAP_THRESH'] = thresh

    if run_aper_phot:
        do_aper_phot(image, tab, extname, ivar_adu, sig_adu=sig_adu_map)

    # add 'image' to set of outputs since it gets modified
    # and this modification won't persist into the GFA_image object
    # when this is being run via multiprocessing
    return tab, detsn, all_detections, image

def pmgstars_forced_phot(xcentroid, ycentroid, image, elg=False):

    assert(len(xcentroid) > 0)
    assert(len(ycentroid) > 0)

    # create the apertures
    # get the fluxes

    print('Attempting to do forced aperture photometry')

    if elg:
        par = common.gfa_misc_params()

        fname = os.path.join(os.environ[par['meta_env_var']],
                             par['exp_kernel_filename'])

        kern = fits.getdata(fname) # non-optimal to repeatedly read this...

        image = ndimage.convolve(image, kern, mode='constant')

    positions = list(zip(xcentroid, ycentroid))

    radius = 3.567 # pixels

    apertures = CircularAperture(positions, r=radius)
    annulus_apertures = CircularAnnulus(positions, r_in=60.0, r_out=65.0)
    annulus_masks = annulus_apertures.to_mask(method='center')

    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(image)
        annulus_data_1d = annulus_data[mask.data > 0]
        # this sigma_clipped_stats call is actually the slow part !!
        _, median_sigclip, std_bg = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)

    bkg_median = np.array(bkg_median)
    phot = aperture_photometry(image, apertures)

    aper_bkg_tot = bkg_median*_get_area_from_ap(apertures[0])

    aper_fluxes = np.array(phot['aperture_sum']) - aper_bkg_tot

    return aper_fluxes
