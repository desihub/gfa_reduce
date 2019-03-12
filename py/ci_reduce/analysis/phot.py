import ci_reduce.analysis.util as util
from astropy.stats import mad_std
from scipy.ndimage.filters import gaussian_filter
import numpy as np
from scipy.ndimage.measurements import label, find_objects
from astropy.table import Table
from photutils.centroids.core import fit_2dgaussian
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import CircularAperture, CircularAnnulus

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

    tab = Table([xcen, ycen, xmin, xmax, ymin, ymax, detmap_peak], names=('xcen_init', 'ycen_init', 'detmap_xmin', 'detmap_xmax', 'detmap_ymin', 'detmap_ymax', 'detmap_peak'))

    return tab

def detection_map(im, fwhm):
    psf_sigma = fwhm/2.355

    smth = gaussian_filter(im, psf_sigma)

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

def do_aper_phot(data, catalog, extname, ivar_adu):
    # catalog should be the catalog with refined centroids 
    # for **one CI camera**

    print('Attempting to do aperture photometry')
    positions = list(zip(catalog['xcentroid'], catalog['ycentroid']))
    apertures = CircularAperture(positions, r=get_nominal_fwhm_pix(extname)/2.0)
    annulus_apertures = CircularAnnulus(positions, r_in=30.0, r_out=40.0)
    annulus_masks = annulus_apertures.to_mask(method='center')

    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        # this sigma_clipped_stats call is actually the slow part !!
        _, median_sigclip, std_bg = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)

    bkg_median = np.array(bkg_median)
    phot = aperture_photometry(data, apertures, 
                               error=aper_phot_unc_map(ivar_adu))
    phot['annulus_median'] = bkg_median
    phot['aper_bkg'] = bkg_median * apertures.area()
    phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']

    # still need to actually append some of this aperture photometry
    # info to the catalog in the form of new columns

    catalog['aper_sum_bkgsub'] = phot['aper_sum_bkgsub']
    catalog['aper_bkg'] = phot['aper_bkg']
    catalog['aperture_sum_err'] = phot['aperture_sum_err']
    # is .area() result a vector or scalar ??
    catalog['annulus_area_pix'] = annulus_apertures.area()

def get_nominal_fwhm_pix(extname):
    # this is a nominal FWHM for use as an initial guess
    # when creating the detection map, when nothing is yet known
    # about the FWHM of the particular data being analyzed
    # could revisit this and/or clean it up later on
    nominal_fwhm_pix = (9.3766944 if (extname == 'CIC') else 10.192287)

    return nominal_fwhm_pix

def refine_centroids(tab, image, bitmask):
    # input table tab gets augmented with additional columns

    print('Attempting to refine initial centroids')
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

    no_centroid_refinement = np.zeros(nobj, dtype=bool)
    centroid_refinement_fail = np.zeros(nobj, dtype=bool)

    for i in range(nobj):
        ix_guess = int(round(tab[i]['xcen_init']))
        iy_guess = int(round(tab[i]['ycen_init']))
        min_edge_dist = util.min_edge_dist_pix(ix_guess, iy_guess)

        # if the centering box extends outside of image, don't
        # try to do any refinement of the centroid
        if min_edge_dist < half:
            xcentroid[i] = tab[i]['xcen_init']
            ycentroid[i] = tab[i]['ycen_init']
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

        if (not np.isfinite(_xcentroid)) or (not np.isfinite(_ycentroid)):
            xcentroid[i] = tab[i]['xcen_init']
            ycentroid[i] = tab[i]['ycen_init']
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
            if (gfit.x_stddev.value <= 0) or (gfit.y_stddev.value <= 0):
                ellipticity[i] = -1
            else:
                ellipticity[i] = 1.0 - gfit.y_stddev.value/gfit.x_stddev.value
            # gfit pos angle is in radians CC from +X
            # use atan to confine output pos_angle value to [-90, 90] deg
            pos_angle[i] = (180.0/np.pi)*np.arctan(np.tan(gfit.theta.value))

    tab['no_centroid_refinement'] = no_centroid_refinement.astype(int)
    tab['centroid_refinement_fail'] = centroid_refinement_fail.astype(int)
    tab['xcentroid'] = xcentroid
    tab['ycentroid'] = ycentroid
    tab['sig_major_pix'] = sig_major_pix
    tab['sig_minor_pix'] = sig_minor_pix
    tab['ellipticity'] = ellipticity
    tab['pos_angle'] = pos_angle

def add_metadata_columns(tab, bitmask):
    # input table tab gets modified

    min_edge_dist = [util.min_edge_dist_pix(c[0], c[1]) for c in zip(tab['xcentroid'], tab['ycentroid'])]
    tab['min_edge_dist_pix'] = min_edge_dist

    # there's probably a scipy function for this but w/e

    xmin = util.ci_pixel_xmin(pix_center=True)
    xmax = util.ci_pixel_xmax(pix_center=True)
    ymin = util.ci_pixel_ymin(pix_center=True)
    ymax = util.ci_pixel_ymax(pix_center=True)

    ixs = [int(min(max(np.round(t['xcentroid']), xmin), xmax)) for t in tab]
    iys = [int(min(max(np.round(t['ycentroid']), ymin), ymax)) for t in tab]

    tab['dq_flags'] = bitmask[iys, ixs]

def get_source_list(image, bitmask, extname, ivar_adu, thresh=5):
    filler_value = np.median(image)

    # should do something like djs_maskinterp instead
    image[bitmask != 0] = filler_value
 
    nominal_fwhm_pix = get_nominal_fwhm_pix(extname)

    detsn = detection_map(image, nominal_fwhm_pix)

    slices = detect_sources(detsn, thresh)

    tab = slices_to_table(slices, detsn, extname)

    refine_centroids(tab, image, bitmask)

    do_aper_phot(image, tab, extname, ivar_adu)

    add_metadata_columns(tab, bitmask)

    return tab
