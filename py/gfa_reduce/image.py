import gfa_reduce.common as common
import gfa_reduce.imred.dq_mask as dq_mask
import gfa_reduce.analysis.sky as sky
import gfa_reduce.analysis.segment as segment
import gfa_reduce.analysis.phot as phot # may not be needed anymore?
from gfa_reduce.gfa_wcs import nominal_tan_wcs
import numpy as np
import astropy.io.fits as fits
from astropy import wcs
from astropy.stats import mad_std
from scipy.stats import scoreatpercentile
import scipy.ndimage as ndimage
import gfa_reduce.analysis.util as util
from gfa_reduce.analysis.djs_photcen import djs_photcen
import os
from gfa_reduce.analysis.radprof import _atv_radplotf
from gfa_reduce.analysis.splinefwhm import _atv_splinefwhm

class PSF:
    def __init__(self, cube, im_header, cube_index):
        self.psf_image = np.median(cube, 2) # seems to work even for nstars = 1

        self.cube_index = cube_index

        sh = self.psf_image.shape
        assert(sh[0] == sh[1])

        self.sidelen = sh[0]

        # number of radius values for radial profile
        self.nrad = (self.sidelen // 2) + 1 # not clear that this is 100% equivalent to what's in radprof.py

        self.profile_radius_pix = np.zeros(self.nrad, dtype='float32')
        self.radial_profile = np.zeros(self.nrad, dtype='float32')

        self.im_header = im_header # header of the full-frame single-camera
                                   # GFA image

        bgmask = util._stamp_radius_mask(sh[0])
        self.psf_image -= np.median(self.psf_image[bgmask])

        self.psf_image /= np.max(self.psf_image)

        
        self.extname = im_header['EXTNAME']
        self.nstars = cube.shape[2]
        self.cube = cube # maybe get rid of this eventually to save memory

        self.psf_centroiding_flag = 0

        self.cbox = 7
        self.flux_weighted_centroid()

        # don't send the centroids to _aperture_corr_fac since
        # I want the aperture correction factor to have any
        # average off-centering baked in to correct for any such
        # off-centering in the catalog aperture fluxes
        self.aper_corr_fac = util._aperture_corr_fac(self.psf_image)

        self.fiber_fracflux, _, psf_tot_flux = util._fiber_fracflux(self.psf_image,
                                                                    x_centroid=self.xcen_flux_weighted,
                                                                    y_centroid=self.ycen_flux_weighted)

        if self.fiber_fracflux < 0.5:
            self.cbox += (4.0/3.0)*10*(0.5 - max(self.fiber_fracflux, 0))
            assert(self.cbox >= 7)

            self.flux_weighted_centroid() # should i also send the initial djs_photcen (x_start, y_start) here ?

            self.fiber_fracflux, _, psf_tot_flux = util._fiber_fracflux(self.psf_image,
                                                                        x_centroid=self.xcen_flux_weighted,
                                                                        y_centroid=self.ycen_flux_weighted)

        self.elg_convolution()

        _, elg_fiber_flux, __ = util._fiber_fracflux(self.smoothed_psf_image,
                                                     x_centroid=self.xcen_flux_weighted,
                                                     y_centroid=self.ycen_flux_weighted)

        # normalize FIBER_FRACFLUX_ELG based on the total PSF flux *before* convolving
        self.fiber_fracflux_elg = elg_fiber_flux/psf_tot_flux

        if (np.abs(self.xcen_flux_weighted - (self.sidelen // 2)) > 1) or (np.abs(self.ycen_flux_weighted - (self.sidelen // 2)) > 1):
            self.xcen_flux_weighted = float(self.sidelen // 2)
            self.ycen_flux_weighted = float(self.sidelen // 2)
            self.psf_centroiding_flag = 1

        self.fit_moffat_fwhm()

        radii, profile = _atv_radplotf(self.psf_image, self.xcen_flux_weighted, self.ycen_flux_weighted)
        self.radprof_fwhm_asec = _atv_splinefwhm(radii, profile)*0.205

        self.profile_radius_pix = radii
        self.radial_profile = profile

        asymmetry_ratio, asymmetry_numerator, asymmetry_denominator = util._asymmetry_score(self.psf_image,
                                                                                            _xcen=self.xcen_flux_weighted,
                                                                                            _ycen=self.ycen_flux_weighted)
        self.psf_asymmetry_ratio = np.float32(asymmetry_ratio)
        self.psf_asymmetry_numerator = np.float32(asymmetry_numerator)
        self.psf_asymmetry_denominator = np.float32(asymmetry_denominator)

    def psf_image_header(self, hdu):
        
         hdu.header['EXTNAME'] = self.extname
         hdu.header['PETALLOC'] = common.gfa_extname_to_gfa_number(self.extname)
         hdu.header['NSTARS'] = self.nstars
         hdu.header['FIBFRAC'] = self.fiber_fracflux if np.isfinite(self.fiber_fracflux) else 0.0 # ??
         hdu.header['EXPID'] = self.im_header['EXPID']
         hdu.header['CBOX'] = self.cbox
         hdu.header['CFAILED'] = self.psf_centroiding_failed

         if self.cube_index is not None:
             hdu.header['CUBE_IND'] = self.cube_index

         # guide*.fits.fz GUIDE? extensions apparently don't have NIGHT
         # or any other date-related information available
         if 'NIGHT' in self.im_header:
             hdu.header['NIGHT'] = self.im_header['NIGHT']

    def to_hdu(self, primary=False):
         f = (fits.PrimaryHDU if primary else fits.ImageHDU)
         hdu = f(self.psf_image)

         self.psf_image_header(hdu)
         
         return hdu

    def cube_to_hdu(self, primary=False):
         f = (fits.PrimaryHDU if primary else fits.ImageHDU)
         hdu = f(self.cube)

         self.psf_image_header(hdu)

         return hdu
     
    def flux_weighted_centroid(self):
        x_start = y_start = self.sidelen // 2

        print('djs_photcen using cbox = ', self.cbox, ' x_start = ', x_start, ' y_start = ', y_start)
        
        xcen, ycen, q = djs_photcen(x_start, y_start, self.psf_image,
                                    cbox=self.cbox,
                                    cmaxiter=10, cmaxshift=0.0,
                                    ceps=0.0)

        self.xcen_flux_weighted = xcen
        self.ycen_flux_weighted = ycen
        self.psf_centroiding_failed = q

    def fit_moffat_fwhm(self):
        res = util._fit_moffat2d(self.xcen_flux_weighted, self.ycen_flux_weighted, self.psf_image)

        # check for success of minimization ?
        self.moffat_fwhm_pix = res.x[0]

        # could do a more detailed job of this later...
        self.moffat_fwhm_asec = self.moffat_fwhm_pix*0.205

    def elg_convolution(self):
        par = common.gfa_misc_params()
        fname = os.path.join(os.environ[par['meta_env_var']],
                             par['exp_kernel_filename'])

        kern = fits.getdata(fname)

        smth = ndimage.convolve(self.psf_image, kern, mode='constant')

        self.smoothed_psf_image = smth

class Overscan:
    """Object to encapsulate single-camera worth of overscan and prescan"""

    def __init__(self, image):
        # image should be a 2D numpy array with dimensions
        # 2248 x 1032 in the case of DESI GFA cameras

        par = common.gfa_misc_params()

        sh = image.shape
        assert(sh[0] == par['height_with_prescan_overscan'])
        assert(sh[1] == par['width_with_prescan_overscan'])

        amps = common.valid_amps_list()

        self.overscan_cutouts = {}
        self.prescan_cutouts = {}

        for amp in amps:
            bdy = common.overscan_bdy_coords(amp)
            self.overscan_cutouts[amp] = image[bdy['y_l']:bdy['y_u'], bdy['x_l']:bdy['x_u']]
            bdy = common.prescan_bdy_coords(amp)
            self.prescan_cutouts[amp] = image[bdy['y_l']:bdy['y_u'], bdy['x_l']:bdy['x_u']]

        self.n_badpix_overscan = self.count_badpixels()
        self.n_badpix_prescan = self.count_badpixels(prescan=True)

        # still per-amp but summing prescan and overscan counts together
        self.n_badpix = dict([(amp, self.n_badpix_overscan[amp] + self.n_badpix_prescan[amp]) for amp in amps])

        # including all amps and lumping together prescan and overscan
        self.n_badpix_all = np.sum([n for n in self.n_badpix.values()])

        # units are raw ADU
        self.overscan_medians = dict([(amp, np.median(self.overscan_cutouts[amp])) for amp in amps])

        # units are raw ADU
        self.prescan_medians = dict([(amp, np.median(self.prescan_cutouts[amp])) for amp in amps])

    def count_badpixels(self, thresh=10000, prescan=False):
        amps = common.valid_amps_list()
        
        if prescan:
            return dict([(amp, int(np.sum(self.prescan_cutouts[amp] > thresh))) for amp in amps])
        else:
            return dict([(amp, int(np.sum(self.overscan_cutouts[amp] > thresh))) for amp in amps])

    def bad_amps_list(self):
        bad_amps = []
        for k, v in self.n_badpix.items():
            if v >= 10:
                bad_amps.append(k)

        return bad_amps

class GFA_image:
    """Single GFA image from one GFA exposure"""

    def __init__(self, image, header, cube_index=None, store_detmap=False, coadd_index_range=None):
        if cube_index is None:
            self.image = image.astype('float32')
            self.nframe = 1
        elif cube_index == -1:
            cube_index_start = coadd_index_range[0]
            cube_index_end = coadd_index_range[1]
            self.image = np.mean(image[cube_index_start:(cube_index_end + 1), :, :].astype('float32'),
                                 axis=0)
            self.nframe = image.shape[0]
        else:
            self.image = image[cube_index, :, :].astype('float32')
            self.nframe = 1

        self.coadd_index_range = coadd_index_range

        self.overscan = Overscan(self.image)
        self.remove_overscan()

        self.store_detmap = store_detmap
        self.detmap = None
        self.full_detlist = None
        
        par = common.gfa_misc_params()
        # important that this saturation threshold be done on truly raw image..
        self.satmask = (self.image > par['sat_thresh_adu']).astype('byte')
        
        self.cube_index = cube_index
        self.header = header
        self.extname = self.header['EXTNAME'].replace(' ', '')
        self.header['CONTRAST'] = 0.0 # typically overwritten with actual value

        # for sequences that include spectro and guiding components
        # REQTIME in guide cube headers can sometimes be the spectro
        # REQTIME (see examples from early March 2020), so I want to avoid
        # getting the wrong time for dark current scaling and transparency
        # estimation by making it impossible to grab REQTIME from the header
        # for guide cubes (think that REQTIME in the guide binary tables
        # should be fine though, could check more exhaustively)
        if (cube_index is not None) and ('REQTIME' in self.header):
            del self.header['REQTIME']

        # similarly, for guide cubes, i want to get MJD-OBS of each frame
        # and not for the start of the guide sequence, so i want to
        # make it impossible to accidentally grab MJD-OBS from a guide cube
        # image extension header (this was not a problem until
        # 20200314 when MJD-OBS showed up in these guide cube image
        # extensions, whereas it had previously not been present)
        if (cube_index is not None) and ('MJD-OBS' in self.header):
            del self.header['MJD-OBS']

        self.initialize_wcs()

        # lazily compute bitmask image only as requested
        self.bitmask = None

        # may want to be more clever about this to allow for the possibility
        # of loading in already reduced images in the future
        self.bias_subtracted = False
        self.dark_subtracted = False
        self.flatfielded = False

        self.var_e_sq = None
        # the _adu in ivar_adu is meant to indicate the units are 1/(ADU^2)
        self.ivar_adu = None
        self.sky_mag = None
        self.sky_mag_per_amp = None
        self.segmap = None # may want to get rid of this entirely
        self.empirical_bg_sigma = None
        self.sky_level_adu = None
        self.sky_level_adu_per_amp = None
        self.bintable_row = None

        # record CCD temperature used for dark current removal
        self.t_c_for_dark = None

        # should become a boolean indicating whether CCD temperature
        # used for dark removal was a guess (1) or is thought to
        # be correct (0)
        self.t_c_for_dark_is_guess = None
        
        # record exposure time used for dark current removal
        self.time_s_for_dark = None

        self.psf = None

        self.max_cbox = 31
        
    def create_dq_mask(self, dark_image):
        if self.bitmask is not None:
            return

        d = common.mask_bit_dict()
        thresh = scoreatpercentile(dark_image, 99.5)
        if self.bitmask is None:
            self.bitmask = ((dark_image > thresh)*(2**d['HOTDARK'])).astype('byte')
        else:
            self.bitmask += ((dark_image > thresh)*(2**d['HOTDARK'])).astype('byte')
        self.bitmask += self.satmask*(2**d['SATUR'])

        self.bitmask = self.bitmask.astype('byte') # just to make sure...
        del self.satmask

    def update_bitmask_flat(self, flatfield):
        # doing this to avoid having to keep flatfield images in memory

        thresh = 0.6 # very little thought put into this choice..
        d = common.mask_bit_dict()
        if self.bitmask is None:
            self.bitmask = ((flatfield < thresh)*(2**d['FLATBAD'])).astype('byte')
        else:
            self.bitmask += ((flatfield < thresh)*(2**d['FLATBAD'])).astype('byte')
        
    def calc_variance_e_squared(self):
        # at this stage the image ought to have been bias subtracted
        # but not flatfielded or dark subtracted

        assert(self.bias_subtracted)

        assert((not self.dark_subtracted) and (not self.flatfielded))

        gain = common.gfa_camera_gain(self.extname)

        var_e_sq = (common.gfa_camera_readnoise(self.extname)**2 + \
                    self.image*(self.image >= 0)*gain)
        var_e_sq /= self.nframe

        assert(np.sum(var_e_sq <= 0) == 0)
        self.var_e_sq = var_e_sq

    def calc_variance_adu(self, flatfield=None):

        gain = common.gfa_camera_gain(self.extname)

        variance_adu_sq = self.var_e_sq/(gain**2)

        del self.var_e_sq
        self.var_e_sq = None

        if flatfield is not None:
            variance_adu_sq *= (flatfield**2)

        assert(np.sum(variance_adu_sq <= 0) == 0)

        # note that I'm not currently taking into account uncertainty due
        # to the uncertainty on the flatfield; not sure if doing so would
        # be desirable eventually

        ivar_adu = 1.0/variance_adu_sq

        assert(self.bitmask is not None)

        # zero out inverse variance of pixels with anything flagged in
        # data quality bitmask
        ivar_adu *= (self.bitmask == 0)

        self.ivar_adu = ivar_adu

    def to_hdu(self, primary=False, flavor=''):
        # convert this image to an HDU

        # currently expect flavor to be one of
        #     REDUCED - reduced image
        #     BITMASK - data quality bitmask
        #     INVVAR - inverse variance image
        #     DETMAP - detection significance map

        # if no flavor is specified then assume ".image" attribute is desired
        # data for this HDU

        f = (fits.PrimaryHDU if primary else fits.ImageHDU)

        if (flavor == '') or (flavor == 'REDUCED'):
            hdu = f(self.image.astype('float32'), header=self.header)
        elif (flavor == 'BITMASK'):
            hdu = f(self.bitmask.astype('int'), header=self.header)
            hdu.header = dq_mask.add_dq_bitmask_header_cards(hdu.header)
        elif (flavor == 'INVVAR'):
            hdu = f(self.ivar_adu.astype('float32'), header=self.header)
        elif (flavor == 'DETMAP'):
            assert(self.detmap is not None)
            hdu = f(self.detmap.astype('float32'), header=self.header)

        hdu.header['FLAVOR'] = flavor

        gain = common.gfa_camera_gain(self.extname)
        hdu.header['GAINA'] = (gain, '[e-/ADU] assumed gain')

        hdu.header['BUNIT'] = common.reduced_flavor_to_bunit(flavor)

        petal_loc = common.gfa_extname_to_gfa_number(self.extname)
        hdu.header['PETALLOC'] = (petal_loc, 'petal number')

        return hdu

    def are_pixels_calibrated(self, flatfielding_on=True):
        if flatfielding_on:
            result = (self.bias_subtracted and self.dark_subtracted and 
                      self.flatfielded)
        else:
            result = (self.bias_subtracted and self.dark_subtracted)

        return result

    def estimate_sky_level(self, careful_sky=False, flatfielding_on=True):
        # do something dumb for now, return to this later with something
        # more sophisticated, possibly involving segmentation
        # and/or finding the mode

        assert(self.are_pixels_calibrated(flatfielding_on=flatfielding_on))

        if careful_sky:
            if self.segmap is None:
                self.set_segmap()

            self.sky_level_adu = np.median(self.image[self.segmap.array == 0])
        else:
            self.sky_level_adu = np.median(self.image)
            self.sky_level_adu_per_amp = {}
            for amp in common.valid_amps_list():
                bdy = common.amp_bdy_coords(amp)
                self.sky_level_adu_per_amp[amp] = np.median(self.image[bdy['y_l']:bdy['y_u'], bdy['x_l']:bdy['x_u']])
                

        return self.sky_level_adu

    def compute_segmap(self):
        print('Attempting to compute segmentation map for ' + self.extname)

        segmap = segment.segmentation_map(self.image, self.extname)
        return segmap

    def set_segmap(self):
        self.segmap = self.compute_segmap()

    def compute_empirical_bg_sigma(self, careful_sky=False):

        if careful_sky:
            if self.segmap is None:
                self.set_segmap()

            # this could go wrong in pathological case that
            # segmap is nonzero for all pixels
            return mad_std(self.image[self.segmap.array == 0])
        else:
            return mad_std(self.image)

    def set_empirical_bg_sigma(self, careful_sky=False):
        if self.empirical_bg_sigma is None:
            self.empirical_bg_sigma = self.compute_empirical_bg_sigma(careful_sky=careful_sky)

    def estimate_sky_mag(self, careful_sky=False, flatfielding_on=True):
        # calculate sky brightness in mag per sq asec
        # this is meant to be run on the reduced image in ADU

        assert(self.are_pixels_calibrated(flatfielding_on=flatfielding_on))

        sky_adu_per_pixel = self.estimate_sky_level(careful_sky=careful_sky,
                                                    flatfielding_on=flatfielding_on)

        acttime = self.time_s_for_dark
        
        sky_mag = sky.adu_to_surface_brightness(sky_adu_per_pixel, 
                                                acttime, self.extname)

        if self.sky_level_adu_per_amp is not None:
            amps = common.valid_amps_list()
            self.sky_mag_per_amp = [sky.adu_to_surface_brightness(self.sky_level_adu_per_amp[amp], acttime, self.extname) for amp in amps]

        print(self.extname + ' sky mag per square asec AB : ' +  
              '{:.3f}'.format(sky_mag))

        self.sky_mag = sky_mag

        # calculate sky mag for just upper 8 rows (lumping together amps G and H)
        self.sky_mag_upper()

        return sky_mag

    def sky_mag_upper(self):
        # only upper 8 rows of image -- see DESI-5334
        self.sky_level_adu_upper = np.median(self.image[1024:1032, 0:2047])
        self.sky_mag_upper = sky.adu_to_surface_brightness(self.sky_level_adu_upper, 
                                                           self.time_s_for_dark, self.extname)

    def catalog_add_radec(self, catalog):
        # use wcs attribute to convert from pixel to world coordinates
        # be careful about 1-indexed versus 0-indexed convention
        # also be careful about any swapping of x and y pixel coordinates
        
        ra, dec = self.wcs.all_pix2world(catalog['xcentroid'],
                                         catalog['ycentroid'], 0)

        catalog['ra'] = ra
        catalog['dec'] = dec

        return catalog

    def ingest_cataloging_results(self, tab, detmap, alldet, image):

        # image gets modified slightly during source detection stage
        # specifically via djs_maskinterp attempting to interpolate
        # over bad pixels
        self.image = image

        # always store alldet since it shouldn't be consuming any
        # appreciable amount of memory
        if len(alldet) > 0:
            alldet['extname'] = self.extname
        self.full_detlist = alldet # should be an astropy Table

        if self.store_detmap:
            self.detmap = detmap

        del detmap
        
        n_sources = (len(tab) if tab is not None else 0)

        print('Found ' + str(n_sources) + ' sources in ' +
              self.extname + ' image')

        if tab is None:
            return tab

        tab = self.catalog_add_radec(tab)

        mjd_obs = self.try_retrieve_meta_keyword('MJD-OBS')
        if mjd_obs is None:
            print('could not find MJD-OBS header keyword !!!')
        else:
            tab['mjd_obs'] = mjd_obs

        util.add_ampname_to_catalog(tab)
        util.sanity_check_catalog(tab)
        return tab

    def ingest_dark_current_results(self, dark_image):
        self.image = self.image - dark_image
        self.dark_subtracted = True
        self.create_dq_mask(dark_image)

    def initialize_wcs(self):
        telra = self.header['SKYRA']
        teldec = self.header['SKYDEC']

        print('Attempting to initialize WCS guess for ' + self.extname)
        self.wcs = nominal_tan_wcs(telra, teldec, self.extname)

    def remove_overscan(self):
        sh = self.image.shape
        if sh[1] == 2248:
            _image = np.zeros((1032, 2048), dtype=float)
            _image[:, 0:1024] = self.image[:, 50:1074]
            _image[:, 1024:2048] = self.image[:, 1174:2198]
            self.image = _image

    def update_wcs(self, d):
        # d is a dictionary with xshift_best, yshift_best
        # for this EXTNAME

        assert(d['extname'] == self.extname)

        self.wcs.wcs.crpix = self.wcs.wcs.crpix + np.array([d['xshift_best'], d['yshift_best']])

        # also want to update the header

        new_wcs_header_cards = self.wcs.to_header()
        new_wcs_header_cards['CONTRAST'] = d['contrast']

        for k,v in new_wcs_header_cards.items():
            if k == 'LATPOLE':
                continue
            self.header[k] = v

        self.header['CD1_1'] = self.header['PC1_1']
        self.header['CD2_1'] = self.header['PC2_1']
        self.header['CD1_2'] = self.header['PC1_2']
        self.header['CD2_2'] = self.header['PC2_2']

        del self.header['PC1_1']
        del self.header['PC2_1']
        del self.header['PC1_2']
        del self.header['PC2_2']

    def try_retrieve_meta_keyword(self, keyword, placeholder=None):
        # examples are MJD-OBS and GCCDTEMP, which
        # are found in different places in the raw data depending
        # on gfa*.fits.fz versus guide*.fits.fz

        # because guider cube metadata has evolved over time, won't
        # always be guaranteed to get e.g., GCCDTEMP at all
        
        # first look in the image header (could be dangerous for EXPTIME
        # in the case of guider cubes)

        if self.bintable_row is not None:
            if self.cube_index != -1:
                bintable_has_keyword = keyword in self.bintable_row.array.dtype.names
            else:
                bintable_has_keyword = keyword in self.bintable_row.colnames

        if keyword in self.header.keys():
            return self.header[keyword]
        elif (self.bintable_row is not None) and bintable_has_keyword:
            return self.bintable_row[keyword]
        else:
            print('could not find ' + keyword + ' !!')
            return placeholder
            
        
    def extract_psf_cutouts(self, __catalog, sidelen=51):

        # handle case of entire exposure with no retained sources
        if __catalog is None:
            return
        
        # sidelen should be an integer...
        assert(np.round(sidelen) == sidelen)

        half = sidelen // 2
        bad_amps = self.overscan.bad_amps_list()

        _catalog = __catalog[__catalog['camera'] == self.extname]
        if len(_catalog) == 0:
            return None
        
        keep = util.use_for_fwhm_meas(_catalog, bad_amps=bad_amps, no_sig_major_cut=True) & (_catalog['min_edge_dist_pix'] > (half + 0.5)) & (_catalog['aper_sum_bkgsub_3'] > 0)
        
        if np.sum(keep) == 0:
            return None

        n = np.sum(keep)
        #cube = np.zeros((sidelen, sidelen, n))

        catalog = _catalog[keep]

        assert(np.sum(catalog['aper_sum_bkgsub_3'] <= 0) == 0)

        cutouts = []
        for i in range(n):
            ixcentroid = int(np.round(catalog['xcentroid'][i]))
            iycentroid = int(np.round(catalog['ycentroid'][i]))
            cutout = self.image[(iycentroid-half):(iycentroid+half+1),
                                (ixcentroid-half):(ixcentroid+half+1)]

            # hack to try removing saturated sources
            if np.sum(cutout >= 30000.0) > 1:
                continue

            dx = np.round(catalog['xcentroid'][i]) - catalog['xcentroid'][i]
            dy = np.round(catalog['ycentroid'][i]) - catalog['ycentroid'][i]

            cutout = util._shift_stamp(cutout, dx, dy)

            # background subtract
            bgmask = util._stamp_radius_mask(sidelen)

            bg = np.median(cutout[bgmask])

            cutout -= bg
            
            cutout = cutout/catalog['aper_sum_bkgsub_3'][i]

            cutouts.append(cutout)

        ncutouts = len(cutouts)
        
        if ncutouts == 0:
            return None

        cube = np.zeros((sidelen, sidelen, ncutouts))
        for i, cutout in enumerate(cutouts):
            cube[:, :, i] = cutout
        
        return cube

    def create_psf(self, catalog, sidelen=51):
        cube = self.extract_psf_cutouts(catalog, sidelen=sidelen)
        print('computing PSF for ' + self.extname)
        if cube is None:
            self.psf = None
            print("WARNING: did not find any PSF 'stars' for " + self.extname)
        else:
            self.psf = PSF(cube, self.header, self.cube_index)

    def compute_zeropoint(self, ps1_matched_catalog):
        if ps1_matched_catalog is None:
            return np.nan

        if self.psf is None:
            return np.nan

        # require ps1 median_1_ r flux to be > 0
        # require correct extname
        # detmap_peak >= 10 (kind of like S/N > 10)
        # ang_sep_deg < 2.0/3600.0

        # require GFA _3 flux > 0
        # require something about minimum edge distance
        # cut on dq_flags
        
        good = (ps1_matched_catalog['median_1_'] > 0) & (ps1_matched_catalog['camera'] == self.extname) & (ps1_matched_catalog['detmap_peak'] >= 10) & (ps1_matched_catalog['ang_sep_deg'] < 2.0/3600.0) & (ps1_matched_catalog['aper_sum_bkgsub_3'] > 0) & (ps1_matched_catalog['min_edge_dist_pix'] >= 10) & (ps1_matched_catalog['dq_flags'] == 0)

        if np.sum(good) == 0:
            return np.nan

        ps1_matched_catalog = ps1_matched_catalog[good]
        r_ps1 = -2.5*np.log10(ps1_matched_catalog['median_1_'])
        m_inst = -2.5*np.log10(ps1_matched_catalog['aper_sum_bkgsub_3']/(self.time_s_for_dark*self.psf.aper_corr_fac))

        zp = np.median(r_ps1 - m_inst)

        # would be good to return metrics
        # regarding how many/which sources were used for determinining the zeropoint
        # what their mag range was, maybe even how well the slope of m_inst vs m_ps1 matches with unity
        
        return zp
