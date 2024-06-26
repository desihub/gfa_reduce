# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
gfa_reduce.exposure
===================

Class that encapsulates a single GFA exposure consisting of multiple single-camera images.
"""
import gfa_reduce.common as common
import gfa_reduce.imred.load_calibs as load_calibs
import gfa_reduce.dark_current as dark_current
import astropy.io.fits as fits
import numpy as np
import gfa_reduce.analysis.util as util
import gfa_reduce.analysis.phot as phot
from multiprocessing import Pool

from desiutil.log import get_logger


class GFA_exposure:
    """Object encapsulating the contents of a single GFA exposure"""

    def __init__(self, image_list, exp_header=None, bintables=None,
                 max_cbox=31, pmgstars=None):
        # images is a dictionary of GFA_image objects
        self.log = get_logger()
        par = common.gfa_misc_params()
        _extnames = [_im.header['EXTNAME'] for _im in image_list]
        _extnames.sort()
        self.images = dict(zip(_extnames,
                               par['n_cameras']*[None]))

        self.dark_current_objs = dict(zip(common.valid_image_extname_list(),
                                          par['n_cameras']*[None]))

        self.assign_image_list(image_list)

        # exposure-level header
        self.exp_header = exp_header

        # hack for 20210106
        if self.exp_header is not None:
            if self.exp_header['SKYRA'] is None:
                self.log.info('REPLACING GUIDER SKYRA WITH REQRA')
                self.exp_header['SKYRA'] = self.exp_header['REQRA']
            if self.exp_header['SKYDEC'] is None:
                self.log.info('REPLACING GUIDER SKYDEC WITH REQDEC')
                self.exp_header['SKYDEC'] = self.exp_header['REQDEC']

        self.pixels_calibrated = None
        self.bintables = bintables
        self.max_cbox = max_cbox
        self.assign_max_cbox() # to the per-camera images ...
        self.pmgstars = pmgstars

        # eventually may get assigned to be the _ccds summary table
        # object for this exposure
        self.ccds = None

    def assign_one_image(self, image):
        extname = (image.header)['EXTNAME']
        self.images[extname] = image

    def assign_image_list(self, image_list):
        for image in image_list:
            self.assign_one_image(image)

    def assign_max_cbox(self):
        for image in self.images.values():
            if image is not None:
                image.max_cbox = self.max_cbox

    def subtract_bias(self):
        self.log.info('Attempting to subtract bias...')
        for extname in self.images.keys():
            if self.images[extname] is not None:
                self.images[extname].image = self.images[extname].image - \
                    load_calibs.read_bias_image(extname)

                amps = common.valid_amps_list()
                for amp in amps:
                    bdy = common.amp_bdy_coords(amp)
                    self.images[extname].image[bdy['y_l']:bdy['y_u'],
                                               bdy['x_l']:bdy['x_u']] -= \
                        self.images[extname].overscan.overscan_medians[amp]

                self.images[extname].bias_subtracted = True
                self.images[extname].calc_variance_e_squared()

    def subtract_dark_current(self, do_dark_rescaling=True, mp=False):
        self.log.info('Attempting to subtract dark current...')

        if mp:
            args = []

        for extname in self.images.keys():
            if self.images[extname] is not None:

                acttime = self.images[extname].try_retrieve_meta_keyword('EXPTIME')
                if (acttime is None) or (acttime == 0):
                    self.log.info('trying REQTIME instead of EXPTIME')
                    acttime = self.images[extname].try_retrieve_meta_keyword('REQTIME')
                if (acttime is None) or (acttime == 0):
                    self.log.critical('could not find an exposure time !!!!')
                    assert(False) # die

                self.images[extname].time_s_for_dark = acttime

                self.images[extname].t_c_for_dark_is_guess = False
                t_c = self.images[extname].try_retrieve_meta_keyword('GCCDTEMP')
                if t_c is None:
                    self.log.info('trying GCOLDTEC instead of GCCDTEMP')
                    t_c = self.images[extname].try_retrieve_meta_keyword('GCOLDTEC')

                if t_c is None:
                    self.log.warning('could not find a CCD temperature !!!!')
                    self.images[extname].t_c_for_dark_is_guess = True
                    t_c = 11.0 # HACK !!!!

                self.images[extname].t_c_for_dark = t_c


                if not mp:
                    dark_image, dc_obj = dark_current.total_dark_image_adu(extname,
                                                                           acttime, t_c,
                                                                           self.images[extname].image,
                                                                           do_dark_rescaling=do_dark_rescaling)

                    self.images[extname].ingest_dark_current_results(dark_image)
                    self.dark_current_objs[extname] = dc_obj
                else:
                    args.append((extname, acttime, t_c, self.images[extname].image, do_dark_rescaling))

        if mp:
            self.log.info('Computing dark scalings for all guide cameras in parallel...')
            nproc = len(args)
            assert(nproc <= 6)
            p = Pool(nproc)
            results = p.starmap(dark_current.total_dark_image_adu, args)

            for i, result in enumerate(results):
                extname = args[i][0]
                self.images[extname].ingest_dark_current_results(result[0])
                self.dark_current_objs[extname] = result[1]

            del args, results

    def apply_flatfield(self):
        self.log.info('Attempting to apply flat field...')
        for extname in self.images.keys():
            if self.images[extname] is not None:
                flatfield = load_calibs.read_flat_image(extname)
                self.images[extname].image = self.images[extname].image / \
                    flatfield
                self.images[extname].flatfielded = True
                self.images[extname].calc_variance_adu(flatfield=flatfield)
                self.images[extname].update_bitmask_flat(flatfield)

    def calibrate_pixels(self, do_dark_rescaling=True, mp=False, do_apply_flatfield=True):
        self.subtract_bias()
        self.subtract_dark_current(do_dark_rescaling=do_dark_rescaling, mp=mp)
        if do_apply_flatfield:
            self.apply_flatfield()
        else:
            self.log.info('Skipping flatfielding')
            for extname in self.images.keys():
                self.images[extname].calc_variance_adu()
        self.pixels_calibrated = True

    def num_images_populated(self):
        return sum( im != None for im in self.images.values() )

    def populated_extnames(self):
        return [k for k,v in self.images.items() if v is not None]

    def to_hdulist(self, flavor=''):
        extname_list = common.valid_image_extname_list()
        hdulist = []
        for extname in extname_list:
            if self.images[extname] is not None:
                hdu = self.images[extname].to_hdu(primary=(len(hdulist) == 0),
                                                  flavor=flavor)
                hdulist.append(hdu)

        return fits.HDUList(hdulist)

    def estimate_all_sky_mags(self, careful_sky=False, flatfield_applied=True):
        for im in self.images.values():
            if im is not None:
                im.estimate_sky_mag(careful_sky=careful_sky,
                                    flatfielding_on=flatfield_applied)

    def estimate_all_sky_sigmas(self, careful_sky=False):
        for im in self.images.values():
            if im is not None:
                im.set_empirical_bg_sigma(careful_sky=careful_sky)
                self.log.info('empirical ' + im.header['EXTNAME'] +
                              ' background sigma = ' +
                              '{:.2f}'.format(im.empirical_bg_sigma) + ' ADU')

    def all_source_catalogs(self, mp=False, run_aper_phot=True,
                            det_sn_thresh=5, skip_2dg=False):
        tables = dict(zip(self.images.keys(), len(self.images.keys())*[None]))

        if not mp:
            for extname, im in self.images.items():
                if im is not None:
                    # tab is a culled and augmented list of sources including e.g.,
                    # refined centroids and photometry

                    # alldet is just the initial, raw list of all detections with
                    # no culling applied

                    tab, detmap, alldet, image = phot.get_source_list(im.image,
                                                                      im.bitmask,
                                                                      im.extname,
                                                                      im.ivar_adu,
                                                                      max_cbox=im.max_cbox,
                                                                      run_aper_phot=run_aper_phot,
                                                                      thresh=det_sn_thresh,
                                                                      skip_2dg=skip_2dg)

                    tables[extname] = im.ingest_cataloging_results(tab, detmap,
                                                                   alldet, image)

        else:
            args = []
            for extname, im in self.images.items():
                if im is not None:
                    args.append((im.image, im.bitmask, im.extname, im.ivar_adu, im.max_cbox, run_aper_phot, det_sn_thresh, skip_2dg))

            self.log.info('Running source cataloging for all guide cameras in parallel...')
            nproc = len(args)
            assert(nproc <= 6)
            p = Pool(nproc)
            results = p.starmap(phot.get_source_list, args)

            for i, result in enumerate(results):
                extname = args[i][2] # hopefully this indexing doesn't change...
                tables[extname] = self.images[extname].ingest_cataloging_results(*result)

        # do I also want to store the tables as an attribute belonging to
        # this exposure object?
        return tables

    def update_wcs(self, astr):

        # don't crash for case when no sources were retained in entire
        # exposure
        if astr is None:
            return

        for a in astr:
            extname = a['extname']
            self.images[extname].update_wcs(a)

    def recompute_catalog_radec(self, cat):

        # don't crash for case when no sources were retained in entire
        # exposure
        if cat is None:
            return

        extnames = np.unique(cat['camera'])

        for extname in extnames:
            cat[cat['camera'] == extname] = self.images[extname].catalog_add_radec(cat[cat['camera'] == extname])

    def set_bintable_rows(self):
        for image in self.images.values():
            if image is None:
                continue

            extname = image.header['EXTNAME'].strip()
            if image.cube_index is None:
                image.bintable_row = None
            elif image.cube_index == -1:
                image.bintable_row = util.average_bintable_metadata(self.bintables[extname][image.coadd_index_range[0]:(image.coadd_index_range[1]+1)])
            else:
                image.bintable_row = self.bintables[extname][image.cube_index]

    def purge_zero_exptime(self):
        # last frame of a DESI sequence guider cube can apparently end
        # up with exactly zero EXPTIME
        # try to remove such cases from downstream processing
        # an example is NIGHT = 20210218, EXPID = 76803, CUBE_INDEX = 137, EXTNAME = GUIDE3

        extnames = list(self.images.keys())
        for extname in extnames:
            if self.images[extname].bintable_row is not None:
                if self.images[extname].cube_index == -1:
                    has_exptime_column = 'EXPTIME' in self.images[extname].bintable_row.colnames
                else:
                    has_exptime_column = 'EXPTIME' in self.images[extname].bintable_row.array.columns.names

                if has_exptime_column:
                    if self.images[extname].bintable_row['EXPTIME'] == 0:
                        self.log.warning('DISCARDING ' + extname + ' BECAUSE IT HAS ZERO EXPOSURE TIME')
                        del self.images[extname]
                        if self.pmgstars is not None:
                            self.pmgstars = self.pmgstars[np.where(self.pmgstars['GFA_LOC'] != extname)[0]]

        if len(self.images) == 0:
            self.log.critical('ALL CAMERAS HAVE ZERO EXPOSURE TIME??')
            assert(False)

    def try_retrieve_header_card(self, keyword, placeholder=None):
        if (keyword in self.exp_header.keys()) and (self.exp_header[keyword] is not None):
            return self.exp_header[keyword]
        else:
            self.log.error('could not find ' + keyword + ' in exposure-level header !!')
            return placeholder

    def compute_psfs(self, catalog):
        for image in self.images.values():
            if image is None:
                continue
            image.create_psf(catalog)

    def assign_input_filename(self, fname_in):
        self.fname_in = fname_in

    def _trim_pmgstars(self):
        # remove PMGSTARS table rows that pertain to a camera
        # not actually present in the guide cube
        # example is GUIDE3 in EXPID = 83995 on night 20210408

        keep = np.zeros(len(self.pmgstars), dtype=bool)
        for extname in self.images.keys():
            keep[self.pmgstars['GFA_LOC'] == extname] = True

        self.pmgstars = self.pmgstars[keep]

    def pmgstars_dq_flags(self):
        # gather dq_flags for PMGSTARS table source positions

        extnames = self.images.keys()

        dq_flags = np.zeros(len(self.pmgstars), dtype='uint8')

        for extname in extnames:
            mask = (self.pmgstars['GFA_LOC'] == extname)
            dq_flags[mask] = util.get_dq_flags(self.pmgstars[mask],
                                               self.images[extname].bitmask)

        self.pmgstars['dq_flags'] = dq_flags

    def pmgstars_airmass_exptime(self):
        # add column for airmass to PMGSTARS table, where
        # there is a unique airmass value per GFA camera

        airmass_per_gfa = np.full(len(self.pmgstars), np.nan)
        time_s_for_dark = np.full(len(self.pmgstars), np.nan)

        for row in self.ccds:
            mask = (self.pmgstars['GFA_LOC'] == row['extname'])
            if np.sum(mask) == 0:
                continue
            airmass_per_gfa[mask] = row['airmass_per_gfa']
            time_s_for_dark[mask] = row['time_s_for_dark']

        self.pmgstars['airmass_per_gfa'] = airmass_per_gfa
        self.pmgstars['time_s_for_dark'] = time_s_for_dark

        assert(np.sum(np.isnan(self.pmgstars['airmass_per_gfa'])) == 0)

    def pmgstars_zp_clear(self):
        # for each row of the PMGSTARS table, compute the
        # zeropoint (1 ADU/second) under clear conditions for the
        # relevant GFA camera's airmass

        zp_clear_adu_per_s = np.full(len(self.pmgstars), np.nan)

        for i, row in enumerate(self.pmgstars):
            zp_clear_adu_per_s[i] = util.zp_photometric_at_airmass(row['GFA_LOC'], row['airmass_per_gfa'])

        self.pmgstars['zp_clear_adu_per_s'] = zp_clear_adu_per_s

    def empty_pmgstars_ccds_columns(self):
        # for the case when user specifes --pmgstars option
        # but no PMGSTARS table exists, put these placeholder values
        # into the _ccds table

        self.ccds['n_pmgstars_all'] = 0
        self.ccds['n_pmgstars_retained'] = 0
        self.ccds['fiberfac'] = np.nan
        self.ccds['fiberfac_elg'] = np.nan
        self.ccds['fiberfac_bgs'] = np.nan

    def pmgstars_forcedphot(self):

        # think this may only happen for junk test data taken during daytime?
        # example is 76356 from observing night 20210217
        if self.pmgstars is None:
            self.empty_pmgstars_ccds_columns()
            return

        # driver for PMGSTARS forced photometry

        x, y = util.row_col_to_xy(self.pmgstars)

        self.pmgstars['expid'] = self.exp_header['EXPID']
        self.pmgstars['cube_index'] = self.ccds['cube_index'][0]
        self.pmgstars['night'] = self.ccds['night'][0]
        self.pmgstars['xcentroid'] = x
        self.pmgstars['ycentroid'] = y
        self.pmgstars['min_edge_dist_pix'] = [util.min_edge_dist_pix(c[0], c[1]) for c in zip(x, y)]

        # decide which stars to retain for the forced photometry analysis

        self._trim_pmgstars()
        self.pmgstars_dq_flags()

        good = (self.pmgstars['median_1_'] > 0) & (self.pmgstars['ang_sep_deg'] < 2.0/3600.0) & (self.pmgstars['min_edge_dist_pix'] >= 10) & (self.pmgstars['dq_flags'] == 0)

        self.pmgstars['good'] = good.astype(int)

        self.pmgstars_airmass_exptime()

        self.pmgstars_zp_clear()

        r_ps1 = -2.5*np.log10(self.pmgstars['median_1_'])

        self.pmgstars['r_mag_ps1_median'] = r_ps1

        clear_total_flux_adu_pred = np.power(10.0, (self.pmgstars['zp_clear_adu_per_s'] - r_ps1)/2.5)*self.pmgstars['time_s_for_dark']

        clear_total_flux_adu_pred[np.logical_not(good)] = np.nan

        self.pmgstars['clear_total_flux_adu_pred'] = clear_total_flux_adu_pred

        par = common.gfa_misc_params()

        self.pmgstars['fiber_flux_nominal_adu_pointsource'] = clear_total_flux_adu_pred*par['fracflux_nominal_pointsource']
        self.pmgstars['fiber_flux_nominal_adu_elg'] = clear_total_flux_adu_pred*par['fracflux_nominal_elg']
        self.pmgstars['fiber_flux_nominal_adu_bgs'] = clear_total_flux_adu_pred*par['fracflux_nominal_bgs']

        # these will be the forced photometry aperture fluxes in ADU
        aper_fluxes = np.full(len(self.pmgstars), np.nan)
        aper_fluxes_elg = np.full(len(self.pmgstars), np.nan)
        aper_fluxes_bgs = np.full(len(self.pmgstars), np.nan)

        # get list of extnames that will actually require
        # forced aperture photometry (i.e. those with at least
        # one PMGSTARS row that passes all quality cuts)
        extnames = np.unique(self.pmgstars[good]['GFA_LOC'])

        for extname in extnames:
            mask = ((self.pmgstars['GFA_LOC'] == extname) & good)

            assert(np.sum(mask) > 0)

            # do i need to send in the dq_flags image here? don't think so...
            fluxes = phot.pmgstars_forced_phot(self.pmgstars[mask]['xcentroid'],
                                               self.pmgstars[mask]['ycentroid'],
                                               self.images[extname].image)

            aper_fluxes[mask] = fluxes

            fluxes_elg = phot.pmgstars_forced_phot(self.pmgstars[mask]['xcentroid'],
                                                   self.pmgstars[mask]['ycentroid'],
                                                   self.images[extname].image,
                                                   elg=True)

            fluxes_bgs = phot.pmgstars_forced_phot(self.pmgstars[mask]['xcentroid'],
                                                   self.pmgstars[mask]['ycentroid'],
                                                   self.images[extname].image,
                                                   bgs=True)

            aper_fluxes_elg[mask] = fluxes_elg
            aper_fluxes_bgs[mask] = fluxes_bgs

            # call phot.pmgstars_forced_phot, filling in aper_fluxes, aper_fluxes_elg
            # will there be a second call to phot.pmgstars_forced_phot
            # or will one call include both ELG-smoothed version and
            # non-smoothed version?
            #    seems like it may be more efficient for one call to do
            #    both
            # what will be the inputs/outputs to phot.pmgstars_forced_phot ?
            #    definitely need xcentroid, ycentroid - what about extname?

        self.pmgstars['fiber_flux_adu_forced'] = aper_fluxes
        self.pmgstars['fiber_flux_adu_forced_elg'] = aper_fluxes_elg
        self.pmgstars['fiber_flux_adu_forced_bgs'] = aper_fluxes_bgs

        # then calculate the fiber-sized aperture throughput factors relative to nominal
        fiberfac = np.array(self.pmgstars['fiber_flux_adu_forced']/self.pmgstars['fiber_flux_nominal_adu_pointsource'])
        self.pmgstars['fiberfac'] = fiberfac

        fiberfac_elg = np.array(self.pmgstars['fiber_flux_adu_forced_elg']/self.pmgstars['fiber_flux_nominal_adu_elg'])
        fiberfac_bgs = np.array(self.pmgstars['fiber_flux_adu_forced_bgs']/self.pmgstars['fiber_flux_nominal_adu_bgs'])

        self.pmgstars['fiberfac_elg'] = fiberfac_elg
        self.pmgstars['fiberfac_bgs'] = fiberfac_bgs

        # now augment the CCDs table with per-camera
        # fiberfac, fiberfac_elg information

        self.ccds['n_pmgstars_all'] = [int(np.sum(self.pmgstars['GFA_LOC'] == extname)) for extname in self.ccds['camera']]
        self.ccds['n_pmgstars_retained'] = [int(np.sum((self.pmgstars['GFA_LOC'] == extname) & good)) for extname in self.ccds['camera']]

        self.ccds['fiberfac'] = [np.nanmedian(fiberfac[(self.pmgstars['GFA_LOC'] == extname) & good]) for extname in self.ccds['camera']]
        self.ccds['fiberfac_elg'] = [np.nanmedian(fiberfac_elg[(self.pmgstars['GFA_LOC'] == extname) & good]) for extname in self.ccds['camera']]
        self.ccds['fiberfac_bgs'] = [np.nanmedian(fiberfac_bgs[(self.pmgstars['GFA_LOC'] == extname) & good]) for extname in self.ccds['camera']]
