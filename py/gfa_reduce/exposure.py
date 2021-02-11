import gfa_reduce.common as common
import gfa_reduce.imred.load_calibs as load_calibs
import gfa_reduce.dark_current as dark_current
import astropy.io.fits as fits
import numpy as np
import gfa_reduce.analysis.util as util
import gfa_reduce.analysis.phot as phot
from multiprocessing import Pool

class GFA_exposure:
    """Object encapsulating the contents of a single GFA exposure"""

    def __init__(self, image_list, exp_header=None, bintables=None,
                 max_cbox=31, pmgstars=None):
        # images is a dictionary of GFA_image objects

        par = common.gfa_misc_params()
        self.images = dict(zip(common.valid_image_extname_list(), 
                               par['n_cameras']*[None]))

        self.dark_current_objs = dict(zip(common.valid_image_extname_list(), 
                                          par['n_cameras']*[None]))

        self.assign_image_list(image_list)

        # exposure-level header
        self.exp_header = exp_header

        # hack for 20210106
        if self.exp_header is not None:
            if self.exp_header['SKYRA'] is None:
                print('REPLACING GUIDER SKYRA WITH REQRA')
                self.exp_header['SKYRA'] = self.exp_header['REQRA']
            if self.exp_header['SKYDEC'] is None:
                print('REPLACING GUIDER SKYDEC WITH REQDEC')
                self.exp_header['SKYDEC'] = self.exp_header['REQDEC']

        self.pixels_calibrated = None
        self.bintables = bintables
        self.max_cbox = max_cbox
        self.assign_max_cbox() # to the per-camera images ...
        self.pmgstars = pmgstars
        
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
        print('Attempting to subtract bias...')
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
        print('Attempting to subtract dark current...')

        if mp:
            args = []

        for extname in self.images.keys():
            if self.images[extname] is not None:
                
                acttime = self.images[extname].try_retrieve_meta_keyword('EXPTIME')
                if (acttime is None) or (acttime == 0):
                    print('trying REQTIME instead of EXPTIME')
                    acttime = self.images[extname].try_retrieve_meta_keyword('REQTIME')
                if (acttime is None) or (acttime == 0):
                    print('could not find an exposure time !!!!')
                    assert(False) # die

                self.images[extname].time_s_for_dark = acttime
                
                self.images[extname].t_c_for_dark_is_guess = False
                t_c = self.images[extname].try_retrieve_meta_keyword('GCCDTEMP')
                if t_c is None:
                    print('trying GCOLDTEC instead of GCCDTEMP')
                    t_c = self.images[extname].try_retrieve_meta_keyword('GCOLDTEC')
                    
                if t_c is None:
                    print('could not find a CCD temperature !!!!')
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
            print('Computing dark scalings for all guide cameras in parallel...')
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
        print('Attempting to apply flat field...')
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
            print('Skipping flatfielding')
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
                print('empirical ' + im.header['EXTNAME'] + 
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

            print('Running source cataloging for all guide cameras in parallel...')
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

    def try_retrieve_header_card(self, keyword, placeholder=None):
        if (keyword in self.exp_header.keys()) and (self.exp_header[keyword] is not None):
            return self.exp_header[keyword]
        else:
            print('could not find ' + keyword + ' in exposure-level header !!')
            return placeholder

    def compute_psfs(self, catalog):
        for image in self.images.values():
            if image is None:
                continue
            image.create_psf(catalog)

    def assign_input_filename(self, fname_in):
        self.fname_in = fname_in

    def pmgstars_dq_flags(self):
        # gather dq_flags for PMGSTARS table source positions

        extnames = np.unique(self.pmgstars['GFA_LOC'])

        dq_flags = np.zeros(len(self.pmgstars), dtype='uint8')

        for extname in extnames:
            mask = (self.pmgstars['GFA_LOC'] == extname)
            dq_flags[mask] = util.get_dq_flags(self.pmgstars[mask],
                                               self.images[extname].bitmask)

        self.pmgstars['dq_flags'] = dq_flags

    def pmgstars_forcedphot(self):
        # driver for PMGSTARS forced photometry

        x, y = util.row_col_to_xy(self.pmgstars)

        self.pmgstars['xcentroid'] = x
        self.pmgstars['ycentroid'] = y
        self.pmgstars['min_edge_dist_pix'] = [util.min_edge_dist_pix(c[0], c[1]) for c in zip(x, y)]

        # decide which stars to retain for the forced photometry analysis

        self.pmgstars_dq_flags()
        
        good = (self.pmgstars['median_1_'] > 0) & (self.pmgstars['ang_sep_deg'] < 2.0/3600.0) & (self.pmgstars['min_edge_dist_pix'] >= 10) & (self.pmgstars['dq_flags'] == 0)

        self.pmgstars['good'] = good.astype(int)
