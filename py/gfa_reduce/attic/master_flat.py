import ci_reduce.common as common
import ci_reduce.imred.load_calibs as load_calibs
import astropy.io.fits as fits
import glob
import numpy as np
import os

# only CIN had its lens cap off in the available engineering data
# so use the CIN flat for all cameras

# ignore dark current here, given the low dark current rates quoted in
# DESI-3358 and the very short integration times (1-5 seconds)

engdir = '/project/projectdirs/desi/users/ameisner/CI/FORDK'
ci_extname = 'CIN'

def get_flat_frames_names(setnum):
    assert((setnum == 1) or (setnum == 2))
    set_string = 'set' + str(setnum)

    flist = glob.glob(engdir + '/*/*.fits.fz')
    flist = np.array(flist)
    nobin_pos = np.array([f.find('NOBIN') for f in flist])
    flist = flist[nobin_pos != -1]

    acttime = np.zeros(len(flist))
    for i, f in enumerate(flist):
        h = fits.getheader(f)
        acttime[i] = h['ACTTIME']

    flist = flist[acttime != 0]
    flist = flist[[(f.find(set_string) != -1) for f in flist]]
    return flist

def read_1flat(fname):

    assert(os.path.exists(fname))
    # again, only working with CIN for now, since others have no relevant data
    flat, header = fits.getdata(fname, extname=ci_extname, header=True)
    flat = flat.astype('float')
    flat -= load_calibs.read_bias_image(ci_extname)

    return flat, header

def countrate_1flat(fname):
    # compute the count rate for one flat field image and also
    # the corresponding inverse variance

    # flat is already bias-subtracted
    flat, header = read_1flat(fname)

    t = header['ACTTIME']

    # ADU/s
    countrate = flat/t

    gain = common.ci_camera_gain(ci_extname)

    # should probably also include a term for the random noise in the
    # master bias, which currently isn't all that small because only
    # 4 bias frames are available

    flat_var_adu_sq = (common.ci_camera_readnoise(ci_extname)**2 + \
                       flat*gain)/(gain**2)

    # make countrate_var a scalar to avoid per-pixel uncertainties that
    # would be biased high (low) for pixels that randomly scattered high (low)
    countrate_var = np.median(flat_var_adu_sq)/(t**2)
    countrate_ivar = 1.0/countrate_var
    
    return countrate, countrate_ivar

def estimate_flatfield(setnum):
    assert((setnum == 1) or (setnum == 2))

    flist = get_flat_frames_names(setnum)

    par = common.ci_misc_params()
    countrate_tot = np.zeros((par['height_pix_native'],
                              par['width_pix_native']))
    countrate_wt = 0.0
    for f in flist:
        countrate, countrate_ivar = countrate_1flat(f)
        countrate_tot += countrate*countrate_ivar
        countrate_wt += countrate_ivar

    countrate_est = countrate_tot/countrate_wt

    # normalize flat field to a median value of 1
    countrate_med = np.median(countrate_est)

    countrate_norm = countrate_est/countrate_med
    countrate_norm_ivar = countrate_wt*(countrate_med**2)

    # note that countrate_norm is an image whereas countrate_norm_ivar is
    # a scalar value
    return countrate_norm, countrate_norm_ivar
     
# could imagine adding a "second round" step to the flat field determination
# process where I use the initial "naive" flat field to mask outliers in
# the individual contributing frames

# might be able to use this to identify any cosmics in the CIN light
# exposures and see what they look like

def create_master_flat():
    # can't assume illumination level was the same for set1, set2
    # "light" exposures, so use each set separately to build a flat field
    # then combine the two by doing a weighted average, and sum the
    # inverse variances

    flat1, ivar1 = estimate_flatfield(1)
    flat2, ivar2 = estimate_flatfield(2)

    ivar = ivar1 + ivar2
    flat = (flat1*ivar1 + flat2*ivar2)/ivar

    # for good measure renormalize very slightly so that median is exactly 1

    med = np.median(flat)
    flat = flat/med
    ivar = ivar/med

    return flat, ivar

def master_flat_header_cards(hdu, ci_extname, ivar):
    h = hdu.header

    # not sure if this is otherwise considered an allowable flavor, but
    # I think it's important to emphasize that these are dimensionless images
    # as opposed to FLAVOR = "LIGHT" images which have units of ADU
    h['FLAVOR'] = 'FLAT'
    h['EXTNAME'] = ci_extname
    h['IVAR'] = ivar

    return h

def write_master_flat():
    
    par = common.ci_misc_params()
    outname = par['master_flat_filename']
    outname = os.path.join(os.environ[par['etc_env_var']], outname)

    assert(not os.path.exists(outname))

    flat, ivar = create_master_flat()

    flat = flat.astype('float32')
    ci_extnames = common.valid_image_extname_list()

    hdus = []
    for ci_extnum in range(len(ci_extnames)):
        ci_extname = common.ci_extnum_to_extname(ci_extnum, fz=False)
        print('Assembling HDU for: ' + ci_extname)
        if len(hdus) == 0:
            hdu = fits.PrimaryHDU(flat)
        else:
            hdu = fits.ImageHDU(flat)
        hdu.header = master_flat_header_cards(hdu, ci_extname, ivar)

        hdus.append(hdu)

    hdul = fits.HDUList(hdus)

    hdul.writeto(outname)
