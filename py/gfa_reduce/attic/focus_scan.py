#!/usr/bin/env python

import astropy.io.fits as fits
import numpy as np
import os
import matplotlib.pyplot as plt
import argparse

def _write_coeff(outdir, first_expid, focus_z, coeff, fwhm_asec, expid_used):
    from astropy.table import Table

    tab = Table()

    assert(len(coeff) == 3)
    assert(len(focus_z) >= 3)
    assert(len(fwhm_asec) == len(focus_z))
    assert(len(expid_used) == len(focus_z))

    tab['first_expid'] = [first_expid]
    tab['poly_coeff'] = [coeff]
    tab['n_fwhm_meas'] = [len(focus_z)]
    tab['z_min_fit'] = [np.min(focus_z)]
    tab['z_max_fit'] = [np.max(focus_z)]
    tab['meas_min_fwhm_asec'] = [np.min(fwhm_asec)]

    zmin = -coeff[1]/(2*coeff[0]) if (coeff[0] > 0) else np.nan

    min_fwhm_fit_asec = (coeff[0]*(zmin**2) + coeff[1]*zmin + coeff[2]) if np.isfinite(zmin) else np.nan

    tab['model_best_z'] = [zmin]
    tab['model_min_fwhm_asec'] = [min_fwhm_fit_asec]

    tab['minimum_scanned'] = [int((zmin >= np.min(focus_z)) and (zmin <= np.max(focus_z)))]
    
    # maybe add something about the residual RMS here ?

    fwhm_model_pred = coeff[0]*(np.power(focus_z, 2)) + coeff[1]*focus_z + coeff[2]

    diff = fwhm_asec - fwhm_model_pred

    rms = np.sqrt(np.mean(np.power(diff, 2)))

    tab['resid_rms_asec'] = [rms]

    # try to package in the actual list of (FWHM, expid, focus_z)  values
    n_exp_max = 21 # 3x longer focus sequence than usual
    n_pad = n_exp_max - len(fwhm_asec)

    tab['fwhm_asec_values'] = [np.concatenate([fwhm_asec, np.full(n_pad, np.nan)])]
    tab['expid_values'] = [np.concatenate([expid_used, np.full(n_pad, -1)])]
    tab['focus_z_values'] = [np.concatenate([focus_z, np.full(n_pad, np.nan)])]

    outname = 'poly_coeff-' + str(first_expid).zfill(8) + '.fits'

    outname = os.path.join(outdir, outname)

    tab.write(outname, format='fits')

def _trim_fit_inputs(expid_used, fwhm_pix, focus_z):
    if len(expid_used) <= 7:
        return expid_used, fwhm_pix, focus_z

    indmin = np.argmin(fwhm_pix)

    diff = np.abs(np.arange(len(expid_used)) - indmin)

    keep = (diff <= 4)

    expid_used = expid_used[keep]
    fwhm_pix = fwhm_pix[keep]
    focus_z = focus_z[keep]

    return expid_used, fwhm_pix, focus_z
    
def color_frame(sidelen, delta_pix=0, color='orange'):
    x = [0 + delta_pix, 0 + delta_pix, sidelen -1 - delta_pix, sidelen - 1 - delta_pix, 0 + delta_pix]
    y = [0 + delta_pix, sidelen - 1 - delta_pix, sidelen - 1 - delta_pix, 0 + delta_pix, 0 + delta_pix]

    plt.plot(x, y, c=color, linewidth=2)

def focus_plots(night, expids,
                basedir='/n/home/datasystems/users/ameisner/reduced/focus',
                outdir='/n/home/desiobserver/focus_scan', no_popups=False,
                dont_plot_centroid=False, n_stars_min=-1, skip_low_n_stamps=False,
                flag_bad_denoising=True, extnames_exclude=[], write_coeff=False):

    # do i also want a separate boolean keyword arg to leave out low N cases
    # from the parabola fits as well??

    plt.figure(1, figsize=(12.0*(len(expids)/7.0), 9))
    extnames = ['GUIDE0', 'GUIDE2', 'GUIDE3', 'GUIDE5', 'GUIDE7', 'GUIDE8']

    focus_z = []
    fwhm_pix = []
    expid_used = []

    # PSF stamps plot
    plt.subplots_adjust(hspace=0.01, wspace=0.01)
    n_stamps_plotted = 0
    for i, expid in enumerate(expids):
        fname = basedir + '/' + night + '/' + str(expid).zfill(8) + '/gfa-' + str(expid).zfill(8) + '_psfs.fits'
        print(fname)
        fname_ccds = fname.replace('_psfs.fits', '_ccds.fits')
        if not os.path.exists(fname):
            continue
        ccds = fits.getdata(fname_ccds)

        if np.sum(np.isfinite(ccds['PSF_FWHM_PIX'])) != 0:
            fwhm_pix.append(np.median(ccds['PSF_FWHM_PIX'][np.isfinite(ccds['PSF_FWHM_PIX'])]))
            focus_z.append(float(ccds[0]['FOCUS'].split(',')[2]))
            expid_used.append(expid)
            
        hdul = fits.open(fname)
        extnames_present = [hdu.header['EXTNAME'] for hdu in hdul]
        for j, extname in enumerate(extnames):
            if (extname not in extnames_present) or (extname in extnames_exclude):
                continue
            print(i, j)

            im, h = fits.getdata(fname, extname=extname, header=True)

            n_stars = h['NSTARS']

            if n_stars < n_stars_min:
                print('expid = ', h['EXPID'], ' ; extname = ', h['EXTNAME'], ' has too few contributing sources')
                if skip_low_n_stamps:
                    continue

            plt.subplot(6, len(expids), len(expids)*j + i +  1)
            plt.xticks([])
            plt.yticks([])

            plt.imshow(im, interpolation='nearest', origin='lower', cmap='gray_r', vmin=0.01)
            n_stamps_plotted += 1

            if n_stars < n_stars_min:
                color_frame(im.shape[0])

            ccd = ccds[ccds['EXTNAME'].replace(' ', '') == extname]

            assert(len(ccd) == 1)

            if flag_bad_denoising and (ccd[0]['NPIX_BAD_TOTAL'] >= 10):
                delta_pix = 1 if (n_stars < n_stars_min) else 0
                color_frame(im.shape[0], delta_pix=delta_pix, color='#4ff222')

            plt.text(5, 44, str(expid) + '; ' + extname, color='r', fontsize=9)
            plt.text(10, 3.5, 'z = ' + str(int(float(ccds[0]['FOCUS'].split(',')[2]))), color='r')
            
            if np.isfinite(ccds[j]['XCENTROID_PSF']) and np.isfinite(ccds[j]['YCENTROID_PSF']) and (not dont_plot_centroid):
                plt.scatter([ccds[j]['XCENTROID_PSF']], [ccds[j]['YCENTROID_PSF']], marker='.', c='r')

    
    expid_min = int(np.min(expids))

    print(focus_z)
    print(fwhm_pix)

    fwhm_pix = np.array(fwhm_pix)
    focus_z = np.array(focus_z)
    expid_used = np.array(expid_used)

    expid_used, fwhm_pix, focus_z = _trim_fit_inputs(expid_used, fwhm_pix, focus_z)
    
    if n_stamps_plotted > 0:
        plt.savefig(os.path.join(outdir, 'stamps_focus_scan-' + str(expid_min).zfill(8)+'.png'), bbox_inches='tight')
    else:
        print('WARNING : NO GOOD PSF MODELS TO PLOT FOR ANY IMAGES !!!')

    # doesn't make sense to fit a parabola to < 3 data points...
    if len(focus_z) < 3:
        print('NOT ENOUGH EXPOSURES AVAILABLE TO FIT A PARABOLA')
        return

    asec_per_pix = 0.205
    
    if (n_stamps_plotted == 0) and (np.min(fwhm_pix*asec_per_pix) < 0.4 ):
        return
    
    plt.figure(200)


    fwhm_asec = fwhm_pix*asec_per_pix
    plt.scatter(focus_z, fwhm_asec)
    plt.xlabel('focus z (micron)')
    plt.ylabel('FWHM (asec)')

    coeff = np.polyfit(focus_z, fwhm_asec, 2)

    xsamp = np.arange(np.min(focus_z), np.max(focus_z))
    ysamp = coeff[0]*(np.power(xsamp, 2)) + coeff[1]*xsamp + coeff[2]

    plt.title('focus scan starting with EXPID = ' + str(expid_min))

    
    plt.plot(xsamp, ysamp)

    
    yrange = [np.min(fwhm_asec), np.max(fwhm_asec)]
    x_text = (np.max(focus_z) - np.min(focus_z))*0.4 + np.min(focus_z)
    plt.text(x_text, yrange[0] + 0.8*(yrange[1]-yrange[0]), 'best FWHM (meas) : ' + '{:.2f}'.format(np.min(fwhm_asec)))

    # only calculate the model-based minimum FWHM value if parabola opens upward...
    if coeff[0] > 0:
        zmin = -coeff[1]/(2*coeff[0])
        min_fwhm_fit_asec = coeff[0]*(zmin**2) + coeff[1]*zmin + coeff[2]
        plt.text(x_text, yrange[0] + 0.7*(yrange[1]-yrange[0]), 'best FWHM (fit) : ' + '{:.2f}'.format(min_fwhm_fit_asec))
        plt.text(x_text, yrange[0] + 0.9*(yrange[1]-yrange[0]), 'best focus : ' + str(int(np.round(zmin))))
    else:
        print('WARNING: BEST-FIT PARABOLA DOES NOT OPEN UPWARD')
        plt.text(x_text, yrange[0] + 0.9*(yrange[1]-yrange[0]), 'best focus : N/A')
    
    plt.savefig(os.path.join(outdir, 'fit_focus_scan-' + str(expid_min).zfill(8) + '.png'), bbox_inches='tight')
    if not no_popups:
        plt.show()

    if write_coeff:
        _write_coeff(outdir, np.min(expids), focus_z, coeff, fwhm_asec, expid_used)
    
def _test():
    night = '20200131'
    expids = 45446 + np.arange(7)

    focus_plots(night, expids, basedir='/project/projectdirs/desi/users/ameisner/GFA/run/psf_flux_weighted_centroid', outdir='.')

def _test_missing_cam():
    night = '20200131'
    expids = 45485 + np.arange(7)

    focus_plots(night, expids, basedir='/project/projectdirs/desi/users/ameisner/GFA/run/psf_flux_weighted_centroid')

if __name__ == "__main__":
    descr = 'GFA focus sequence plots/analysis'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('first_expid', type=int, nargs=1)

    parser.add_argument('night', type=str, nargs=1)
    
    parser.add_argument('--basedir', default='/n/home/datasystems/users/ameisner/reduced/focus',
                        type=str, help='base directory for GFA reductions')

    parser.add_argument('--outdir', default='/n/home/desiobserver/focus_scan', 
                        type=str, help='output directory for plot PNGs')

    parser.add_argument('--no_popups', default=False, action='store_true',
                        help='write PNGs without popping up plot windows')

    parser.add_argument('--dont_plot_centroid', default=False, action='store_true',
                        help='do not overplot centroid location on each postage stamp')

    args = parser.parse_args()

    expids = args.first_expid + np.arange(7, dtype=int)

    print(expids)
    print(args.night[0])
    print(args.basedir)
    
    outdir = args.outdir if os.path.exists(args.outdir) else '.'
    focus_plots(args.night[0], expids, basedir=args.basedir, outdir=outdir, no_popups=args.no_popups,
                dont_plot_centroid=args.dont_plot_centroid)
