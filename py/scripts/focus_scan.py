#!/usr/bin/env python

import astropy.io.fits as fits
import numpy as np
import os
import matplotlib.pyplot as plt
import argparse

def focus_plots(night, expids,
               basedir='/n/home/datasystems/users/ameisner/reduced/realtime'):

    plt.figure(figsize=(12, 9))
    extnames = ['GUIDE0', 'GUIDE2', 'GUIDE3', 'GUIDE5', 'GUIDE7', 'GUIDE8']

    focus_z = []
    fwhm_pix = []

    # PSF stamps plot
    plt.subplots_adjust(hspace=0.01, wspace=0.01)
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
            
        hdul = fits.open(fname)
        extnames_present = [hdu.header['EXTNAME'] for hdu in hdul]
        for j, extname in enumerate(extnames):
            if extname not in extnames_present:
                continue
            print(i, j)
            plt.subplot(6, len(expids), len(expids)*j + i +  1)
            plt.xticks([])
            plt.yticks([])
            im = fits.getdata(fname, extname=extname)
            plt.imshow(im, interpolation='nearest', origin='lower', cmap='gray_r', vmin=0.01)
            plt.text(5, 44, str(expid) + '; ' + extname, color='r', fontsize=9)

            if np.isfinite(ccds[j]['XCENTROID_PSF']) and np.isfinite(ccds[j]['YCENTROID_PSF']):
                plt.scatter([ccds[j]['XCENTROID_PSF']], [ccds[j]['YCENTROID_PSF']], marker='.', c='r')

    
    expid_min = int(np.min(expids))

    print(focus_z)
    print(fwhm_pix)

    plt.savefig('stamps_focus_scan-' + str(expid_min).zfill(8)+'.png', bbox_inches='tight')
    plt.cla()

    plt.figure()
    
    asec_per_pix = 0.205

    focus_z = np.array(focus_z)
    fwhm_asec = np.array(fwhm_pix)*asec_per_pix
    plt.scatter(focus_z, fwhm_asec)
    plt.xlabel('focus z (micron)')
    plt.ylabel('FWHM (asec)')

    coeff = np.polyfit(focus_z, fwhm_asec, 2)

    xsamp = np.arange(np.min(focus_z), np.max(focus_z))
    ysamp = coeff[0]*(np.power(xsamp, 2)) + coeff[1]*xsamp + coeff[2]

    plt.title('focus scan starting with EXPID = ' + str(expid_min))

    
    plt.plot(xsamp, ysamp)

    zmin = -coeff[1]/(2*coeff[0])

    min_fwhm_fit_asec = coeff[0]*(zmin**2) + coeff[1]*zmin + coeff[2]
    
    yrange = [np.min(fwhm_asec), np.max(fwhm_asec)]
    plt.text(focus_z[2], yrange[0] + 0.8*(yrange[1]-yrange[0]), 'best FWHM (meas) : ' + '{:.2f}'.format(np.min(fwhm_asec)))
    plt.text(focus_z[2], yrange[0] + 0.7*(yrange[1]-yrange[0]), 'best FWHM (fit) : ' + '{:.2f}'.format(min_fwhm_fit_asec))
    plt.text(focus_z[2], yrange[0] + 0.9*(yrange[1]-yrange[0]), 'best focus : ' + str(int(np.round(zmin))))
    
    plt.savefig('fit_focus_scan-' + str(expid_min) + '.png', bbox_inches='tight')
    
def _test():
    night = '20200131'
    expids = 45446 + np.arange(7)

    focus_plots(night, expids, basedir='/project/projectdirs/desi/users/ameisner/GFA/run/psf_flux_weighted_centroid')

def _test_missing_cam():
    night = '20200131'
    expids = 45485 + np.arange(7)

    focus_plots(night, expids, basedir='/project/projectdirs/desi/users/ameisner/GFA/run/psf_flux_weighted_centroid')

if __name__ == "__main__":
    descr = 'GFA focus sequence plots/analysis'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('first_expid', type=int, nargs=1)

    parser.add_argument('night', type=str, nargs=1)
    
    parser.add_argument('--basedir', default='/n/home/datasystems/users/ameisner/reduced/realtime',
                        type=str, help='base directory for GFA reductions')

    args = parser.parse_args()

    expids = args.first_expid + np.arange(7, dtype=int)

    print(expids)
    print(args.night[0])
    print(args.basedir)
    
    focus_plots(args.night[0], expids, basedir=args.basedir)