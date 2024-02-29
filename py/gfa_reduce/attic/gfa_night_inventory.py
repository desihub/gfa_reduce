#!/usr/bin/env python3

import argparse
import astropy.io.fits as fits
import glob

def _printouts(night, full_filename=False, imagecam=False, guider=False):
    basedir = '/exposures/desi'
    nightdir = basedir + '/' + night
    pattern = 'gfa*.fits.fz' if not guider else 'guide-????????.fits.fz'
    flist = glob.glob(nightdir + '/*/' + pattern)

    flist.sort()

    extname = 'GFA' if not guider else 'GUIDER'

    n_frames_tot = 0
    for f in flist:
        h = fits.getheader(f, extname=extname)
        fp = (f.replace(basedir + '/', '') + '   ') if not full_filename else f + '   '
        if not guider:
            print(fp, h['FLAVOR'].ljust(10, ' '), '{:.1f}'.format(h['EXPTIME']).ljust(10, ' '), h['PROGRAM'].ljust(40), '' if not imagecam else h['IMAGECAM'])
        else:
            n_frames_tot += h['FRAMES']
            print(fp, str(h['FRAMES']).ljust(10, ' '), h['PROGRAM'].ljust(40))

    if guider:
        print('total # of guider frames = ' +str(n_frames_tot))

if __name__ == "__main__":
    descr = 'print inventory of GFA data for a night'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('night', type=str, nargs=1,
                        help='NIGHT string')

    parser.add_argument('--full_filename', default=False,
                        action='store_true', 
                        help='print full gfa*.fz filename path')

    parser.add_argument('--imagecam', default=False,
                        action='store_true',
                        help='also print IMAGECAM')

    parser.add_argument('--guider', default=False, 
                        action='store_true', 
                        help='guide-????????.fits.fz rather than gfa*.fits.fz')

    args = parser.parse_args()

    night = args.night[0]

    _printouts(night, full_filename=args.full_filename, imagecam=args.imagecam, guider=args.guider)
