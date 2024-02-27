from focus_scan import focus_plots
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt

fname_inv = '/global/homes/a/ameisner/gfa/pro/focus_scan_inventory-thru_20200315.fits'

def plot_one_row(ind, inv=None, outdir=None, basedir=None,
                 dont_plot_centroid=False, n_stars_min=3,
                 skip_low_n_stamps=True, extnames_exclude=[]):

    if inv is None:
        inv = fits.getdata(fname_inv)

    row = inv[ind]

    night = row['NIGHT']
    expids = row['EXPIDS'][row['EXPIDS'] != -1]

    if basedir is None:
        basedir = '/global/cfs/cdirs/desi/users/ameisner/GFA/reduced/v0017'

    if outdir is None:
        outdir = '/global/cscratch1/sd/ameisner/focus_plots'
        if dont_plot_centroid:
            outdir = outdir + '_no_centroid'

    print('Working on index : ', ind)
    focus_plots(night, expids, basedir=basedir, outdir=outdir,
                no_popups=True, dont_plot_centroid=dont_plot_centroid,
                n_stars_min=n_stars_min,
                skip_low_n_stamps=skip_low_n_stamps,
                extnames_exclude=extnames_exclude, write_coeff=True)

def _loop(indstart=0, nproc=None, skip_low_n_stamps=True):

    inv = fits.getdata(fname_inv)

    if nproc is None:
        nproc = len(inv)

    indend = min(indstart + nproc, len(inv))
    
    for ind in np.arange(indstart, indend):
        plot_one_row(ind, inv=inv, skip_low_n_stamps=skip_low_n_stamps)
        plt.cla()
