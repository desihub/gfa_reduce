import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

fname = '/global/cscratch1/sd/ameisner/M57_clean.fits'

im = fits.getdata(fname)

# linear stretch, for whole image
#vmin = -38.6
#vmax = 811.35

im[im < 25.0] = 25.0 #11.5
im[im > 3100.0] = 3100.0 # 3082.0

_im = np.log10(im)

# y: 325 -> 1700
# x: 800 -> 2300

_im = _im[325:1700, 800:2249]

plt.imshow(_im, cmap='GnBu_r', interpolation='nearest')

plt.gca().get_xaxis().set_visible(False)
plt.gca().get_yaxis().set_visible(False)

plt.gca().set_axis_off()

plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
plt.margins(0,0)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())

dpi = 400

plt.savefig('/project/projectdirs/cosmo/www/temp/ameisner/DESI_CI_ring_nebula-GnBu_r-'+ str(dpi) + 'dpi.png', dpi=dpi, bbox_inches='tight', pad_inches=-0.02)
