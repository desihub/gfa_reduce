import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

fname = '/global/cscratch1/sd/ameisner/firstlight_clean.long.fits'


im = fits.getdata(fname)

# linear stretch, for whole image
#vmin = -38.6
#vmax = 811.35

im[im < 109] = 109
im[im > 60000] = 60000

_im = np.log10(im)

plt.imshow(_im, cmap='viridis')

plt.gca().get_xaxis().set_visible(False)
plt.gca().get_yaxis().set_visible(False)

plt.gca().set_axis_off()

plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
plt.margins(0,0)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())

dpi = 400

plt.savefig('/project/projectdirs/cosmo/www/temp/ameisner/DESI_CI_firstlight-viridis-'+ str(dpi) + 'dpi.png', dpi=dpi, bbox_inches='tight', pad_inches=-0.02)
