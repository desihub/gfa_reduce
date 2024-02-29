import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os

fname = os.path.join(os.environ['CI_REDUCE_ETC'], 'desi-tiles.fits')

tab = fits.getdata(fname)

tab = tab[(tab['IN_DESI'] == 1) & (tab['PASS'] == 0)]

plt.figure(figsize=(5.5, 2.25))

plt.scatter(tab['RA'], tab['DEC'], s=2, edgecolor='none')

plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')

plt.xlim((365, -5))

plt.title('CI data challenge pointings')

plt.savefig('ci_data_challenge_pointings.png', bbox_inches='tight')
