import astropy.io.fits as fits
import glob
import numpy as np
import os

# directories that are after klaus' e-mail where he said that 
# CCDTEMP had been added, which was: Sun, Mar 31, 1:59 PM
dirnames = ['20190331', '20190401', '20190402', '20190403']
extnames = ['CIE', 'CIN', 'CIC', 'CIS', 'CIW']

def get_flist():
    flist = []
    for d in dirnames:
        pattern = '/project/projectdirs/desi/spectro/data/' + d + '/*/ci*.fits.fz'
        #print(pattern)
        l = glob.glob(pattern)
        flist = flist + l

    flist.sort()
    expid = np.array([int(os.path.basename(f)[3:11]) for f in flist])

    flist = np.array(flist)
    # 3093 is the first exposure for which Klaus said ccdtemp had been added
    flist = flist[expid >= 3093]
    return flist

# assume for now that all of CIE, CIS, CIW, CIN, CIC extensions are
# always present

def summarize_ccdtemps_1exp(fname):

    hdul = fits.open(fname)
    if len(hdul) != 7:
        return

    print('='*5 + ' ' +  fname + ' ' + '='*5)
    for extname in extnames:
        h = fits.getheader(fname, extname=extname)
        cards = [t[0] for t in h.cards]
        has_ccdtemp = ('CCDTEMP' in cards)
        has_string = ('has' if has_ccdtemp else 'no')
        if has_ccdtemp:
            str_extra = ', CCDTEMP = ' + str(h['CCDTEMP'])
        else:
            str_extra = ''
        print(extname + ': ' + has_string + ' CCDTEMP')

def _loop():
    flist = get_flist()
    for f in flist:
        summarize_ccdtemps_1exp(f)

if __name__=="__main__":
    _loop()
