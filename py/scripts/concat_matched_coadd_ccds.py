import astropy.io.fits as fits
from astropy.table import Table, vstack, hstack
import glob
import numpy as np
import os

basedir = '/global/cfs/cdirs/desi/users/ameisner/GFA/reduced/v0022_matched'

def _nights_list(night_min, night_max, basedir=basedir):

    dirs = glob.glob(basedir + '/????????')

    _dirs = np.array([os.path.split(d)[-1] for d in dirs])

    dirs = np.array(dirs)
    dirs = dirs[(_dirs >= night_min) & (_dirs <= night_max)]

    nights = [os.path.split(d)[-1] for d in dirs]
    return nights

def _concat(night='20201214'):

    pattern = basedir + '/' + night + '/????????/*ccds*.fits'

    flist = glob.glob(pattern)

    tables = []
    for i, f in enumerate(flist):
        print(i, f)
        tab = fits.getdata(f)
        tab = Table(tab)
        tab['fwhm_asec'] = tab['psf_fwhm_asec']
    
        tables.append(tab)

    result = vstack(tables)

    # for matched coadds during SV1 this is always true...
    result['spectro_expid'] = result['expid']

    for name in result.colnames:
        result.rename_column(name, name.upper())

    return result

def _concat_many_nights(night_min='20201214', night_max='20210203'):

    nights = _nights_list(night_min, night_max)

    tables = []
    for night in nights:
        table = _concat(night=night)
        tables.append(table)

    result = vstack(tables)

    sind = np.argsort(result['EXPID'] + 0.1*result['PETAL_LOC'])

    result = result[sind]
    
    return result
