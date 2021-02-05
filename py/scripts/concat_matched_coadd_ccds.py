import astropy.io.fits as fits
from astropy.table import Table, vstack, hstack
import glob
import numpy as np
import os
import argparse

basedir = '/global/cfs/cdirs/desi/users/ameisner/GFA/reduced/v0022_matched'

def cube_index_median(tab, extra_cuts=False):

    # make quality cuts

    bad = (tab['MAX_Q1'] == 0) | (tab['MAX_Q2'] == 0) | (tab['MAX_Q3'] == 0) | (tab['MAX_Q4'] == 0) | (tab['NPIX_BAD_TOTAL'] >= 10)

    is_coadd = np.max(tab['CUBE_INDEX']) == -1

    contrast_thresh = 1.85 if is_coadd else 2.0

    if extra_cuts:
        bad = np.logical_or(bad, tab['CONTRAST'] < contrast_thresh)
        bad = np.logical_or(bad, tab['N_SOURCES_FOR_PSF'] < 3)


    print(len(tab))
    tab = tab[np.logical_not(bad)]
    print(len(tab))
    
    _id = np.array([str(tab[i]['EXPID']).zfill(8) + str(tab[i]['CUBE_INDEX']).zfill(5) for i in range(len(tab))])

    ids_u = np.unique(_id)

    print(len(ids_u), ' !!!!!!')

    rows = []
    for id_u in ids_u:
        _tab = tab[_id == id_u]

        row = Table()

        for colname in ['EXPID', 'CUBE_INDEX', 'NIGHT', 'EXPTIME', 'FNAME_RAW', 'SKYRA', 'SKYDEC', 'PROGRAM', 'MOON_ILLUMINATION', 'MOON_ZD_DEG', 'MOON_SEP_DEG', 'KTERM']:
            row[colname] = [_tab[0][colname]]

        for colname in ['MJD', 'FWHM_ASEC', 'TRANSPARENCY', 'SKY_MAG_AB', 'FIBER_FRACFLUX', 'FIBER_FRACFLUX_ELG', 'AIRMASS', 'RADPROF_FWHM_ASEC']:
            row[colname] = [np.nanmedian(_tab[colname])]

        rows.append(row)

    print(len(rows))
    result = vstack(rows)

    return result

def _nights_list(night_min, night_max, basedir=basedir):

    dirs = glob.glob(basedir + '/????????')

    _dirs = np.array([os.path.split(d)[-1] for d in dirs])

    dirs = np.array(dirs)
    dirs = dirs[(_dirs >= night_min) & (_dirs <= night_max)]

    nights = [os.path.split(d)[-1] for d in dirs]
    return nights

def _concat(night='20201214', basedir=basedir):

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

def _concat_many_nights(night_min='20201214', night_max='99999999',
                        basedir=basedir):

    nights = _nights_list(night_min, night_max, basedir=basedir)

    tables = []
    for night in nights:
        table = _concat(night=night, basedir=basedir)
        tables.append(table)

    result = vstack(tables)

    sind = np.argsort(result['EXPID'] + 0.1*result['PETAL_LOC'])

    result = result[sind]

    med = cube_index_median(result)
    _med = cube_index_median(result, extra_cuts=True)
    
    return result, med, _med

def _write_many_nights(night_min='20201214', night_max='99999999',
                       basedir=basedir):

    result, med, _med = _concat_many_nights(night_min=night_min,
                                            night_max=night_max,
                                            basedir=basedir)
    
    hdul = fits.HDUList(hdus=[fits.PrimaryHDU(),
                              fits.BinTableHDU(data=result),
                              fits.BinTableHDU(data=med),
                              fits.BinTableHDU(data=_med)])

    night = str(np.max(result['NIGHT']))
    outname = 'offline_matched_coadd_ccds_SV1-thru_' + night + '.fits'
    hdul.writeto(outname)

if __name__=="__main__":
    descr = 'gather gfa_reduce _ccds table outputs'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('--basedir', default=basedir, type=str,
                        help='input directory')

    parser.add_argument('--night_min', default='20201214', type=str,
                        help='first observing night')

    parser.add_argument('--night_max', default='99999999', type=str,
                        help='last observing night')

    args = parser.parse_args()

    _write_many_nights(night_min=args.night_min, night_max=args.night_max,
                       basedir=args.basedir)

