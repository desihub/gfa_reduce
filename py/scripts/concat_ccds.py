#!/usr/bin/env python

import astropy.io.fits as fits
from astropy.table import Table, vstack, hstack
import glob
import numpy as np
import os
import argparse
from multiprocessing import Pool
import datetime

def _get_default_basedir(acq=False):
    basedir = os.environ['GFA_REDUX_BASE'] + '_' + \
              ('matched' if not acq else 'acq')

    return basedir

def cube_index_median(tab, extra_cuts=False, matched_coadd=True):

    # make quality cuts

    bad = (tab['MAX_Q1'] == 0) | (tab['MAX_Q2'] == 0) | (tab['MAX_Q3'] == 0) | (tab['MAX_Q4'] == 0) | (tab['NPIX_BAD_TOTAL'] >= 10)

    is_coadd = np.max(tab['CUBE_INDEX']) == -1

    contrast_thresh = 1.85 if is_coadd else 2.0

    if extra_cuts:
        bad = np.logical_or(bad, tab['CONTRAST'] < contrast_thresh)
        bad = np.logical_or(bad, tab['N_SOURCES_FOR_PSF'] < 3)

    tab = tab[np.logical_not(bad)]
    
    _id = np.array([str(tab[i]['EXPID']).zfill(8) + str(tab[i]['CUBE_INDEX']).zfill(5) + tab[i]['EXTNAME'] for i in range(len(tab))])

    sind = np.argsort(_id)
    _id = _id[sind]
    tab = tab[sind]
    del sind

    _id = np.array([str(tab[i]['EXPID']).zfill(8) + str(tab[i]['CUBE_INDEX']).zfill(5) for i in range(len(tab))])

    ids_u, ind_u = np.unique(_id, return_index=True)

    cols_to_median = ['MJD', 'FWHM_ASEC', 'TRANSPARENCY', 'SKY_MAG_AB',
                      'FIBER_FRACFLUX', 'FIBER_FRACFLUX_ELG',
                      'FIBER_FRACFLUX_BGS', 'AIRMASS', 'RADPROF_FWHM_ASEC']

    if matched_coadd:
        cols_to_median += ['FIBERFAC', 'FIBERFAC_ELG', 'FIBERFAC_BGS']
    
    rows = []
    for i, ind in enumerate(ind_u):
        ind_lower = ind
        ind_upper = ind_u[i+1] if (i != (len(ind_u)-1)) else len(tab)

        _tab = tab[ind_lower:ind_upper]

        row = Table()

        for colname in ['EXPID', 'CUBE_INDEX', 'NIGHT', 'EXPTIME',
                        'FNAME_RAW', 'SKYRA', 'SKYDEC', 'PROGRAM',
                        'MOON_ILLUMINATION', 'MOON_ZD_DEG', 'MOON_SEP_DEG',
                        'KTERM', 'FRACFLUX_NOMINAL_POINTSOURCE',
                        'FRACFLUX_NOMINAL_ELG', 'FRACFLUX_NOMINAL_BGS']:
            row[colname] = [_tab[0][colname]]

        for colname in cols_to_median:
            row[colname] = [np.nanmedian(_tab[colname])]

        row['MINCONTRAST'] = np.min(_tab['CONTRAST'])
        row['MAXCONTRAST'] = np.max(_tab['CONTRAST'])

        rows.append(row)

    result = vstack(rows)

    return result

def _nights_list(night_min, night_max, basedir=None, acq=False):

    if basedir is None:
        basedir = _get_default_basedir(acq=acq)

    dirs = glob.glob(basedir + '/????????')

    _dirs = np.array([os.path.split(d)[-1] for d in dirs])

    dirs = np.array(dirs)
    dirs = dirs[(_dirs >= night_min) & (_dirs <= night_max)]

    nights = [os.path.split(d)[-1] for d in dirs]

    nights.sort()
    return nights

def _read_one_ccds_table(fname):
    print('READING : ' + fname)
    tab = fits.getdata(fname)
    tab = Table(tab)
    tab['fwhm_asec'] = tab['psf_fwhm_asec']

    return tab

def _concat(night='20201214', basedir=None, acq=False, user_basedir=None,
            workers=1):

    if basedir is None:
        basedir = _get_default_basedir(acq=acq)

    # give precedence to user-generated version of this night's
    # gfa_reduce outputs, if those exist
    if user_basedir is not None:
        user_nightdir = os.path.join(user_basedir, night)
        if os.path.exists(user_nightdir):
            basedir = user_basedir

    suffix = '-0000_ccds.fits' if acq else '_ccds--0001.fits'
    pattern = basedir + '/' + night + '/????????/*' + suffix

    flist = glob.glob(pattern)

    if len(flist) == 0:
        return None

    print('initializing pool with ' + str(workers) + ' workers')
    p = Pool(workers)
    
    tables = p.map(_read_one_ccds_table, flist)

    p.close()
    p.join()

    result = vstack(tables)

    # for matched coadds during SV1 this is always true...
    if not acq:
        result['spectro_expid'] = result['expid']

    for name in result.colnames:
        result.rename_column(name, name.upper())

    return result

def _concat_many_nights(night_min='20201214', night_max='99999999',
                        basedir=None, acq=False, user_basedir=None, workers=1):

    if basedir is None:
        basedir = _get_default_basedir(acq=acq)

    print('reading in _ccds tables...')

    nights = _nights_list(night_min, night_max, basedir=basedir, acq=acq)

    if user_basedir is not None:
        nights_user = _nights_list(night_min, night_max, basedir=user_basedir,
                                   acq=acq)
        nights = list(set(nights_user + nights))
        nights.sort()

    print('first night is : ' + str(np.min(np.array(nights, dtype=int))))
    print('last night is : ' + str(np.max(np.array(nights, dtype=int))))

    tables = []
    for i, night in enumerate(nights):
        print('Working on night ' + night + ' (' + str(i+1) + ' of ' +
              str(len(nights)) + ')')
        table = _concat(night=night, basedir=basedir, acq=acq,
                        user_basedir=user_basedir, workers=workers)
        if table is None:
            continue
        tables.append(table)

    if tables == []:
        return -1, -1, -1

    result = vstack(tables)

    sind = np.argsort(result['EXPID'] + 0.1*result['PETAL_LOC'])

    result = result[sind]

    print('assembling per-frame median extension with minimal quality cuts...')
    med = cube_index_median(result, matched_coadd=(not acq))
    print('assembling per-frame median extension with additional quality cuts...')

    _med = cube_index_median(result, extra_cuts=True, matched_coadd=(not acq))
    
    return result, med, _med

def _write_many_nights(night_min='20201214', night_max='99999999',
                       basedir=None, acq=False, phase='SV1',
                       outdir='.', user_basedir=None, workers=1):

    if not os.path.exists(outdir):
        print('output directory does not exists ... quitting')
        return

    if basedir is None:
        basedir = _get_default_basedir(acq=acq)

    result, med, _med = _concat_many_nights(night_min=night_min,
                                            night_max=night_max,
                                            basedir=basedir, acq=acq,
                                            user_basedir=user_basedir,
                                            workers=workers)

    cube_index = 0 if acq else -1
    if (np.sum(result['CUBE_INDEX'] != cube_index) > 0) or (np.sum(med['CUBE_INDEX'] != cube_index) > 0) or (np.sum(_med['CUBE_INDEX'] != cube_index) > 0):
        print('WARNING: wrong CUBE_INDEX detected')
    
    hdul = fits.HDUList(hdus=[fits.PrimaryHDU(),
                              fits.BinTableHDU(data=result),
                              fits.BinTableHDU(data=med),
                              fits.BinTableHDU(data=_med)])

    night = str(np.max(result['NIGHT']))
    flavor = 'matched_coadd' if not acq else 'acq'
    outname = 'offline_' + flavor + '_ccds_' + phase + '-thru_' + \
              night + '.fits'

    outname = os.path.join(outdir, outname)

    if os.path.exists(outname):
        print('summary file already exists !!!')
        return

    print('attempting to write multi-extension FITS output to ' + outname)
    hdul.writeto(outname)
    os.system('chgrp desi ' + outname)
    os.system('chmod ug-w ' + outname)
    os.system('chmod a+r ' + outname)
    print('done')

def _latest_summary_file(outdir='.'):
    
    files = glob.glob(outdir + '/offline_matched_coadd_ccds_SV3-thru_????????.fits')

    files = np.array([os.path.splitext(f)[0] for f in files])

    nights = np.array([int(f.split("_")[-1]) for f in files])    

    latest_summary_file = files[np.argmax(nights)]+'.fits'
    
    latest_night = str(nights[np.argmax(nights)])
    
    return latest_summary_file, latest_night

def _next_obsnight(night='00000000'):

    assert(len(night) == 8)

    d = datetime.date.fromisoformat(night[0:4]+'-'+night[4:6]+'-'+night[6:8])

    delta = datetime.timedelta(days=1)

    next_obsnight = datetime.date.isoformat(d+delta)

    next_obsnight = next_obsnight.replace("-","")
    
    return next_obsnight

def _load_concas_ccds_file(file):
    
    hdul=fits.open(file)

    result = Table(hdul[1].data)
    med = Table(hdul[2].data)
    _med = Table(hdul[3].data)

    return result, med, _med

def _append_many_nights(night_min='20201214', night_max='99999999',
                       basedir=None, acq=False, phase='SV1',
                       outdir='.', user_basedir=None, workers=1):

    if not os.path.exists(outdir):
        print('output directory does not exists ... quitting')
        return

    # Find the most recent daily summary file in the output dir
    latest_summary_file, latest_night = _latest_summary_file(outdir=outdir)
    print('The latest daily summary file includes data through '+latest_night)

    next_obsnight = _next_obsnight(night=latest_night)

    # Check for new data
    nights = _nights_list(night_min, night_max, basedir=basedir, acq=acq)

    if len(nights) == 0:
        print('Error: no reduced GFA data found in '+basedir)
        return
    
    if nights[-1] < next_obsnight:
        print('Quitting... No new GFA data on or after '+next_obsnight)
        return

    # Get the content of the most recent daily summary file
    print('Loading '+latest_summary_file)
    existingResult, existingMed, existing_Med = _load_concas_ccds_file(latest_summary_file)

    # Read in new data 

    if basedir is None:
        basedir = _get_default_basedir(acq=acq)

    result, med, _med = _concat_many_nights(night_min=next_obsnight,
                                            night_max=nights[-1],
                                            basedir=basedir, acq=acq,
                                            user_basedir=user_basedir,
                                            workers=workers)

    if result == -1:
        print('No new GFA image found!')
        return
    
    cube_index = 0 if acq else -1
    if (np.sum(result['CUBE_INDEX'] != cube_index) > 0) or (np.sum(med['CUBE_INDEX'] != cube_index) > 0) or (np.sum(_med['CUBE_INDEX'] != cube_index) > 0):
        print('WARNING: wrong CUBE_INDEX detected')

    # Combine old and new records
    new_result = vstack([existingResult, result])
    new_med = vstack([existingMed, med])
    new__med = vstack([existing_Med, _med])

    # save
    
    hdul = fits.HDUList(hdus=[fits.PrimaryHDU(),
                              fits.BinTableHDU(data=new_result),
                              fits.BinTableHDU(data=new_med),
                              fits.BinTableHDU(data=new__med)])

    night = str(np.max(result['NIGHT']))
    flavor = 'matched_coadd' if not acq else 'acq'
    outname = 'offline_' + flavor + '_ccds_' + phase + '-thru_' + \
              night + '.fits'

    outname = os.path.join(outdir, outname)

    if os.path.exists(outname):
        print('summary file already exists !!!')
        return

    print('attempting to write multi-extension FITS output to ' + outname)
    hdul.writeto(outname)
    os.system('chgrp desi ' + outname)
    os.system('chmod ug-w ' + outname)
    os.system('chmod a+r ' + outname)
    print('done')

if __name__=="__main__":
    descr = 'gather gfa_reduce _ccds table outputs'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('--acq', default=False, action='store_true',
                        help='gather acquisition image _ccds tables')

    parser.add_argument('--basedir', default=None, type=str,
                        help='input directory')

    parser.add_argument('--night_min', default='20201214', type=str,
                        help='first observing night')

    parser.add_argument('--night_max', default='99999999', type=str,
                        help='last observing night')

    parser.add_argument('--phase', default='SV1', type=str,
                        help='survey phase (SV1, SV2, ...)')

    parser.add_argument('--outdir', default='.', type=str,
                        help='directory in which to write output file')

    parser.add_argument('--my_redux_dir', default=None, type=str,
                        help='base directory for your gfa_reduce outputs')

    parser.add_argument('--workers', default=1, type=int,
                        help='number of concurrent read-in processes')

    args = parser.parse_args()

    if (args.workers < 1) or (args.workers > 32):
        print('bad number of workers specified')
    
    _write_many_nights(night_min=args.night_min, night_max=args.night_max,
                       basedir=args.basedir, acq=args.acq, phase=args.phase,
                       outdir=args.outdir, user_basedir=args.my_redux_dir,
                       workers=args.workers)

