#!/usr/bin/env python

import os
import glob
from astropy.table import Table
import astropy.io.fits as fits
import numpy as np
import random
import argparse

from gfa_reduce.common import expid_from_filename
import gfa_reduce.gfa_red as gfa_red

out_basedir = 'gfa_redux'

basedir = os.environ['DESI_SPECTRO_DATA']
# lostfound = os.environ['DESI_LOST_FOUND']
lostfound = os.path.join(os.environ['DESI_ROOT'], 'spectro', 'staging', 'lost+found')

def _search_dirs():
    return [basedir, lostfound]

def _spectro_list(night):
    # night should be a string

    search_dirs = _search_dirs()

    result = []
    for _basedir in search_dirs:
        dir = os.path.join(_basedir, night)

        if not os.path.exists(dir):
            continue

        pattern = dir + '/' + '????????' + '/desi-????????.fits.fz'

        flist = glob.glob(pattern)

        result = result + flist

    return result

# this is the list of guide-????????.fits.fz files corresponding to a
# desi-????????.fits.fz spectro file
def _guider_list(spectro_flist):

    flist_pred = [s.replace('desi-', 'guide-') for s in spectro_flist]

    result = []

    search_dirs = _search_dirs()

    spectro_flist_matched = []

    for i, f in enumerate(flist_pred):
        if os.path.exists(f):
            result.append(f)
            spectro_flist_matched.append(spectro_flist[i])
        elif os.path.exists(f.replace(search_dirs[0], search_dirs[1])):
            result.append(f.replace(search_dirs[0], search_dirs[1]))
            spectro_flist_matched.append(spectro_flist[i])
        elif os.path.exists(f.replace(search_dirs[1], search_dirs[0])):
            result.append(f.replace(search_dirs[1], search_dirs[0]))
            spectro_flist_matched.append(spectro_flist[i])

    if len(result) == 0:
        print('no guide cubes with corresponding spectra ???')

    assert(len(spectro_flist_matched) == len(result))

    return result, spectro_flist_matched

def _acq_list(night):
    # night should be a string

    search_dirs = _search_dirs()

    result = []
    for _basedir in search_dirs:
        dir = os.path.join(_basedir, night)

        if not os.path.exists(dir):
            continue

        pattern = dir + '/' + '????????' + '/guide-????????-0000.fits.fz'

        flist = glob.glob(pattern)

        result = result + flist

    return result

def _one_command(fname, night, out_basedir=out_basedir,
                 background=False, mjdrange=None, fieldmodel=False,
                 pmgstars=True, make_exp_outdir=True,
                 log_prefix='coadd', acq=False):

    # assume that if mjdrange is not None, then it will be a two element list
    # [mjdmin, mjdmax]

    assert(os.path.exists(fname))
    assert(os.path.exists(out_basedir))

    expid = expid_from_filename(fname)

    outdir = os.path.join(out_basedir, night + '/' + str(expid).zfill(8))

    if make_exp_outdir:
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    cmd = 'python -u ' + gfa_red.__file__ + ' ' + fname + ' --outdir ' + \
          outdir + ' --skip_image_outputs'

    if not acq:
        cmd += ' --cube_index -1 '
    else:
        cmd += ' '

    if mjdrange is not None:
        _extra = '--mjdmin ' + str(mjdrange[0]) + ' --mjdmax ' + \
                 str(mjdrange[1]) + ' '
        cmd += _extra

    if fieldmodel:
        cmd += '--fieldmodel '

    if pmgstars:
        cmd += '--pmgstars '

    cmd += '&> ' + log_prefix + '-' + str(expid).zfill(8) + '.log'

    return cmd

def _all_coadd_commands(flist, flist_spectro, night, out_basedir=out_basedir,
                        background=False, match_spectro_mjd=True,
                        fieldmodel=False, pmgstars=True,
                        make_exp_outdirs=True):
    # just loop over _one_command

    cmds = []
    for i, f in enumerate(flist):
        if match_spectro_mjd:
            fname_spectro = flist_spectro[i]
            assert(os.path.exists(fname_spectro))
            h = fits.getheader(fname_spectro, extname='SPEC')

            if (h['MJD-OBS'] is None) or (h['EXPTIME'] is None):
                continue

            mjdmin = h['MJD-OBS']
            mjdmax = h['MJD-OBS'] + h['EXPTIME']/(3600.0*24.0)
            mjdrange = [mjdmin, mjdmax]
        else:
            mjdrange = None

        cmd = _one_command(f, night, out_basedir=out_basedir,
                           background=background, mjdrange=mjdrange,
                           fieldmodel=fieldmodel, pmgstars=pmgstars,
                           make_exp_outdir=make_exp_outdirs)
        cmds.append(cmd)

    return cmds

def _gen_acq_commands(night='20201214', out_basedir=out_basedir,
                      background=False, fieldmodel=False,
                      make_exp_outdirs=True):

    flist_acq = _acq_list(night)

    if len(flist_acq) == 0:
        print('no acquisition images to process')
        return []

    cmds = []
    for f in flist_acq:
        cmd = _one_command(f, night, out_basedir=out_basedir,
                           background=background,
                           fieldmodel=fieldmodel, pmgstars=False,
                           make_exp_outdir=make_exp_outdirs,
                           log_prefix='acq', acq=True)
        cmds.append(cmd)

    return cmds

def _gen_coadd_commands(night='20201214', out_basedir=out_basedir,
                        background=False, match_spectro_mjd=True,
                        fieldmodel=False, pmgstars=True,
                        make_exp_outdirs=True):

    flist_spectro = _spectro_list(night)

    flist, flist_spectro_matched = _guider_list(flist_spectro)

    night_dir = os.path.join(out_basedir, night)

    if not os.path.exists(night_dir):
        print('attempting to make nightly subdirectory ' + \
              night_dir)
        os.mkdir(night_dir)

    cmds = _all_coadd_commands(flist, flist_spectro_matched, night,
                               match_spectro_mjd=match_spectro_mjd,
                               out_basedir=out_basedir, fieldmodel=fieldmodel,
                               pmgstars=pmgstars,
                               make_exp_outdirs=make_exp_outdirs)

    return cmds

def _launch_scripts(night, match_spectro_mjd=True, out_basedir=out_basedir,
                    fieldmodel=False, pmgstars=True, make_exp_outdirs=True):

    if not os.path.exists(out_basedir):
        print('the base output directory ' + out_basedir + \
              ' does not exist, will attempt to create it now')
        # might be good to check that out_basedir is a one-element string
        os.mkdir(out_basedir)

    cmds = _gen_coadd_commands(night=night, match_spectro_mjd=match_spectro_mjd,
                               out_basedir=out_basedir, fieldmodel=fieldmodel,
                               pmgstars=pmgstars,
                               make_exp_outdirs=make_exp_outdirs)

    print(str(len(cmds)) + ' matched coadd jobs to run')

    cmds_acq = _gen_acq_commands(night=night, out_basedir=out_basedir,
                                 fieldmodel=fieldmodel,
                                 make_exp_outdirs=make_exp_outdirs)

    print(str(len(cmds_acq)) + ' acquisition images to run')

    cmds = cmds + cmds_acq

    random.seed(99)
    random.shuffle(cmds)

    # number of concurrent processes to run on one Cori node
    # could consider pushing this higher
    n_scripts_max = 20

    chunksize = int(np.ceil(float(len(cmds))/float(n_scripts_max)))

    print('packaging ' + str(chunksize) + ' jobs per worker')

    n_scripts = int(np.ceil(float(len(cmds))/float(chunksize)))

    fnames = []
    for i in range(n_scripts):

        fname = 'chunk_' + str(i).zfill(3) + '_' + night + '.sh'

        assert(not os.path.exists(fname))

        indstart = i*chunksize
        indend = min((i + 1)*chunksize, len(cmds))

        with open(fname, 'wb') as f:
            print('writing chunk script ' + str(i+1) + ' ' +
                  ' of ' + str(n_scripts) + ' ' + fname)
            for cmd in cmds[indstart:indend]:
                f.write((cmd + '\n').encode('ascii'))

        f.close()

        fnames.append(fname)

    # create the launch script

    launch_name = 'launch_' + night + '.sh'
    with open(launch_name, 'wb') as f:
        print('writing overall launch script ' + launch_name)
        for fname in fnames:
            cmd = './' + fname + ' &\n'
            f.write(cmd.encode('ascii'))

    f.close()

    for f in fnames:
        os.system('chmod a+rx ' + f)

    os.system('chmod a+rx ' + launch_name)

if __name__ == "__main__":
    descr = 'generate nightly gfa_reduce processing launch scripts'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('night', type=str, nargs=1,
                        help="observing night")

    parser.add_argument('--out_basedir', default=out_basedir, type=str,
                        help='base output directory')

    parser.add_argument('--fieldmodel', default=False, action='store_true',
                        help='fit desimeter FieldModel')

    parser.add_argument('--skip_exp_outdirs', default=False,
                        action='store_true',
                        help="don't pre-generate per EXPID output directories")

    parser.add_argument('--skip_pmgstars', default=False, action='store_true',
                        help="skip PMGSTARS forced photometry")

    args = parser.parse_args()

    make_exp_outdirs = not args.skip_exp_outdirs
    pmgstars = not args.skip_pmgstars

    _launch_scripts(args.night[0], match_spectro_mjd=True,
                    out_basedir=args.out_basedir,
                    fieldmodel=args.fieldmodel, pmgstars=pmgstars,
                    make_exp_outdirs=make_exp_outdirs)
