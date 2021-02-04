import os
import glob
from astropy.table import Table
import astropy.io.fits as fits
import numpy as np
import random

from gfa_reduce.common import expid_from_filename

out_basedir = '/global/cscratch1/sd/ameisner/matched_coadd'

basedir = '/global/cfs/cdirs/desi/spectro/data'

def _spectro_list(night):
    # night should be a string

    dir = os.path.join(basedir, night)

    assert(os.path.exists(dir))

    pattern = dir + '/' + '????????' + '/desi-????????.fits.fz'

    flist = glob.glob(pattern)

    return flist

# this is the list of guide-????????.fits.fz files corresponding to a desi-????????.fits.fz spectro file
def _guider_list(spectro_flist):
    
    flist_pred = [s.replace('desi-', 'guide-') for s in spectro_flist]
    
    result = []

    for f in flist_pred:
        if os.path.exists(f):
            result.append(f)

    if len(result) == 0:
        print('no guide cubes with corresponding spectra ???')
        assert(False)

    return result

def _one_coadd_command(fname, night, out_basedir=out_basedir,
                       background=False, mjdrange=None, fieldmodel=True):

    # assume that if mjdrange is not None, then it will be a two element list [mjdmin, mjdmax]

    assert(os.path.exists(fname))
    assert(os.path.exists(out_basedir))

    expid = expid_from_filename(fname)

    outdir = os.path.join(out_basedir, night + '/' + str(expid).zfill(8))

    cmd = 'python -u /global/homes/a/ameisner/gfa_reduce/py/gfa_reduce/gfa_red.py ' + fname + ' --outdir ' + outdir + ' --skip_image_outputs --cube_index -1 '

    if mjdrange is not None:
        _extra = '--mjdmin ' + str(mjdrange[0]) + ' --mjdmax ' + str(mjdrange[1]) + ' '
        cmd += _extra

    if fieldmodel:
        cmd += '--fieldmodel '

    cmd += '&> coadd-' + str(expid).zfill(8) + '.log'

    return cmd


def _all_coadd_commands(flist, night, out_basedir=out_basedir,
                        background=False, match_spectro_mjd=False,
                        fieldmodel=True):
    # just loop over _one_coadd_command

    cmds = []
    for f in flist:
        if match_spectro_mjd:
            fname_spectro = f.replace('guide-', 'desi-')
            assert(os.path.exists(fname_spectro))
            h = fits.getheader(fname_spectro, extname='SPEC')
            mjdmin = h['MJD-OBS']
            mjdmax = h['MJD-OBS'] + h['EXPTIME']/(3600.0*24.0)
            mjdrange = [mjdmin, mjdmax]
        else:
            mjdrange = None

        cmd = _one_coadd_command(f, night, out_basedir=out_basedir, background=background, mjdrange=mjdrange, fieldmodel=fieldmodel)
        print(cmd)
        cmds.append(cmd)

    return cmds

def _commands(night='20201214', out_basedir=out_basedir, background=False,
              match_spectro_mjd=False, fieldmodel=True):

    flist_spectro = _spectro_list(night)

    flist = _guider_list(flist_spectro)

    cmds = _all_coadd_commands(flist, night, match_spectro_mjd=match_spectro_mjd, out_basedir=out_basedir, fieldmodel=fieldmodel)

    night_dir = os.path.join(out_basedir, night)

    if not os.path.exists(night_dir):
        os.mkdir(night_dir)

    return cmds

def _launch_scripts(night, chunksize=8, match_spectro_mjd=False,
                    out_basedir=out_basedir, fieldmodel=True):

    # hack for output base directory ...
    if match_spectro_mjd:
        out_basedir = out_basedir.replace('coadd_mode', 'matched_coadd')

    # eventually propagate all keywords
    cmds = _commands(night=night, match_spectro_mjd=match_spectro_mjd, out_basedir=out_basedir, fieldmodel=fieldmodel)

    random.seed(99)
    random.shuffle(cmds)

    n_scripts = int(np.ceil(float(len(cmds))/float(chunksize)))

    if n_scripts > 20:
         print('do not run more than 20 GFA reductions in parallel !!!')
         assert(False)

    fnames = []
    for i in range(n_scripts):

        fname = 'coadd_chunk_' + str(i).zfill(3) + '.sh'

        assert(not os.path.exists(fname))

        indstart = i*chunksize
        indend = min((i + 1)*chunksize, len(cmds))

        with open(fname, 'wb') as f:
       	    for cmd in cmds[indstart:indend]:
                print(cmd)
                f.write((cmd + '\n').encode('ascii'))

        f.close()

        fnames.append(fname)
        print('~'*80)

    # create the launch script
    
    with open('launch.sh', 'wb') as f:
        for fname in fnames:
            cmd = './' + fname + ' &\n'
            f.write(cmd.encode('ascii'))

    f.close()

    for f in fnames:
        os.system('chmod a+rx ' + f)

    os.system('chmod a+rx launch.sh')
