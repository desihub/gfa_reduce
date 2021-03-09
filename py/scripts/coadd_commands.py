import os
import glob
from astropy.table import Table
import astropy.io.fits as fits
import numpy as np
import random

from gfa_reduce.common import expid_from_filename
import gfa_reduce.gfa_red as gfa_red

out_basedir = 'your_output_directory'

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
                       background=False, mjdrange=None, fieldmodel=False,
                       pmgstars=True):

    # assume that if mjdrange is not None, then it will be a two element list [mjdmin, mjdmax]

    assert(os.path.exists(fname))
    assert(os.path.exists(out_basedir))

    expid = expid_from_filename(fname)

    outdir = os.path.join(out_basedir, night + '/' + str(expid).zfill(8))

    cmd = 'python -u ' + gfa_red.__file__ + ' ' + fname + ' --outdir ' + outdir + ' --skip_image_outputs --cube_index -1 '

    if mjdrange is not None:
        _extra = '--mjdmin ' + str(mjdrange[0]) + ' --mjdmax ' + str(mjdrange[1]) + ' '
        cmd += _extra

    if fieldmodel:
        cmd += '--fieldmodel '

    if pmgstars:
        cmd += '--pmgstars '

    cmd += '&> coadd-' + str(expid).zfill(8) + '.log'

    return cmd


def _all_coadd_commands(flist, night, out_basedir=out_basedir,
                        background=False, match_spectro_mjd=True,
                        fieldmodel=False, pmgstars=True):
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

        cmd = _one_coadd_command(f, night, out_basedir=out_basedir, background=background, mjdrange=mjdrange, fieldmodel=fieldmodel, pmgstars=pmgstars)
        cmds.append(cmd)

    return cmds

def _commands(night='20201214', out_basedir=out_basedir, background=False,
              match_spectro_mjd=True, fieldmodel=False, pmgstars=True):

    flist_spectro = _spectro_list(night)

    flist = _guider_list(flist_spectro)

    cmds = _all_coadd_commands(flist, night, match_spectro_mjd=match_spectro_mjd, out_basedir=out_basedir, fieldmodel=fieldmodel, pmgstars=pmgstars)

    night_dir = os.path.join(out_basedir, night)

    if not os.path.exists(night_dir):
        print('attempting to make nightly subdirectory ' + \
              night_dir)
        os.mkdir(night_dir)

    return cmds

def _launch_scripts(night, chunksize=8, match_spectro_mjd=True,
                    out_basedir=out_basedir, fieldmodel=False,
                    pmgstars=True):

    if not os.path.exists(out_basedir):
        print('the base output directory ' + out_basedir + \
              ' does not exist, will attempt to create it now')
        # might be good to check that out_basedir is a one-element string
        os.mkdir(out_basedir)

    # eventually propagate all keywords
    cmds = _commands(night=night, match_spectro_mjd=match_spectro_mjd, out_basedir=out_basedir, fieldmodel=fieldmodel, pmgstars=pmgstars)

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
            print('writing chunk script ' + str(i+1) + ' ' +
                  ' of ' + str(n_scripts) + ' ' + fname)
       	    for cmd in cmds[indstart:indend]:
                f.write((cmd + '\n').encode('ascii'))

        f.close()

        fnames.append(fname)

    # create the launch script

    launch_name = 'launch.sh'
    with open(launch_name, 'wb') as f:
        print('writing overall launch script ' + launch_name)
        for fname in fnames:
            cmd = './' + fname + ' &\n'
            f.write(cmd.encode('ascii'))

    f.close()

    for f in fnames:
        os.system('chmod a+rx ' + f)

    os.system('chmod a+rx launch.sh')
