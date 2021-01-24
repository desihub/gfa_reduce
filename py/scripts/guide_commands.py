import os
import glob
import _util
from astropy.table import Table
import astropy.io.fits as fits
import numpy as np
import random

from gfa_reduce.common import expid_from_filename

basedir = '/exposures/desi'

def _raw_file_list(night):

    dir = os.path.join(basedir, night)
    assert(os.path.exists(dir))

    pattern = dir + '/????????/guide-????????.fits.fz'

    flist = glob.glob(pattern)

    # randomize order ?

    return flist

def _gather_metadata(flist):
    
    domshutu = []
    imagecam = []
    ncameras = []
    flavor = []
    night = []
    exptime = []
    program = []
    expid = []
    pmcover = []
    frames = []
    _flist = []

    for f in flist:
        print(f)    
        assert(os.path.exists(f))
        h = fits.getheader(f, extname='GUIDER')

        # HACK !!!!!!!!!!
        if 'DOMSHUTU' not in h:
            continue

        _flist.append(f)
        domshutu.append(h['DOMSHUTU'].strip())
        imagecam.append(h['GUIDECAM'].strip())
        flavor.append(h['FLAVOR'].strip().lower())
        ncameras.append(_util.n_guide_cameras(h['GUIDECAM'].strip()))
        night.append(str(h['NIGHT']))
        # I think GUIDTIME used to be EXPTIME in the pre-covid time period
        exptime.append(h['GUIDTIME'])
        program.append(h['PROGRAM'])
        expid.append(h['EXPID'])
        pmcover.append(h['PMCOVER'].strip())
        frames.append(h['FRAMES'])

    t = Table()
    t['expid'] = expid
    t['domshutu'] = domshutu
    t['guidecam'] = imagecam
    t['ncameras'] = ncameras
    t['flavor'] = flavor
    t['night'] = night
    t['fname_raw'] = _flist
    t['exptime'] = exptime
    t['program'] = program
    t['pmcover'] = pmcover
    t['frames'] = frames

    return t

def _apply_cuts(meta, no_domshutu_requirement=False):

    # meta should be a table

    assert(len(meta) > 0)

    good = np.ones(len(meta), dtype=bool)

    good = np.logical_and(good, meta['exptime'] > 0)
    good = np.logical_and(good, meta['flavor'] == 'science')
    good = np.logical_and(good, meta['ncameras'] > 0)
    good = np.logical_and(good, meta['pmcover'] == 'open')

    domshutu_good = np.ones(len(meta), dtype=bool)
    if not no_domshutu_requirement:
        good = np.logical_and(good, meta['domshutu'] == 'open')

    if np.sum(good) == 0:
        print('NO EXPOSURES PASS METADATA CUTS')
        assert(False)

    return meta[good]

def _one_guide_command(fname, night, cube_index, out_basedir='/n/home/datasystems/users/ameisner/reduced/v0021_guide', background=False):
    assert(os.path.exists(fname))
    assert(os.path.exists(out_basedir))

    expid = expid_from_filename(fname)

    outdir = os.path.join(out_basedir, night + '/' + str(expid).zfill(8))

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    cmd = 'nohup python -u /n/home/datasystems/users/ameisner/latest/gfa_reduce/py/gfa_reduce/gfa_red.py ' + fname + ' --outdir ' + outdir + ' --skip_image_outputs --cube_index ' + str(cube_index) + ' &> guide-' + str(expid).zfill(8) + '-' + str(cube_index).zfill(8) + '.log'

    return cmd

def _all_guide_commands(meta, out_basedir='/n/home/datasystems/users/ameisner/reduced/v0021_guide', background=False):
    # just loop over _one_guide_command

    assert(len(meta) > 0)

    cmds = []
    for row in meta:
        for cube_index in range(row['frames']):
            cmd = _one_guide_command(row['fname_raw'], row['night'], cube_index, out_basedir=out_basedir, background=background)
            print(cmd)
            cmds.append(cmd)

    return cmds

def _commands(night='20201214', out_basedir='/n/home/datasystems/users/ameisner/reduced/v0021_guide', no_domshutu_requirement=False, background=False, min_expid=-1, max_expid=10000000):

    assert(os.path.exists(out_basedir))

    night_outdir = os.path.join(out_basedir, night)
    if not os.path.exists(night_outdir):
        os.mkdir(night_outdir)

    flist = _raw_file_list(night)

    meta = _gather_metadata(flist)

    meta = meta[(meta['expid'] >= min_expid) & (meta['expid'] <= max_expid)]

    assert(len(meta) > 0)

    meta = _apply_cuts(meta, no_domshutu_requirement=no_domshutu_requirement)

    cmds = _all_guide_commands(meta, out_basedir=out_basedir, background=background)

    return cmds

def _launch_scripts(night, chunksize=207, min_expid=-1, max_expid=10000000, no_domshutu_requirement=False):

    # eventually propagate all keywords
    cmds = _commands(night=night, min_expid=min_expid, max_expid=max_expid, no_domshutu_requirement=no_domshutu_requirement)

    random.shuffle(cmds)

    n_scripts = int(np.ceil(float(len(cmds))/float(chunksize)))

    if n_scripts > 20:
         print('do not run more than 20 guide cube reductions in parallel !!!')
         assert(False)

    fnames = []
    for i in range(n_scripts):

        fname = 'guide_chunk_' + str(i).zfill(3) + '.sh'

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
    
    cmd_wait = 'sleep 10\n'.encode('ascii')

    with open('launch.sh', 'wb') as f:
        for i, fname in enumerate(fnames):
            cmd = 'nohup ./' + fname + ' &\n'
            f.write(cmd.encode('ascii'))
            if i != (len(fnames)-1):
                f.write(cmd_wait)

    f.close()

    for f in fnames:
        os.system('chmod a+rx ' + f)

    os.system('chmod a+rx launch.sh')
