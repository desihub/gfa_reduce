#!/usr/bin/env python

import sys, os, time
from datetime import datetime
import multiprocessing as mp
import argparse
import glob
from gfa_reduce.gfa_red import _proc
from gfa_reduce.common import expid_from_filename
import numpy as np
import astropy.io.fits as fits
import json
import gfa_reduce


# revised from Aaron Meisner's gfa_realtime.py

default_out_basedir = os.environ['DEFAULT_REDUX_DIR']
# dts_raw = os.environ['DTS_RAW']
dts_raw = os.environ['DESI_SPECTRO_DATA']


class ProcItem:
    # need to add MJDMIN, MJDMAX here for matched coaddition case
    def __init__(self, fname_raw, night, cube_index=None, mjdmin=None, mjdmax=None):
        self.fname_raw = fname_raw
        self.cube_index = cube_index
        self.is_guider_cube = cube_index is not None
        self.night = night
        self.mjdmin = mjdmin
        self.mjdmax = mjdmax

    def _cube_index_string(self):
        if self.is_guider_cube:
            return str(self.cube_index)
        else:
            return ''

def _check_flavor_json(gfa_image_fname):
    gfa_json_fname = gfa_image_fname.replace('gfa-', 'request-')
    gfa_json_fname = gfa_json_fname.replace('.fits.fz', '.json')

    print(gfa_json_fname)
    assert(os.path.exists(gfa_json_fname))

    with open(gfa_json_fname) as json_file:
        data = json.load(json_file)

    return data['FLAVOR']

def _is_flavor_science(gfa_image_fname):
    return _check_flavor_json(gfa_image_fname).lower() == 'science'

def _set_indir(night, indir_base):
    indir = os.path.join(indir_base, night)

    print('SEARCHING FOR FILES IN : ' + indir)

    if not os.path.exists(indir):
        print('WARNING: INPUT DIRECTORY DOES NOT CURRENTLY EXIST')

    return indir

def _set_night_basedir_out(night, out_basedir):
    night_basedir_out = os.path.join(out_basedir, night)
    if not os.path.exists(night_basedir_out):
        os.mkdir(night_basedir_out)

    return night_basedir_out

#- Function to run for each worker.
#- Listens on Queue q for filenames to process.
def _run(workerid, q, out_basedir, focus):
    print('Worker {} ready to go'.format(workerid))
    # pause to allow the queue to be filled
    # time.sleep(5)
    while not q.empty():
       # print('Worker {} checking queue status'.format(workerid))
       # print('Is queue empty?'+str(q.empty()))
       # print('Items in the queue?'+str(q.qsize()))
       #  image = q.get(block=True)
        image = q.get_nowait()
        filename = image.fname_raw
        print('Worker {} processing {}'.format(workerid, filename))
        sys.stdout.flush()
        #- Do something with that filename
        outdir = os.path.join(os.path.join(out_basedir, image.night),
                              str(expid_from_filename(filename)).zfill(8))

        # if args.guider:
            ### print('sleeping for 1 minute to avoid bad DTS links')
            ### time.sleep(60.0) # hack to deal with bad DTS links

        try:
            if not focus:
                _proc(filename, outdir=outdir, realtime=True,
                      cube_index=image.cube_index, skip_image_outputs=True,
                      skip_raw_imstats=False, pmgstars=True, mjdmin=image.mjdmin, mjdmax=image.mjdmax)
            else:
                _proc(filename, outdir=outdir, realtime=True,
                      cube_index=image.cube_index, skip_image_outputs=True,
                      skip_raw_imstats=True, skip_astrometry=True,
                      no_ps1_xmatch=True, no_gaia_xmatch=True,
                      do_sky_mag=False, skip_2d_gaussians=True)
        except:
            print('PROCESSING FAILURE: ' + image.fname_raw + '   ' + \
                  image._cube_index_string())
        print('Worker {} done with {}'.format(workerid, filename))
        sys.stdout.flush()

    print('No files in the queue, Terminating worker {}'.format(workerid))


def _gfa_single_night(night='20210405', numworkers=8,
                      out_basedir=default_out_basedir, guider=True, focus=False,
                      indir=dts_raw):
    t0 = time.time()

    print('running on host : ' + os.environ['HOSTNAME'])
    print('PATH TO gfa_reduce IS : ')
    print(gfa_reduce.__file__)

    assert(os.path.exists(indir))

    print('OBSERVING NIGHT IS : ' + night)

    night_indir = _set_indir(night, indir)

    print('BASE OUTPUT DIRECTORY : ')
    print(out_basedir)
    assert(os.path.exists(out_basedir))

    night_basedir_out = _set_night_basedir_out(night, out_basedir)

#- Create communication queue to pass files to workers
    q = mp.Queue()


#- Track what files have already been added to queue.
#- TODO: Upon startup, this could compare against files in output dir
#- and only load input files haven't already been processed.
    exp_outdirs = glob.glob(night_basedir_out + '/????????')
    prefix = 'desi-' if guider else 'gfa-'
# this will not work correctly for cases where some subset of slices of a guide cube have been processed...
    known_files = set([night_indir + '/' + os.path.split(d)[-1] + '/' + prefix + os.path.split(d)[-1] + '.fits.fz' for d in exp_outdirs])

    print('Number of known files = ', len(known_files))


    pattern = '????????/gfa*.fits.fz' if not guider else '????????/desi-????????.fits.fz'


    glob_pattern = os.path.join(night_indir, pattern)

    flist = glob.glob(glob_pattern)
    flist.sort()
    flist = np.array(flist)
    expids = np.array([expid_from_filename(f) for f in flist])

    num_processed = 0
    for filename in flist:
        if filename not in known_files:
            if guider or _is_flavor_science(filename):
                print('Found file {} on the server'.format(filename))
                sys.stdout.flush()
                if not guider:
                    image = ProcItem(filename, night, cube_index=None)
                    q.put(image)
                    print('Server putting {} in the queue'.format(filename))
                    sys.stdout.flush()
                else:
                # should put in a pause here
                    fname_guide = filename.replace('desi-', 'guide-')
                    if not os.path.exists(fname_guide):
                        print('spectro file has no corresponding guide cube : ' + filename)
                        known_files.add(filename)
                        continue

                    h_spec = fits.getheader(filename, extname='SPEC')

                    mjdmin = h_spec['MJD-OBS']
                    mjdmax = h_spec['MJD-OBS'] + h_spec['EXPTIME']/(3600.0*24.0)

                    h = fits.getheader(fname_guide, extname='GUIDER')
                    outdir = os.path.join(night_basedir_out,
                                        str(expid_from_filename(fname_guide)).zfill(8))

                    if (h['FRAMES'] <= 1):
                        print('too few frames for matched coaddition : ' + fname_guide)
                        known_files.add(filename)
                        continue
                    __cube_index = -1
                    os.mkdir(outdir) # avoids race condition in _proc ...
                    image = ProcItem(fname_guide, night, cube_index=__cube_index, mjdmin=mjdmin, mjdmax=mjdmax)
                    q.put(image)
                    print('Server putting {} in the queue'.format(filename))
                    sys.stdout.flush()
                    num_processed = num_processed+1
            else:
                print('skipping ' + filename + ' ; NOT flavor=science')
            known_files.add(filename)

# wait for a few seconds to allow the queue to be filled
    time.sleep(5)

    procs = []

#- Start workers
    for i in range(numworkers):
        p = mp.Process(target=_run, args=(i, q, out_basedir, focus))
        procs.append(p)
        p.start()

# wait for the queue to be emptied
    while(not q.empty()):
        time.sleep(5)

# wait for all child-process to finish
    for proc in procs:
        proc.join()

    dt = time.time() - t0

    print('gfa_single_night took ' + '{:.2f}'.format(dt) + ' seconds')
    print('Number of parallel processes used: ', len(procs))
    print('Number of GFA guide images processed: ', num_processed)

    return num_processed


if __name__=="__main__":
    descr = 'reduce a single-night GFA guide data'
    parser = argparse.ArgumentParser(description=descr)

    parser = argparse.ArgumentParser(usage = "{prog} [options]")
    parser.add_argument("--night", type=str, default='20210405',
                        help="NIGHT string")
    parser.add_argument("-n", "--numworkers", type=int,  default=1,
                        help="number of workers")
    parser.add_argument("--out_basedir", type=str, default=default_out_basedir,
                        help="base output directory for GFA reductions")
    parser.add_argument("--guider", default=False, action='store_true',
                        help="process guide-????????.fits.fz files instead of gfa-????????.fz files")
    parser.add_argument("--focus", default=False, action='store_true',
                        help="optimize for focus scan analysis")
    parser.add_argument("--indir", type=str, default=dts_raw,
                        help="base input directory to watch for new exposures")

    args = parser.parse_args()

    if (args.numworkers < 1) or (args.numworkers > 32):
        print('bad numworkers specified')

    _gfa_single_night(night=args.night, numworkers=args.numworkers,
                      out_basedir=args.out_basedir, guider=args.guider, focus=args.focus,
                      indir=args.indir)
