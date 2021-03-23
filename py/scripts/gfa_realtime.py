#!/usr/bin/env python

import sys, os, time
import multiprocessing as mp
import argparse
import glob
from gfa_reduce.gfa_red import _proc
from gfa_reduce.common import expid_from_filename
import numpy as np
import astropy.io.fits as fits
import json
import gfa_reduce
import gfa_reduce.analysis.util as util

# using Stephen Bailey's "multirunner" template as the basis for this script
# https://raw.githubusercontent.com/sbailey/multirunner/master/multirunner.py

default_out_basedir = os.environ['DEFAULT_REDUX_DIR']
dts_raw = os.environ['DTS_RAW']

parser = argparse.ArgumentParser(usage = "{prog} [options]")
parser.add_argument("--night", type=str, default='now',
                    help="NIGHT string")
parser.add_argument("-n", "--numworkers", type=int,  default=1,
                    help="number of workers")
parser.add_argument("-w", "--waittime", type=int, default=5,
                    help="wait time between directory checks")
parser.add_argument("-e", "--expid_min", type=int, default=-1,
                    help="start with this EXPID value")
parser.add_argument("--out_basedir", type=str, default=default_out_basedir,
                    help="base output directory for GFA reductions")
parser.add_argument("--guider", default=False, action='store_true',
                    help="process guide-????????.fits.fz files instead of gfa-????????.fz files")
parser.add_argument("--focus", default=False, action='store_true',
                    help="optimize for focus scan analysis")
parser.add_argument("--indir", type=str, default=dts_raw,
                    help="base input directory to watch for new exposures")

args = parser.parse_args()

class ProcItem:
    def __init__(self, fname_raw, night, cube_index=None):
        self.fname_raw = fname_raw
        self.cube_index = cube_index
        self.is_guider_cube = cube_index is not None
        self.night = night

    def _cube_index_string(self):
        if self.is_guider_cube:
            return str(self.cube_index)
        else:
            return ''

def check_flavor_json(gfa_image_fname):
    gfa_json_fname = gfa_image_fname.replace('gfa-', 'request-')
    gfa_json_fname = gfa_json_fname.replace('.fits.fz', '.json')

    print(gfa_json_fname)
    assert(os.path.exists(gfa_json_fname))
 
    with open(gfa_json_fname) as json_file:
        data = json.load(json_file)

    return data['FLAVOR']

def is_flavor_science(gfa_image_fname):
    return check_flavor_json(gfa_image_fname).lower() == 'science'

print('running on host : ' + os.environ['HOSTNAME'])

print('PATH TO gfa_reduce IS : ')
print(gfa_reduce.__file__)

assert(os.path.exists(args.indir))

if args.night != 'now':
    night = args.night
else:
    night = util.get_obs_night_now()

print('OBSERVING NIGHT IS : ' + night)

def set_indir(night):
    global args

    indir = os.path.join(args.indir, night)

    print('WATCHING FOR NEW FILES IN : ' + indir)

    if not os.path.exists(indir):
        print('WARNING: INPUT DIRECTORY DOES NOT CURRENTLY EXIST')
    
    return indir

indir = set_indir(night)

guider = args.guider

print('BASE OUTPUT DIRECTORY : ')
print(args.out_basedir)
assert(os.path.exists(args.out_basedir))

out_basedir = args.out_basedir

def set_night_basedir_out(night):
    global out_basedir

    night_basedir_out = os.path.join(out_basedir, night)
    if not os.path.exists(night_basedir_out):
        os.mkdir(night_basedir_out)

    return night_basedir_out

night_basedir_out = set_night_basedir_out(night)

#- Create communication queue to pass files to workers
q = mp.Queue()

#- Function to run for each worker.
#- Listens on Queue q for filenames to process.
def run(workerid, q):
    global out_basedir
    print('Worker {} ready to go'.format(workerid))
    while True:
        image = q.get(block=True)
        filename = image.fname_raw
        print('Worker {} processing {}'.format(workerid, filename))
        sys.stdout.flush()
        #- Do something with that filename
        outdir = os.path.join(os.path.join(out_basedir, image.night),
                              str(expid_from_filename(filename)).zfill(8))

        if args.guider and (image.cube_index < 12):
            print('sleeping for 1 minute to avoid bad DTS links')
            time.sleep(60.0) # hack to deal with bad DTS links

        try:
            if not args.focus:
                _proc(filename, outdir=outdir, realtime=True,
                      cube_index=image.cube_index, skip_image_outputs=True,
                      skip_raw_imstats=True)
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

#- Start workers
for i in range(args.numworkers):
    p = mp.Process(target=run, args=(i, q))
    p.start()

#- Track what files have already been added to queue.
#- TODO: Upon startup, this could compare against files in output dir
#- and only load input files haven't already been processed.
exp_outdirs = glob.glob(night_basedir_out + '/????????')
prefix = 'guide-' if guider else 'gfa-'
# this will not work correctly for cases where some subset of slices of a guide cube have been processed...
known_files = set([indir + '/' + os.path.split(d)[-1] + '/' + prefix + \
                   os.path.split(d)[-1] + '.fits.fz' for d in exp_outdirs])

print('Number of known files = ', len(known_files))

#- Periodically check for any new files that may have appeared and add them
#- to the queue for a worker to process.

pattern = '????????/gfa*.fits.fz' if not guider else '????????/guide-????????.fits.fz'

while(True):

    if args.night == 'now':
        night_now = util.get_obs_night_now(verbose=True)
        if night_now != night:
            print('UPDATING OBSERVING NIGHT FROM ' + night + \
                  ' TO ' + night_now)
            night = night_now
            indir = set_indir(night)
            night_basedir_out = set_night_basedir_out(night)
    
    glob_pattern = os.path.join(indir, pattern)

    flist = glob.glob(glob_pattern)
    flist.sort()
    flist = np.array(flist)
    expids = np.array([expid_from_filename(f) for f in flist])
    flist = flist[expids >= args.expid_min]
    for filename in flist:
        if filename not in known_files:
            if guider or is_flavor_science(filename):
                print('Server putting {} in the queue'.format(filename))
                sys.stdout.flush()
                if not guider:
                    image = ProcItem(filename, night, cube_index=None)
                    q.put(image)
                else:
                    h = fits.getheader(filename, extname='GUIDER')
                    outdir = os.path.join(night_basedir_out,
                                          str(expid_from_filename(filename)).zfill(8))
                    os.mkdir(outdir) # avoids race condition in _proc ...
                    __cube_index = 1 if (h['FRAMES'] > 1) else 0
                    image = ProcItem(filename, night, cube_index=__cube_index)
                    q.put(image)
            else:
                print('skipping ' + filename + ' ; NOT flavor=science')
            known_files.add(filename)

    time.sleep(args.waittime)
