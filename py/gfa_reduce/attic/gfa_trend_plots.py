#!/usr/bin/env python

"""
Multirunner example

Stephen Bailey
Lawrence Berkeley National Lab
March 2019
"""

import sys, os, time
import multiprocessing as mp
import argparse
import glob
import astropy.io.fits as fits
from astropy.table import Table, vstack
import numpy as np
import matplotlib.pyplot as plt
from ci_reduce.common import expid_from_filename
from ci_reduce.io import realtime_raw_read

parser = argparse.ArgumentParser(usage = "{prog} [options]")
parser.add_argument("-i", "--indir", type=str,  help="input directory")
parser.add_argument("-n", "--numworkers", type=int,  default=1, help="number of workers")
parser.add_argument("-w", "--waittime", type=int, default=5, help="wait time between directory checks")
args = parser.parse_args()

#- Create communication queue to pass files to workers
q = mp.Queue()

all_ccds = None

def _seeing_plot():
    global all_ccds
    expid_u = np.unique(all_ccds['EXPID'])

    seeing = []
    expids = []
    for expid in expid_u:
        keep = (all_ccds['EXPID'] == expid) & ((all_ccds['n_sources_for_shape'] >= 5))
        if np.sum(keep) == 0:
            continue
        expids.append(expid)
        seeing.append(np.median(all_ccds[keep]['fwhm_asec']))

    seeing = np.array(seeing)
    expids = np.array(expids)
    
    plt.scatter(expids, seeing)
    ymin = 0
    ymax = 2
    plt.scatter(expids[seeing > ymax], seeing[seeing > ymax]*0.0 + 1.95, marker='^')
    plt.xticks(rotation='vertical')
    plt.ylim((ymin, ymax))
    plt.xlabel('EXPID')
    plt.ylabel('FWHM (asec)')
    plt.savefig('seeing_plots/seeing-' + str(max(expid_u)).zfill(8) + '.png', dpi=200, bbox_inches='tight')
    plt.cla()
        
def _read_ccds_1exp(fname):
    global all_ccds
    print('READING ' + fname)
    hdul = realtime_raw_read(fname)
    ccds = Table(hdul[1].data)
    ccds['EXPID'] = expid_from_filename(fname)
    if all_ccds is None:
        all_ccds = ccds
    else:
        all_ccds = vstack([all_ccds, ccds])

    print('ccds length : ', len(all_ccds))
    _seeing_plot()
    
#- Function to run for each worker.
#- Listens on Queue q for filenames to process.
def run(workerid, q):
    print('Worker {} ready to go'.format(workerid))
    while True:
        filename = q.get(block=True)
        print('Worker {} processing {}'.format(workerid, filename))
        sys.stdout.flush()
        #- Do something with that filename; in this case just sleep
        _read_ccds_1exp(filename)
        print('Worker {} done with {}'.format(workerid, filename))
        sys.stdout.flush()

#- Start workers
for i in range(args.numworkers):
    p = mp.Process(target=run, args=(i, q))
    p.start()

#- Track what files have already been added to queue.
#- TODO: Upon startup, this could compare against files in output dir
#- and only load input files haven't already been processed.
known_files = set()

#- Periodically check for any new files that may have appeared and add them
#- to the queue for a worker to process.
glob_pattern = os.path.join(args.indir, '*/*_ccds.fits')
while(True):
    flist = glob.glob(glob_pattern)
    flist.sort()
    for filename in flist:
        if filename not in known_files:
            print('Server putting {} in the queue'.format(filename))
            sys.stdout.flush()
            q.put(filename)
            known_files.add(filename)

    time.sleep(args.waittime)
