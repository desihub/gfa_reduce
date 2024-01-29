import concat_ccds
import gfa_single_night
import os, time
import numpy as np
import gfa_reduce.analysis.util as util

import argparse

basedir = os.environ['DTS_RAW']
out_basedir = os.environ['DEFAULT_REDUX_DIR']
numworkers = 8


def _gfa_recent_nights(verbose=False):

    nights = concat_ccds._nights_list('20210405', '99999999', basedir=basedir)
    processed_nights = concat_ccds._nights_list('20210405', '99999999', basedir=out_basedir)

    if (verbose==True):
        print('Most recent night with raw GFA data: '+max(nights))
        print('Most recent night with processed GFA data: '+max(processed_nights))
    
    if(max(processed_nights) >= max(nights)):
        print('No new raw GFA data to process')
        return

    _nights = np.array(nights)
    _nights = _nights[(_nights > max(processed_nights))]

    nights = _nights.tolist()

    nights.sort()

    t0 = time.time()
    num_processed = 0
    for night in nights:
        _n = gfa_single_night._gfa_single_night(night=night, numworkers=numworkers, out_basedir=out_basedir, guider=True, focus=False, indir=basedir)
        num_processed = num_processed + _n
        
    dt = time.time() - t0
    
    print('gfa_recent_nights took ' + '{:.2f}'.format(dt) + ' seconds')
    print('Number of GFA guide images processed: ', num_processed)

if __name__=="__main__":
    descr = 'reduce all GFA guide images that have not been processed yet'
    parser = argparse.ArgumentParser(description=descr)

    parser = argparse.ArgumentParser(usage = "{prog} [options]")
    parser.add_argument("--verbose", default=False, action='store_true',
                        help="verbose")

    args = parser.parse_args()
    
    _gfa_recent_nights(verbose=args.verbose)
