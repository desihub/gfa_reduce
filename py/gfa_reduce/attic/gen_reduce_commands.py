#!/usr/bin/env python

import argparse
import glob
import os
import numpy as np

def ci_reduce_1command(fname_in, outdir, backgrounded=False):
    basename = os.path.basename(fname_in)

    pos = basename.find('.fits')

    assert(pos != -1)

    expname = basename[0:pos]
    this_outdir = os.path.join(outdir, expname)

    cmd = 'python -u /global/homes/a/ameisner/ci_reduce/py/ci_reduce/ci_proc.py ' + fname_in + ' --outdir ' + this_outdir + ' &> ' + expname + '.log'

    if backgrounded:
        cmd += ' &'

    return cmd

def ci_reduce_commands(data_dir, outdir, backgrounded=False):
    flist = glob.glob(data_dir + '/*/ci-????????.fits.fz')

    cmds = [ci_reduce_1command(f, outdir, backgrounded=backgrounded) for f in flist]

    return cmds

def write_one_shell_script(n, cmds):
    outname = 'r' + str(n).zfill(3) + '.sh'

    assert(not os.path.exists(outname))

    f = open(outname, 'wb')

    for cmd in cmds:
        f.write((cmd + '\n').encode())

    f.close()

    return outname

def write_launch_script(script_names):
    outname = 'launch.sh'

    assert(not os.path.exists(outname))

    f = open(outname, 'wb')

    for script_name in script_names:
        cmd = '. ' + script_name + ' &'
        f.write((cmd + '\n').encode())

    f.close()

if __name__ == "__main__":
    descr = 'create Python commands to run ci_reduce pipeline'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('data_dir', type=str, nargs=1)

    default_outdir = '/project/projectdirs/desi/users/ameisner/CI/ci_data_challenge/results'

    parser.add_argument('--outdir', default=default_outdir, type=str,
                        help='directory to write outputs in')

    parser.add_argument('--backgrounded', default=False, action='store_true', 
                        help='add ampersand to the end of each command')

    parser.add_argument('--write', default=False, action='store_true',
                        help='write commands to output shell script(s)')

    parser.add_argument('--chunksize', type=int, default=65,
                        help='commands per shell script')

    parser.add_argument('--launch_script', default=False, action='store_true',
                        help='also write out launch script')

    args = parser.parse_args()

    data_dir = args.data_dir[0]
    outdir = args.outdir
    backgrounded = args.backgrounded
    write = args.write

    cmds = ci_reduce_commands(data_dir, outdir, backgrounded=backgrounded)

    if not args.write:
        for cmd in cmds:
            print(cmd)
    else:
        chunksize = args.chunksize
        nchunks = int(np.ceil(float(len(cmds))/chunksize))

        script_names = []
        for i in range(nchunks):
            script_name = write_one_shell_script(i, cmds[i*chunksize:(i+1)*chunksize])
            script_names.append(script_name)

        if args.launch_script:
            write_launch_script(script_names)
