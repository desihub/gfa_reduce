import gen_reduce_commands
import argparse
import os
import astropy.io.fits as fits
import numpy as np

def ci_reduce_commands(data_dir, outdir, tab, backgrounded=False):
    flist = []
    for t in tab:
        flist.append(data_dir + '/' + t['NIGHT'].replace(' ', '') + '/' + str(t['EXPID']).zfill(8) + '/ci-' + str(t['EXPID']).zfill(8) + '.fits.fz')

    #print(flist)

    cmds = [gen_reduce_commands.ci_reduce_1command(f, outdir + '/' + tab[i]['NIGHT'].replace(' ', ''), backgrounded=backgrounded) for i, f in enumerate(flist)]

    return cmds

if __name__ == "__main__":
    descr = 'create Python commands to run ci_reduce pipeline from index file'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('data_dir', type=str, nargs=1)

    default_outdir = '/project/projectdirs/desi/users/ameisner/CI/reduced/v0001'

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

    default_index_name = os.path.join(os.environ['CI_REDUCE_ETC'], 
                                      'post_run_index.fits')

    parser.add_argument('--index_file', default=default_index_name, type=str,
                        help='name of index file to use')

    args = parser.parse_args()

    data_dir = args.data_dir[0]
    outdir = args.outdir
    backgrounded = args.backgrounded
    write = args.write

    tab = fits.getdata(args.index_file)

    tab = tab[((tab['FLAVOR'] == 'SCIENCE') | (tab['FLAVOR'] == 'science')) & (tab['NIGHT'].replace(' ', '') > '20190401')]

    print(len(tab))

    cmds = ci_reduce_commands(data_dir, outdir, tab, backgrounded=backgrounded)

    if not args.write:
        for cmd in cmds:
            print(cmd)
    else:
        chunksize = args.chunksize
        nchunks = int(np.ceil(float(len(cmds))/chunksize))

        script_names = []
        for i in range(nchunks):
            script_name = gen_reduce_commands.write_one_shell_script(i, cmds[i*chunksize:(i+1)*chunksize])
            script_names.append(script_name)

        if args.launch_script:
            gen_reduce_commands.write_launch_script(script_names)
