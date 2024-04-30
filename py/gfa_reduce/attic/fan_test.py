import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import scoreatpercentile
import glob
from ci_reduce.common import expid_from_filename
import imageio

def get_overscan(im, amp):
    sh = im.shape
    assert(len(sh) == 2)
    assert(sh[0] == 1032)
    assert(sh[1] == 2248)

    if amp == 'E':
        xmin = 1074
        xmax = 1124
        ymin = 0
        ymax = 516
        
    if amp == 'F':
        xmin = 1124
        xmax = 1174
        ymin = 0
        ymax = 516
        
    if amp == 'G':
        xmin = 1124
        xmax = 1174
        ymin = 516
        ymax = 1032

    if amp == 'H':
        xmin = 1074
        xmax = 1124
        ymin = 516
        ymax = 1032

    result = im[ymin:ymax, xmin:xmax]
    shr = result.shape
    assert(shr[0] == 516)
    assert(shr[1] == 50)

    return result.astype('float64')

def retrieve_overscan_med(extname, amp):
    tab = fits.getdata('/global/homes/a/ameisner/gfa/pro/bias_medians_20191027.fits')

    mask = (tab['EXTNAME'] == extname) & (tab['AMP'] == amp)

    assert(np.sum(mask) == 1)

    return tab[mask]['OVERSCAN_MED'][0]

def plot_one_overscan(im, extname, amp):
    subim = get_overscan(im, amp)
    #medval = retrieve_overscan_med(extname, amp)
    medval = np.median(subim)
    
    print(subim.shape)
    print(medval)
    subim -= medval
    plt.imshow(subim, cmap='gray', interpolation='nearest', origin='lower',
               vmin=-30, vmax=30, aspect=1.0)
    plt.title(extname + '; ' + amp, fontsize=7)
    
    plt.gca().set_axis_off()
    # plt.show()
    
def _test_one_overscan():
    fname = '/project/projectdirs/desi/spectro/data/20191107/00025304/gfa-00025304.fits.fz'
    extname = 'GUIDE2'
    amp = 'E'
    im = fits.getdata(fname, extname=extname)
    plot_one_overscan(im, extname, amp)
    

def plot_one_exp(fname, expid, fans=True, amp='E', night=None):
    fig = plt.figure(figsize=(8, 8))


    plt.subplots_adjust(hspace=0.01, wspace=0.01)
    # focus1 missing from fan test data set probably due to 
    # michael schubnell's testing
    extnames = ['GUIDE0', 'FOCUS1', 'GUIDE2', 'GUIDE3', 'FOCUS4', 'GUIDE5',
                'FOCUS6', 'GUIDE7', 'GUIDE8', 'FOCUS9']
    fig.tight_layout()
    #fig.subplots_adjust(top=0.7)
    for i in range(len(extnames)):
        plt.subplot(1, 10, i+1)
        im = fits.getdata(fname, extname=extnames[i])
        plot_one_overscan(im, extnames[i], amp)

    if fans is not None:
        title = 'fans ' + ('on' if fans else 'off') + '; '
    else:
        if night is None:
            title = ''
        else:
            title = night + '; '

    plt.suptitle(title + 'expid = ' + str(expid) + ' ; amp ' + amp, y=0.93, color='r')
    #plt.savefig('blat.png', bbox_inches='tight', dpi=400)

def test_one_exp():
    fname = '/project/projectdirs/desi/spectro/data/20191107/00025304/gfa-00025304.fits.fz'
    plot_one_exp(fname)
    
def _loop(amp='E'):
    expids = np.arange(25304, 25364)

    for expid in expids:
        print('Working on: ', expid)
        fan_off = (25324 <= expid) and (25344 > expid)
        fname = '/project/projectdirs/desi/spectro/data/20191107/' + str(expid).zfill(8) + '/gfa-' + str(expid).zfill(8) + '.fits.fz'
        plot_one_exp(fname, expid, fans=(not fan_off), amp=amp)
        plt.savefig('fan_test/expid_' + str(expid) + '_amp' + amp + '.png',
                     dpi=300, bbox_inches='tight')
        plt.cla()

def _is_simulated(f):
    im = fits.getdata(f, extname='GUIDE0')
    sh = im.shape
    
    result = sh[1] == 2048
    return result

def _loop_night(night, amp='E', indstart=None):
    # get the list of file names

    night_dir = '/exposures/desi/' + night

    flist = glob.glob(night_dir + '/*/gfa*.fz')

    for i, f in enumerate(flist):
        if indstart is not None:
            if i < indstart:
                continue

        if _is_simulated(f):
            continue

        expid = expid_from_filename(f)
        plot_one_exp(f, expid, fans=None, amp=amp, night=night)
        plt.savefig('expid_' + str(expid) + '_amp' + amp + '-' + night + '.png',
                     dpi=300, bbox_inches='tight')
        plt.cla()
        
def _loop_control(amp='E'):
    expids = np.arange(24520, 24580)

    for expid in expids:
        print('Working on: ', expid)
        #fan_off = (25324 <= expid) and (25344 > expid)
        fan_off = False # HACK !!
        fname = '/project/projectdirs/desi/spectro/data/20191105/' + str(expid).zfill(8) + '/gfa-' + str(expid).zfill(8) + '.fits.fz'
        plot_one_exp(fname, expid, fans=(not fan_off), amp=amp)
        plt.savefig('fan_test_control/expid_' + str(expid) + '_amp' + amp +
                     '.png',
                     dpi=300, bbox_inches='tight')
        plt.cla()

def gif_one_amp(amp='E', control=False):

    png_dir = 'fan_test' + ('' if not control else '_control')
    flist = glob.glob(png_dir + '/*' + amp + '.png')
    flist.sort()
    images = []
    for i, f in enumerate(flist):
        print(i)
        images.append(imageio.imread(f))
    basename = 'fan_test' if not control else 'control'
    imageio.mimsave('/global/cscratch1/sd/ameisner/' + basename + '_amp' + amp + '.gif', images, duration=0.35)

def gif_one_amp_one_night(night, amp='E', lastn=30):
    flist = glob.glob('*amp' + amp + '*' + night + '.png')

    print(len(flist))
    flist.sort()

    if lastn is not None:
        flist = flist[(-1*lastn):]

    images = []
    for i, f in enumerate(flist):
        print(i)
        images.append(imageio.imread(f))

    outname = night + '_amp' + amp + '.gif'
    imageio.mimsave(outname, images, duration=0.35)
