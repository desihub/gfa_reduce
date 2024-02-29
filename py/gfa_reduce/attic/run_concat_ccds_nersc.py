import concat_ccds
import os

basedir = concat_ccds._get_default_basedir(acq=False)
outdir = os.environ['GFA_SUMMRY_FILE_DIR']
workers = 8
concat_ccds._append_many_nights(night_min='20210405', night_max='99999999', basedir=basedir, acq=False, phase='SV3', outdir=outdir, user_basedir=None, workers=workers)
