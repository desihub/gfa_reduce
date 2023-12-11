import os
import concat_ccds
import time
from datetime import datetime

#outdir = '/n/home/datasystems/users/ameisner/GFA'
outdir = os.environ['GFA_SUMMRY_FILE_DIR']
workers = 6

basedir = concat_ccds._get_default_basedir(acq=False)

last_time = 10000
while True:
    time.sleep(60)
    now = datetime.now()

    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)

    current_time = current_time.replace(':', '')
    current_time = int(current_time[0:4])

    thresh = 700 # 7am KPNO time
    print(last_time, current_time)

    if (current_time >= thresh) and (last_time < thresh):
        concat_ccds._write_many_nights(night_min='20210405', night_max='99999999',
                                       basedir=basedir, acq=False, phase='SV3',
                                       outdir=outdir, user_basedir=None, workers=workers)

    last_time = current_time
