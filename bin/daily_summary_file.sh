#!/bin/bash
source /global/common/software/desi/desi_environment.sh main
module load gfa_reduce/main
python -u ${GFA_REDUCE}/py/gfa_reduce/scripts/gfa_recent_nights.py
python ${GFA_REDUCE}/py/gfa_reduce/scripts/run_concat_ccds_nersc.py
