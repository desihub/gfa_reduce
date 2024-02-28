# gfa_reduce

[![GitHub Actions CI Status](https://github.com/desihub/gfa_reduce/actions/workflows/python-package.yml/badge.svg?branch=desiInstall-compatibility)](https://github.com/desihub/gfa_reduce/actions)
[![Coverage Status](https://coveralls.io/repos/github/desihub/gfa_reduce/badge.svg?branch=desiInstall-compatibility)](https://coveralls.io/github/desihub/gfa_reduce?branch=desiInstall-compatibility)

Reduce and summarize GFA images.

## Directories

### py/gfa_reduce/scripts/

Scripts actively used for GFA reduction.

### py/gfa_reduce/attic/

Unused and unmaintained scripts.

### py/gfa_reduce/analysis/

Core image analysis such as source detection, centroid measurements, flux measurements.

### py/gfa_reduce/imred/

Utilities for converting raw images to reduced images.

### py/gfa_reduce/xmatch/

Cross-matching utilities, such as for matching to Gaia or other external catalogs.

## Top-level py/gfa_reduce Python files

### py/gfa_reduce/gfa_red.py

Primary driver that gets called to run the GFA reduction pipeline end to end.

### py/gfa_reduce/gfa_wcs.py

WCS-related utilities.

### py/gfa_reduce/common.py

Miscellaneous utilities.

### py/gfa_reduce/dark_current.py

Utilities related to dark current.

### py/gfa_reduce/exposure.py

Class that encapsulates a single GFA exposure consisting of multiple single-camera images.

### py/gfa_reduce/image.py

Class that encapsulates a single-camera image (and its metadata) drawn from a single GFA exposure.

### py/gfa_reduce/io.py

Input/output utilities.

## Convention for pixel coordinates

Unless otherwise stated in a particular portion of the code...

Throughout the codebase, pixel coordinate (0, 0) will be the **center** of the lower left pixel (where "lower left" applies in the context of the IDL/DS9 convention for 2D image display). "x" will refer to the "AXIS1" coordinate (long axis for GFA cameras, a.k.a. "width"), and "y" will refer to the "AXIS2" coordinate (short axis for GFA cameras, a.k.a. "height").

## Environment

### Details

* To access auxiliary calibration files, it will be necessary to set the `GFA_REDUCE_META` environment variable to their directory's location
  * examples of auxiliary calibration files include master biases, master flats, static bad pixel masks
  * these auxiliary files will not be checked in to git due to the file sizes involved
  * the authoritative copy of this directory can be found on NERSC at `/global/cfs/cdirs/desi/users/ameisner/GFA/gfa_reduce_meta`
* To enable Gaia cross-matching, the `GAIA_CAT_DIR` environment variable must be set to the full path of "chunks-gaia-dr2-astrom"
  * the recommended path at NERSC is `/global/cfs/cdirs/cosmo/work/gaia/chunks-gaia-dr2-astrom`
  * this directory contains a full-sky set of 12,288 FITS "chunk" files, one per nside = 32 HEALPix pixel, with ring-ordered indexing in equatorial coordinates
* This code is intended to be run at NERSC using the DESI software environment, which is Python 3 based:
  * https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/NERSC

### Example environment configuration script

The following is a bash script that can be used to configure one's environment for running the `gfa_reduce` package at NERSC:
```bash
source /global/cfs/cdirs/desi/software/desi_environment.sh master
export PYTHONPATH=${PYTHONPATH}:/global/homes/a/ameisner/gfa_reduce/py
export GFA_REDUCE_META=/global/cfs/cdirs/desi/users/ameisner/GFA/gfa_reduce_meta
export GAIA_CAT_DIR=/global/cfs/cdirs/cosmo/work/gaia/chunks-gaia-dr2-astrom
export PS_CAT_DIR=/global/cfs/cdirs/desi/users/ameisner/GFA/chunks-fitscat_qz_trim
```
Reading auxiliary files from other places or running a `gfa_reduce` checkout located elsewhere on NERSC would require modifications of these example paths.

it would be recommended to run `gfa_reduce` using a checkout of this repository's code rather than the exact NERSC location given above, as code development may be taking place within that `/global/homes/a/ameisner/gfa_reduce/py` directory

to obtain a checkout of the `gfa_reduce` code
```bash
git clone https://github.com/desihub/gfa_reduce.git
```
`gfa_reduce` also runs on the DESI cluster computers at Kitt Peak, although no configuration example for doing so is provided here

## Examples

### Basic examples of running the pipeline

for a `gfa*.fits.fz` GFA exposure
```
python -u /global/homes/a/ameisner/gfa_reduce/py/gfa_reduce/gfa_red.py /global/cfs/cdirs/desi/spectro/data/20191022/00020012/gfa-00020012.fits.fz --outdir 00020012 --dont_write_invvar --compress_reduced_image
```
for one frame of a `guide*.fits.fz` guider cube
```
python -u /global/homes/a/ameisner/gfa_reduce/py/gfa_reduce/gfa_red.py /global/cfs/cdirs/desi/spectro/data/20200124/00043860/guide-00043860.fits.fz --outdir 00043860 --dont_write_invvar --compress_reduced_image --cube_index 1
```
such `gfa_reduce` invocations print out various logging information, such as warnings and timings, which can be piped to log files for record keeping and/or debugging purposes

### Full help for running the pipeline
```
    gfa_reduce/py/gfa_reduce> python gfa_red.py --help
    usage: gfa_red.py [-h] [--outdir OUTDIR] [--careful_sky] [--no_cataloging]
                      [--no_gaia_xmatch] [--no_ps1_xmatch]
                      [--cube_index CUBE_INDEX] [--skip_image_outputs]
                      [--realtime] [--no_dark_rescaling] [--dont_write_invvar]
                      [--skip_psf_models] [--compress_reduced_image]
                      [--skip_raw_imstats] [--skip_astrometry] [--no_pm_pi_corr]
                      [--write_psf_cubes] [--write_detmap] [--write_full_detlist]
                      fname_in

    run the gfa_reduce pipeline on a GFA exposure

    positional arguments:
      fname_in

    optional arguments:
      -h, --help            show this help message and exit
      --outdir OUTDIR       directory to write outputs in
      --careful_sky         use image segmentation when deriving sky quantities
      --no_cataloging       reduce image without cataloging sources
      --no_gaia_xmatch      skip Gaia cross-match
      --no_ps1_xmatch       skip PS1 cross-match
      --cube_index CUBE_INDEX
                            guide cube index
      --skip_image_outputs  skip writing of full-frame image outputs
      --realtime            avoid crashing on partially written raw files
      --no_dark_rescaling   skip empirical rescaling of dark current
      --dont_write_invvar   don't write out invvar maps
      --skip_psf_models     skip generating per-camera PSF models
      --compress_reduced_image
                            compress reduced image output file
      --skip_raw_imstats    skip computing of raw image pixel statistics
      --skip_astrometry     skip astrometric recalibration
      --no_pm_pi_corr       do not correct Gaia positions for proper motion or parallax
      --write_psf_cubes     write image cubes of sources used to build PSF models
      --write_detmap        write detection map
      --write_full_detlist  write out the initial, full list of detections
      --max_cbox MAX_CBOX   maximum centroiding box size (pixels)
      --fieldmodel          fit and write desimeter field model
      --dont_write_catalog  don't write source catalog
      --dont_write_ccds     don't write CCDs table
      --multiproc           use multiprocessing to decrease wall time
      --skip_aper_phot      don't perform aperture photometry
      --det_sn_thresh DET_SN_THRESH
                            source detection significance threshold
      --skip_flatfield      skip flatfielding during pixel-level detrending
      --search_rad_arcmin SEARCH_RAD_ARCMIN
                            astrometric pattern match search radius (arcmin)
      --skip_sky_mag        skip sky surface brightness estimation
```

## Change Log

See `doc/changes.rst`.
