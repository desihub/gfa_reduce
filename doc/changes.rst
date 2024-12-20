==========
Change Log
==========

*Note*: prior to ``0.30.0`` an alternate tag name convention was used. New
tags have been created to be consistent with ``X.Y.Z`` conventions, and the
corresponding old tag is listed with the date in the titles below.

1.0.4 (unreleased)
------------------

* No changes yet.

1.0.3 (2024-12-03)
------------------

* Add ``EXTNAME`` header keywords to summary file (PR `#27`_).

.. _`#27`: https://github.com/desihub/gfa_reduce/pull/27

1.0.2 (2024-10-23)
------------------

* Debug some processing failures by improving exception handling in
  ``gfa_reduce.scripts.gfa_single_night._run()`` (direct commit to ``main``, `e534bf3`_)
* Add a notebook to find instances of ``NaN`` in GFA summary files
  (direct commit to ``main``, `988ad02`_).
* Allow reprocessing of additional nights out-of-sequence (PR `#25`_).

.. _`e534bf3`: https://github.com/desihub/gfa_reduce/commit/e534bf3e82eb976ed8254f62fdf79f6ab41ad9ef
.. _`988ad02`: https://github.com/desihub/gfa_reduce/commit/988ad0263ff963ce0a3c66138e9b901f9f7abedc

.. _`#25`: https://github.com/desihub/gfa_reduce/pull/25

1.0.1 (2024-05-16)
------------------

* Fix an issue where summary files from a previous phase were not recognized (PR `#23`_).

.. _`#23`: https://github.com/desihub/gfa_reduce/pull/23

1.0.0 (2024-04-30)
------------------

First version for NERSC-only GFA processing, with :command:`desiInstall` support (PR `#20`_).

.. _`#20`: https://github.com/desihub/gfa_reduce/pull/20

0.30.0 (2024-01-09)
-------------------

Final reference tag before merging NERSC compatibility branch (PR `#18`_).

.. _`#18`: https://github.com/desihub/gfa_reduce/pull/18

0.27.0 (v0027, 2021-04-22)
--------------------------

v0027 set of reductions; fixed FIBER_FRACFLUX bug

0.26.0 (v0026, 2021-03-08)
--------------------------

v0026 full SV1 reprocessing (guider frames, matched coadds, acquisition images), including FIBERFAC, _ELG and _BGS quantities

0.22.0 (v0022, 2021-01-30)
--------------------------

version used for SV1 reprocessing labeled v0022

0.10.0 (20200530, 2020-05-30)
-----------------------------

before module rename

0.8.0 (v0008, 2020-03-19)
-------------------------

adapt to deal with a change to the raw guide cube image extension headers

0.7.0 (v0007, 2020-03-08)
-------------------------

changes to PSF-based FWHM fitting and handling of missing EXPTIME values

0.5.2 (v0005.02, 2020-02-12)
----------------------------

minor change to add another dummy extension to ``_catalog`` output

0.5.1 (v0005.01, 2020-02-10)
----------------------------

minor change to logging; same data model as v0005

0.5.0 (v0005, 2020-02-09)
-------------------------

version used for first attempt at full reprocessing of DESI GFA dataset

0.4.0 (v0004, 2020-02-03)
-------------------------

version used for v0004 processing of all ``gfa*.fits.fz`` images from start of commissioning through 20200201

0.3.0 (v0003_on_target, 2020-01-30)
-----------------------------------

used for initial batch of v0003 reductions for first week of on-target spectroscopy, 20200119-20200127

0.2.0 (center, 2019-11-12)
--------------------------

with Python implementation of center added

0.1.0 (pre_gfa, 2019-09-09)
---------------------------

saving what's here before the repo tranforms into a GFA reduction pipeline
