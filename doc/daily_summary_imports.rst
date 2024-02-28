================
Internal imports
================

This is the internal import tree used by :command:`desi_gfa_reduce_daily_summary`.
This exercise was performed on 2024-02-27 to ensure that there were no
circular imports in this relatively complex package.

* :mod:`gfa_reduce.scripts.daily_summary` imports:

  - :mod:`gfa_reduce.scripts.concat_ccds`
  - :mod:`gfa_reduce.scripts.gfa_single_night` imports:

    * :mod:`gfa_reduce` (to get ``__file__``)
    * :mod:`gfa_reduce.common`
    * :mod:`gfa_reduce.gfa_red` imports:

      - :mod:`gfa_reduce.io` imports:

        * :mod:`gfa_reduce.image` imports:

          - :mod:`gfa_reduce.imred.dq_mask` imports:

            * :mod:`gfa_reduce.imred.load_calibs` imports:

              - :mod:`gfa_reduce.common`

            * :mod:`gfa_reduce.common`

          - :mod:`gfa_reduce.gfa_wcs` imports:

            * :mod:`gfa_reduce.common`

          - :mod:`gfa_reduce.analysis.util` imports:

            * :mod:`gfa_reduce.xmatch.gaia` imports:

              - :mod:`gfa_reduce.common`

            * :mod:`gfa_reduce.common`

          - :mod:`gfa_reduce.analysis.sky` imports:

            * :mod:`gfa_reduce.analysis.util` imports:

              - :mod:`gfa_reduce.xmatch.gaia` imports:

                * :mod:`gfa_reduce.common`

              - :mod:`gfa_reduce.common`

          - :mod:`gfa_reduce.analysis.segment` imports:

            * :mod:`gfa_reduce.analysis.util` imports:

              - :mod:`gfa_reduce.xmatch.gaia` imports:

                * :mod:`gfa_reduce.common`

              - :mod:`gfa_reduce.common`

            * :mod:`gfa_reduce.common`

          - :mod:`gfa_reduce.analysis.djs_photcen`
          - :mod:`gfa_reduce.analysis.radprof`
          - :mod:`gfa_reduce.analysis.splinefwhm`
          - :mod:`gfa_reduce.common`

        * :mod:`gfa_reduce.exposure` imports:

          - :mod:`gfa_reduce.analysis.phot` imports:

            * :mod:`gfa_reduce.analysis.util` imports:

              - :mod:`gfa_reduce.xmatch.gaia` imports:

                * :mod:`gfa_reduce.common`

              - :mod:`gfa_reduce.common`

            * :mod:`gfa_reduce.analysis.gaussian`
            * :mod:`gfa_reduce.analysis.djs_maskinterp`
            * :mod:`gfa_reduce.analysis.djs_photcen`
            * :mod:`gfa.reduce.common`

          - :mod:`gfa_reduce.analysis.util` imports:

            * :mod:`gfa_reduce.xmatch.gaia` imports:

              - :mod:`gfa_reduce.common`

            * :mod:`gfa_reduce.common`

          - :mod:`gfa_reduce.imred.load_calibs` imports:

            * :mod:`gfa_reduce.common`

          - :mod:`gfa_reduce.dark_current` imports:

            * :mod:`gfa_reduce.common`

          - :mod:`gfa_reduce.common`

        * :mod:`gfa_reduce.analysis.basic_catalog_stats` imports:

          - :mod:`gfa_reduce.analysis.util` imports:

            * :mod:`gfa_reduce.xmatch.gaia` imports:

              - :mod:`gfa_reduce.common`

            * :mod:`gfa_reduce.common`

        * :mod:`gfa_reduce.analysis.basic_image_stats` imports:

          - :mod:`gfa_reduce.analysis.util` imports:

            * :mod:`gfa_reduce.xmatch.gaia` imports:

              - :mod:`gfa_reduce.common`

            * :mod:`gfa_reduce.common`

        * :mod:`gfa_reduce.analysis.util` imports:

          - :mod:`gfa_reduce.xmatch.gaia` imports:

            * :mod:`gfa_reduce.common`

          - :mod:`gfa_reduce.common`

        * :mod:`gfa_reduce.common`

        * :mod:`gfa_reduce.xmatch.gaia` imports:

          - :mod:`gfa_reduce.common` imports:

        * :mod:`gfa_reduce.gfa_wcs` imports:

          - :mod:`gfa_reduce.common`

        * :mod:`gfa_reduce.common`

      - :mod:`gfa_reduce.analysis.util` imports:

        * :mod:`gfa_reduce.xmatch.gaia` imports:

          - :mod:`gfa_reduce.common`

        * :mod:`gfa_reduce.common`

      - :mod:`gfa_reduce.analysis.recalib_astrom` imports:

        * :mod:`gfa_reduce.analysis.asterisms` imports:

          - :mod:`gfa_reduce.analysis.amm_2dhist`
          - :mod:`gfa_reduce.analysis.center_contrast`
          - :mod:`gfa_reduce.common`
          - :mod:`gfa_reduce.gfa_wcs` imports:

            * :mod:`gfa_reduce.common`

          - :mod:`gfa_reduce.xmatch.gaia` imports:

            * :mod:`gfa_reduce.common`

      - :mod:`gfa_reduce.analysis.dm` imports:

        * :mod:`gfa_reduce.common`

      - :mod:`gfa_reduce.common`
