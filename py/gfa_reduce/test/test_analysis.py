# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test gfa_reduce.analysis.
"""
import unittest
import os
import tempfile
from shutil import rmtree
from ..analysis import (amm_2dhist, asterisms, basic_catalog_stats,
                        basic_image_stats, center_contrast, djs_maskinterp,
                        djs_photcen, dm, gaussian, phot, radprof,
                        recalib_astrom, segment, sky, splinefwhm, util)


class TestAnalysis(unittest.TestCase):
    """Test gfa_reduce.analysis.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        rmtree(cls.tmp)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_something(self):
        pass
