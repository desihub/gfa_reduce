# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test gfa_reduce.dark_current.
"""
import unittest
import os
import tempfile
from shutil import rmtree
from .. import dark_current


class TestDarkCurrent(unittest.TestCase):
    """Test gfa_reduce.dark_current.
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
