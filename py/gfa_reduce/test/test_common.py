# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test gfa_reduce.common.
"""
import unittest
import os
import tempfile
from shutil import rmtree
from ..common import retrieve_git_rev


class TestCommon(unittest.TestCase):
    """Test gfa_reduce.common.
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

    def test_retrieve_git_rev(self):
        """Test git revision in a git checkout.
        """
        foo = retrieve_git_rev()
        self.assertRegex(foo, '[0-9a-f]+')

    def test_retrieve_git_rev_no_checkout(self):
        """Test git revision outside of a git checkout.
        """
        test_file = os.path.join(self.tmp, 'git_rev.txt')
        with open(test_file, 'w') as f:
            f.write('This is a test.')
        with self.assertRaises(RuntimeError) as e:
            foo = retrieve_git_rev(test_file)
        os.remove(test_file)
