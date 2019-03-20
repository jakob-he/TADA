"""Test loading Tad objects from a bed file"""
import unittest
import datetime

from lib.tad import Tad
from lib import utils

class TestTadLoading(unittest.TestCase):
    """Test case for loading Tad data"""

    def test_data_load(self):
        tad_beds = utils.objects_from_file("tests/test_data/test_tads.bed", 'TAD')
        self.assertEqual(len(tad_beds),4)
        self.assertEqual(tad_beds[0].start,800000)
        self.assertEqual(tad_beds[1].end,100000)
        self.assertEqual(tad_beds[2].chr,'chr3')
