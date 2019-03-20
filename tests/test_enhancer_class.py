"""Test loading Enhancer objects from a bed file"""
import unittest
import datetime

from lib.enhancer import Enhancer
from lib import utils

class TestEnhancerLoading(unittest.TestCase):
    """Test case for loading Enhancer data"""

    def test_data_load(self):
        enhancer_beds = utils.objects_from_file("tests/test_data/test_enhancer.bed", 'Enhancer', column_names=['ID'])
        self.assertEqual(len(enhancer_beds),3,'Number of read bed entries is incorrect!')
        self.assertEqual(enhancer_beds[0].start,800000,'The first bed entry is not loaded properly!')
        self.assertEqual(enhancer_beds[1].end,300000,'The second bed entry is not loaded properly!')
        self.assertEqual(enhancer_beds[1].data['ID'],'18',f'Additonal columns are not recognized!')
