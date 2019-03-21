"""Test loading CNV objects from a bed file"""
import unittest
import datetime

from lib.cnv import CNV
from lib import utils

class TestCnvLoading(unittest.TestCase):
    """Test case for loading CNV data"""

    def test_data_load(self):
        cnv_beds = utils.objects_from_file("tests/test_data/test_cnv.bed", 'cnv')
        self.assertEqual(len(cnv_beds),3,'Number of read bed entries is incorrect!')
        self.assertEqual(cnv_beds[0].data['SVTYPE'],'INS','Columns are not read properly!')
        #test if gzipped files can be loaded
        cnv_beds = utils.objects_from_file("tests/test_data/test_cnv.bed.gz", 'cnv')
