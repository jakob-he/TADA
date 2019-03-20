"""Test loading Gene objects from a bed file"""
import unittest
import datetime

from lib.gene import Gene
from lib import utils

class TestGeneLoading(unittest.TestCase):
    """Test case for loading Gene data"""

    def test_data_load(self):
        gene_beds = utils.objects_from_file("tests/test_data/test_genes.bed", 'Gene', column_names=['name'])
        self.assertEqual(len(gene_beds),3,'Number of read bed entries is incorrect!')
        self.assertEqual(gene_beds[0].start,700000,'The first bed entry is not loaded properly!')
        self.assertEqual(gene_beds[1].end,3002000,'The second bed entry is not loaded properly!')
        self.assertEqual(gene_beds[1].data['name'],'GENE2',f'Additonal columns are not recognized!')
