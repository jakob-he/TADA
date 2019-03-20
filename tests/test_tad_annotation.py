"""Test the annotation of sample TADs"""
import unittest

from lib.tad import Tad
from lib.gene import Gene
from lib.enhancer import Enhancer
import lib.utils as utils

class TadAnnotationTest(unittest.TestCase):
    """Test class for the annotation of TADs"""

    def test_annotation(self):
        tad_beds = utils.objects_from_file('tests/test_data/test_tads.bed', 'TAD')
        enhancer_beds = utils.objects_from_file('tests/test_data/test_enhancer.bed', 'Enhancer',['ID'])
        gene_beds = utils.objects_from_file('tests/test_data/test_genes.bed', 'Gene',['name'])

        # create dict with chromsomes as keys
        gene_dict = utils.create_chr_dictionary_from_beds(gene_beds)
        enhancer_dict = utils.create_chr_dictionary_from_beds(enhancer_beds)
        tad_dict = utils.create_chr_dictionary_from_beds(tad_beds)

        #test dict creation
        self.assertEqual(len(tad_dict),4,'The conversion to dictionaries is not working!')

        # Annotate TADs with overlapping enhancer and genes
        annotated_tads = utils.create_annotated_tad_dict(tad_dict, gene_dict, enhancer_dict)

        self.assertEqual(annotated_tads['chr1'][0].count_genes(),2,f'TAD annotation with genes is not working!')
        self.assertEqual(annotated_tads['chr1'][0].count_enhancer(),2),'TAD annoation with enhancer is not working!'
