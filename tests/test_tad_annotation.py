"""Test the annotation of sample TADs"""
import unittest
import lib.utils as utils

class TadAnnotationTest(unittest.TestCase):
    """Test class for the annotation of TADs"""

    def test_annotation(self):
        tad_beds = utils.objects_from_file('tests/test_data/test_tads.bed', 'TAD')
        enhancer_beds = utils.objects_from_file('tests/test_data/test_enhancer.bed', 'Bed',['ID'])
        gene_beds = utils.objects_from_file('tests/test_data/test_genes.bed', 'Bed',['name'])

        # create dict with chromsomes as keys
        gene_dict = utils.create_chr_dictionary_from_beds(gene_beds)
        enhancer_dict = utils.create_chr_dictionary_from_beds(enhancer_beds)
        tad_dict = utils.create_chr_dictionary_from_beds(tad_beds)

        #test dict creation
        self.assertEqual(len(tad_dict),4,'The conversion to dictionaries is not working!')

        # Annotate TADs with overlapping enhancer and genes
        annotated_tads = utils.create_annotated_tad_dict(tad_dict, {'enhancers':enhancer_dict,'genes':gene_dict})

        self.assertEqual(len(annotated_tads['chr1'][0].annotations['genes']),2,'TAD annotation with genes is not working!')
        self.assertEqual(len(annotated_tads['chr1'][0].annotations['enhancers']),2),'TAD annoation with enhancer is not working!'
