"""Test the annotation of CNVs"""
import unittest
import pickle
import lib.utils as utils


class CnvAnnotationTest(unittest.TestCase):
    """Test class for the annotation of CNVs"""

    def test_annotation(self):
        tad_beds = utils.objects_from_file('tests/test_data/test_tads.bed', 'TAD')
        enhancer_beds = utils.objects_from_file('tests/test_data/test_enhancer.bed', 'Bed',column_names=['ID'])
        gene_beds = utils.objects_from_file('tests/test_data/test_genes.bed', 'Bed',column_names=['name'])

        # create dict with chromsomes as keys
        gene_dict = utils.create_chr_dictionary_from_beds(gene_beds)
        enhancer_dict = utils.create_chr_dictionary_from_beds(enhancer_beds)
        tad_dict = utils.create_chr_dictionary_from_beds(tad_beds)

        # Annotate TADs with overlapping enhancer and genes
        annotated_tads = utils.create_annotated_tad_dict(tad_dict, {'enhancers':enhancer_dict,'genes':gene_dict})

        #load cnvs
        cnv_beds = utils.objects_from_file("tests/test_data/test_cnvs_nonpatho.bed", 'cnv')

        #create cnv dict
        cnv_dict = utils.create_chr_dictionary_from_beds(cnv_beds)

        #annotate cnvs
        annotated_cnvs = utils.annotate_cnvs(annotated_tads,cnv_dict)

        self.assertEqual(len(annotated_cnvs['chr2'][0].tads),1,'Annotation of TADs does not work!')
        self.assertEqual(len(annotated_cnvs['chr1'][0].tads[0].annotations['genes']),2,'Genes are not transferred to the CNV object!')
