"""Test the annotation of CNVs"""
import io
import sys

import unittest
import pickle
import pathlib
import yaml

import tadacnv.lib.utils as utils
import tadacnv.lib.preprocessing as preprocessing
from tadacnv.annotate_cnvs import annotate



class CnvAnnotationTest(unittest.TestCase):
    """Test class for the annotation of CNVs"""

    def test_annotation(self):
        capturedOutput = io.StringIO()
        sys.stdout = capturedOutput

        # read config file
        with pathlib.Path("tests/test_config_cnvs.yml").open() as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.Loader)

        output_dir = "tests/test_data"

        annotated_cnvs = annotate(cfg)

        #save annotated cnvs
        output_path = pathlib.Path(output_dir)
        for label, cnvs in annotated_cnvs.items():
            with open(output_path / f'Annotated_{label}.p', "wb") as output:
                pickle.dump(cnvs, output)
        feature_labels = ['Number of affected Genes','Number of affected Enhancers','Boundary Distance', 'Boundary Stability', 'Gene Distance', 'Enhancer Distance', 'DDG2P Distance', 'Gene LOEUF','Enhancer conservation', 'Gene HI', 'CTCF Distance', 'HI LogOdds Score', 'Exon Overlap', 'MPOI']
        feature_df = preprocessing.create_feature_df(cnvs,cfg['FEATURES'],feature_labels,csv=True)
        feature_df.to_csv(output_path / f'Annotated_{label}.csv',sep='\t',header=True,index=False)
        sys.stdout = sys.__stdout__
        self.assertEqual(len(annotated_cnvs['TEST_PATHOGENIC']['chr2'][0].tads),1,'Annotation of TADs does not work!')
        self.assertEqual(len(annotated_cnvs['TEST_PATHOGENIC']['chr2'][0].tads[0].annotations['GENES']),17,'Genes are not transferred to the CNV object!')
