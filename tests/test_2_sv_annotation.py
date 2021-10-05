"""Test the annotation of SVs"""
import io
import sys

import unittest
import pickle
import pathlib
import yaml

import tadasv.lib.utils as utils
import tadasv.lib.preprocessing as preprocessing
from tadasv.annotate_svs import annotate



class CnvAnnotationTest(unittest.TestCase):
    """Test class for the annotation of SVs"""

    def test_annotation(self):
        capturedOutput = io.StringIO()
        sys.stdout = capturedOutput

        # read config file
        with pathlib.Path("tests/test_config_svs.yml").open() as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.Loader)

        output_dir = "tests/test_data"

        annotated_svs = annotate(cfg)

        #save annotated svs
        output_path = pathlib.Path(output_dir)
        for label, svs in annotated_svs.items():
            with open(output_path / f'Annotated_{label}.p', "wb") as output:
                pickle.dump(svs, output)
        feature_labels = ['Stop Codon', 'Start Codon', '5_UTR', '3_UTR', 'Gene Overlap (bp)', 'Enhancer Overlap (bp)', 'DDG2P Overlap (bp)', 'Number of affected Genes', 'Number of affected Enhancers', 'Boundary Distance',
                          'Boundary Stability', 'Gene Distance', 'Enhancer Distance', 'DDG2P Distance', 'Gene LOEUF', 'Enhancer conservation', 'Gene HI', 'CTCF Distance', 'HI LogOdds Score', 'Exon Overlap', 'MPOI']
        feature_df = preprocessing.create_feature_df(svs,cfg['FEATURES'],feature_labels,csv=True)
        feature_df.to_csv(output_path / f'Annotated_{label}.csv',sep='\t',header=True,index=False)
        sys.stdout = sys.__stdout__
        self.assertEqual(len(annotated_svs['TEST_PATHOGENIC']['chr2'][0].tads),1,'Annotation of TADs does not work!')
        self.assertEqual(len(annotated_svs['TEST_PATHOGENIC']['chr2'][0].tads[0].annotations['GENES']),17,'Genes are not transferred to the SV object!')
