"""Test the annotation of sample TADs"""
import unittest
import pickle
import pathlib
import yaml

import tada.lib.utils as utils
import tada.lib.preprocessing as preprocessing
from tada.annotate_tads import annotate

class TadAnnotationTest(unittest.TestCase):
    """Test class for the annotation of TADs"""

    def test_annotation(self):
        utils.blockPrint()
        # read config file
        with pathlib.Path("tests/test_config_tads.yml").open() as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.Loader)

        output_dir = pathlib.Path("tests/test_data")

        annotated_tads = annotate(cfg)

        #save annotated tads
        with open(output_dir / 'Annotated_TADs.p', 'wb') as output:
            pickle.dump(annotated_tads, output)
        utils.enablePrint()
        self.assertEqual(len(annotated_tads['chr1'][0].annotations['GENES']),62,'TAD annotation is not working!')
        self.assertEqual(len(annotated_tads['chr1'][0].annotations['ENHANCERS']),111),'TAD annoation is not working!'
