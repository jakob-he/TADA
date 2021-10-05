"""Test the training of a classifier"""
import io
import sys

import unittest
import pathlib
import pickle
import yaml

import tadasv.lib.utils as utils
import tadasv.lib.preprocessing as preprocessing
from tadasv.lib.classifier import Classifier



class ClassifierTrainingTest(unittest.TestCase):
    """Test class for the training of a classifier"""

    def test_classification(self):
        capturedOutput = io.StringIO()
        sys.stdout = capturedOutput
        # read config file
        with pathlib.Path("tests/test_config_cls.yml").open() as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.Loader)

        output_dir = "tests/test_data"
        # unpickle the annotated SV files
        with pathlib.Path(cfg['SVS']['ANNOTATED']['TEST_NON_PATHOGENIC']).open('rb') as non_pathogenic_svs:
            non_patho_svs = pickle.load(non_pathogenic_svs)

        with pathlib.Path(cfg['SVS']['ANNOTATED']['TEST_NON_PATHOGENIC']).open('rb') as pathogenic_svs:
            patho_svs = pickle.load(pathogenic_svs)

        # create test and training set
        feature_labels = ['Stop Codon', 'Start Codon', '5_UTR', '3_UTR', 'Gene Overlap (bp)', 'Enhancer Overlap (bp)', 'DDG2P Overlap (bp)', 'Number of affected Genes', 'Number of affected Enhancers', 'Boundary Distance',
                          'Boundary Stability', 'Gene Distance', 'Enhancer Distance', 'DDG2P Distance', 'Gene LOEUF', 'Enhancer conservation', 'Gene HI', 'CTCF Distance', 'HI LogOdds Score', 'Exon Overlap', 'MPOI']
        train_set, test_set = preprocessing.create_stratified_training_and_test_set([non_patho_svs,patho_svs],feature_type='extended',labels=feature_labels)
        lr = Classifier(classifier=cfg['CLASSIFIER'])
        lr.train(train_set, output_dir=output_dir, permut_importance=False, gridcv=False)
        lr.test(test_set, save=False, plot=False, output_dir=output_dir)
        sys.stdout = sys.__stdout__
        self.assertEqual(lr.trained,True,'Annotation of TADs does not work!')
