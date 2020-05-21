"""Test the training of a classifier"""
import unittest
import pathlib
import pickle
import yaml

import tada.lib.utils as utils
import tada.lib.preprocessing as preprocessing
from tada.lib.classifier import Classifier



class ClassifierTrainingTest(unittest.TestCase):
    """Test class for the training of a classifier"""

    def test_classification(self):
        utils.blockPrint()
        # read config file
        with pathlib.Path("tests/test_config_cls.yml").open() as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.Loader)

        output_dir = "tests/test_data"
        # unpickle the annotated CNV files
        with pathlib.Path(cfg['CNVS']['ANNOTATED']['TEST_NON_PATHOGENIC']).open('rb') as non_pathogenic_cnvs:
            non_patho_cnvs = pickle.load(non_pathogenic_cnvs)

        with pathlib.Path(cfg['CNVS']['ANNOTATED']['TEST_NON_PATHOGENIC']).open('rb') as pathogenic_cnvs:
            patho_cnvs = pickle.load(pathogenic_cnvs)

        # create test and training set
        feature_labels = ['Number of affected Genes','Number of affected Enhancers','Boundary Distance', 'Boundary Stability', 'Gene Distance', 'Enhancer Distance', 'DDG2P Distance', 'Gene LOEUF','Enhancer conservation', 'Gene HI', 'CTCF Distance', 'HI LogOdds Score', 'Exon Overlap', 'MPOI']
        train_set, test_set = preprocessing.create_stratified_training_and_test_set([non_patho_cnvs,patho_cnvs],feature_type='extended',labels=feature_labels)
        lr = Classifier(classifier=cfg['CLASSIFIER'])
        lr.train(train_set, output_dir=output_dir, permut_importance=False, gridcv=False)
        lr.test(test_set, save=False, plot=False, output_dir=output_dir)
        utils.enablePrint()
        self.assertEqual(lr.trained,True,'Annotation of TADs does not work!')
