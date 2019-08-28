"""Test the training of a classifier"""
import unittest
import pathlib
import pickle

import lib.utils as utils
import lib.preprocessing as preprocessing
from lib.classifier import Classifier



class ClassifierTrainingTest(unittest.TestCase):
    """Test class for the training of a classifier"""

    def test_classification(self):
        # unpickle the annotated CNV files
        with pathlib.Path('tests/test_data/test_cnvs_nonpatho.p').open('rb') as non_pathogenic_cnvs:
            non_patho_cnvs = pickle.load(non_pathogenic_cnvs)

        with pathlib.Path('tests/test_data/test_cnvs_patho.p').open('rb') as pathogenic_cnvs:
            patho_cnvs = pickle.load(pathogenic_cnvs)
        print(patho_cnvs)
        # create test and training set
        train_set, test_set, scaler, imputer = preprocessing.create_stratified_training_and_test_set(non_patho_cnvs,patho_cnvs,feature_type='basic_binary',oneHot=False)
        lr = Classifier(classifier='lr', imputer = imputer, scaler=scaler)
        lr.train(train_set,output_dir='tests/test_data',cv=False)
        lr.test(test_set,save=True,plot=False)
