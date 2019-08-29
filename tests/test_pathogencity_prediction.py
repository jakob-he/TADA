"""Test the pathogencity prediction of CNVs"""
import pathlib
import unittest
import pickle

import numpy as np

import lib.preprocessing as preprocessing
import lib.utils as utils

class TestPathogencityPrediction(unittest.TestCase):
    """Class to test pathogencity prediction"""

    def test_pathogencity_pred(self):
        # unpickle the annotated CNV files
        with pathlib.Path('tests/test_data/test_cnvs_nonpatho.p').open('rb') as non_pathogenic_cnvs:
            non_patho_cnvs = pickle.load(non_pathogenic_cnvs)

        # get feature df
        feature_df = preprocessing.create_feature_df(non_patho_cnvs, 'basic_binary')

        # load model
        model = pickle.load(pathlib.Path('tests/test_data/lr_model.p').open('rb'))

        y_pred_scores = model.clf.predict_proba(feature_df)[:, 1]

        # output CSV file with predictions
        feature_df = preprocessing.create_feature_df(non_patho_cnvs,'basic_binary',csv=True)
        feature_df['Predicted Pathogencity'] = y_pred_scores
        feature_df.to_csv('tests/test_data/predicted_nonpatho.csv',sep='\t',header=True,index=False)
