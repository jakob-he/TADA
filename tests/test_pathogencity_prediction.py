"""Test the pathogencity prediction of CNVs"""
import io
import sys

import pathlib
import unittest
import pickle
import yaml

import numpy as np

import tada.lib.preprocessing as preprocessing
import tada.lib.utils as utils

class TestPathogencityPrediction(unittest.TestCase):
    """Class to test pathogencity prediction"""

    def test_pathogencity_prio(self):
        capturedOutput = io.StringIO()
        sys.stdout = capturedOutput

        # read config file
        with pathlib.Path("tests/test_config_pred.yml").open() as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.Loader)

        # unpickle the annotated CNV files
        with pathlib.Path(cfg['CNVS']['ANNOTATED']['TEST_NON_PATHOGENIC']).open('rb') as non_pathogenic_cnvs:
            non_patho_cnv_dict = pickle.load(non_pathogenic_cnvs)

        # get feature df
        feature_labels = ['Number of affected Genes','Number of affected Enhancers','Boundary Distance', 'Boundary Stability', 'Gene Distance', 'Enhancer Distance', 'DDG2P Distance', 'Gene LOEUF','Enhancer conservation', 'Gene HI', 'CTCF Distance', 'HI LogOdds Score', 'Exon Overlap', 'MPOI']
        feature_df = preprocessing.create_feature_df(non_patho_cnv_dict, cfg['FEATURES'],feature_labels,csv=True)

        # load model
        with pathlib.Path(cfg['PRETRAINED_MODEL']).open('rb') as saved_model:
            model = pickle.load(saved_model)

        # drop genomic location for prediction
        predict_df = feature_df.drop(['CHR', 'START', 'END'],axis=1)

        # predict pathogenicity
        y_pred_scores = model.pipeline.predict_proba(predict_df)[:, 1]

        sys.stdout = sys.__stdout__

        # output CSV file with predictions
        feature_df['Predicted Pathogenicity'] = y_pred_scores
        feature_df.to_csv('tests/test_data/predicted_nonpatho.csv',sep='\t',header=True,index=False)
        self.assertEqual(len(feature_df['Predicted Pathogenicity']),36,'Pathogenicity Prediction is not working!')
