"""Classification of pathogenic and non pathogenic variants with variable settings based on input arguments"""
import argparse
import pathlib
import pickle
import numpy as np
import yaml

# own libraries
import lib.preprocessing as preprocessing
import lib.utils as utils

# plotting
import seaborn as sns
import matplotlib.pyplot as plt



def argparser():
    # parse inputs
    parser = argparse.ArgumentParser('Full classification between pathogenic and non pathogenic variants with variable features. Requires pretrained model. Run predict_variants -h for more details.')
    parser.add_argument('-c', '--config', default = 'config.yml', help='Path to the config file containing TAD,CNV and annotation locations.')
    parser.add_argument('-o','--output', default='./',help='Output location.')
    return parser.parse_args()


def predict(cfg, output):
    # check if pickled preannotated CNVs files are available
    if cfg['CNVS']['ANNOTATED']:
        cnv_dicts = {}
        for label, cnvs in cfg['CNVS']['ANNOTATED'].items():
            with pathlib.Path(cnvs).open('rb') as cnv_dict:
                cnv_dicts[label] = pickle.load(cnv_dict)
    else:
        # check if preannotated TADs are available
        if cfg['TADS']['ANNOTED']:
            tads = pickle.load(pathlib.Path(cfg['TADS']['ANNOTED']).open('rb'))
        else:
            tads = annotate_tads.annotate(cfg)
        cnv_dicts = annotate_cnvs.annotate(cfg)

    # load model
    model = pickle.load(pathlib.Path(cfg['PRETRAINED_MODEL']).open('rb'))

    # get labels depending on the feature type
    if cfg['FEATURES'] == 'extended':
        feature_labels = ['Number of affected Genes','Number of affected Enhancers','Boundary Distance', 'Boundary Stability', 'Gene Distance', 'Enhancer Distance', 'DDG2P Distance', 'Gene LOEUF','Enhancer conservation', 'Gene HI', 'CTCF Distance', 'HI LogOdds Score', 'Exon Overlap', 'MPOI']
    else:
        feature_labels = [f'{annotation} distance' for annotation in cfg['ANNOTATIONS']]

    for label, annotated_cnvs in cnv_dicts.items():
        # get feature df
        feature_df = preprocessing.create_feature_df(annotated_cnvs, cfg['FEATURES'],feature_labels, csv=True)

        # drop genomic location for prediction
        predict_df = feature_df.drop(['CHR', 'START', 'END'],axis=1)

        # compute predictions
        y_pred_scores = model.pipeline.predict(predict_df)

        # save annotations and preditioncs as csv file
        feature_df['Predicted Pathogenicity'] = y_pred_scores
        feature_df.to_csv(output /  f'Annotated_Predicted_{label}.csv',sep='\t',header=True,index=False)


def main():
    args = argparser()

    # read config file
    with pathlib.Path(args.config).open() as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.Loader)

    output = pathlib.Path(args.output)

    predict(cfg, output)

if __name__ == '__main__':
    main()
