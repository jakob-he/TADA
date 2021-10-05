"""Classification of pathogenic and non pathogenic variants with variable settings based on input arguments"""
import argparse
import pathlib
import pickle
import numpy as np
import yaml
import pkg_resources

# own libraries
import tadasv.annotate_tads as annotate_tads
import tadasv.annotate_svs as annotate_svs
import tadasv.lib.preprocessing as preprocessing
import tadasv.lib.utils as utils

# plotting
import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd

def argparser():
    # parse inputs
    parser = argparse.ArgumentParser('Full classification between pathogenic and non pathogenic variants with variable features. Requires pretrained model. Run predict_variants -h for more details.')
    parser.add_argument('-c', '--config', default = 'config.yml', help='Path to the config file containing TAD,SV and annotation locations.')
    parser.add_argument('-d', '--default', action='store_true',
                        help='Use default settings and annptations. This requires the variant path to be set via -v!')
    parser.add_argument(
        '-v', '--variants', help='Path to the SV bed/vcf-file. Only usable in combination with the -d flag!')
    parser.add_argument(
        '-t', '--type', help='Type of variant (Either DEL or DUP). Only usable in combination with the -d flag!')
    parser.add_argument('-o','--output', default='./',help='Output location.')
    return parser.parse_args()


def predict(cfg, output):
    # check if pickled preannotated SVs files are available
    if any(cfg['SVS']['ANNOTATED'].values()):
        sv_dicts = {}
        for label, svs in cfg['SVS']['ANNOTATED'].items():
            with pathlib.Path(svs).open('rb') as sv_dict:
                sv_dicts[label] = pickle.load(sv_dict)
        labeled_sv_dicts = sv_dicts.items()
    else:
        # check if the preannotated TADs are available
        if cfg['TADS']['ANNOTATED']:
            with pathlib.Path(cfg['TADS']['ANNOTATED']).open('rb') as annotated_tads:
                tads = pickle.load(annotated_tads)
        else:
            tads = annotate_tads.annotate(cfg)
        labeled_sv_dicts = annotate_svs.annotate(cfg).items()

    # load model
    model = pickle.load(pathlib.Path(cfg['PRETRAINED_MODEL']).open('rb'))

    # get labels depending on the feature type
    if cfg['FEATURES'] == 'extended':
        feature_labels = ['Stop Codon', 'Start Codon', '5_UTR', '3_UTR', 'Gene Overlap (bp)', 'Enhancer Overlap (bp)', 'DDG2P Overlap (bp)', 'Number of affected Genes', 'Number of affected Enhancers', 'Boundary Distance',
                          'Boundary Stability', 'Gene Distance', 'Enhancer Distance', 'DDG2P Distance', 'Gene LOEUF', 'Enhancer conservation', 'Gene HI', 'CTCF Distance', 'HI LogOdds Score', 'Exon Overlap', 'MPOI']
    else:
        feature_labels = [f'{annotation} distance' for annotation in cfg['ANNOTATIONS']]

    for label, annotated_svs in labeled_sv_dicts:
        # get feature df
        feature_df = preprocessing.create_feature_df(annotated_svs, cfg['FEATURES'],feature_labels, csv=True)
        # drop genomic location for prediction
        predict_df = feature_df.drop(['CHR', 'START', 'END'],axis=1)

        # compute predictions
        y_pred = model.pipeline.predict(predict_df)
        if cfg['CLASSIFIER'] != "ossvm":
            y_pred_scores = model.pipeline.predict_proba(predict_df)[:,1]
            feature_df['Pathogenicity Score'] = y_pred_scores

        # save annotations and preditioncs as csv file
        feature_df['Pathogenicity Label'] = y_pred
        feature_df.to_csv(output /  f'Annotated_Predicted_{label}.csv',sep='\t',header=True,index=False)


def main():
    args = argparser()

    if args.default:
        # depending on the variant type get config file from package
        if args.type.upper() == 'DEL':
            cfg_stream = pathlib.Path(pkg_resources.resource_filename(
                __name__, 'data/config_del_default.yml'))
        elif args.type.upper() == 'DUP':
            cfg_stream = pathlib.Path(pkg_resources.resource_filename(
                __name__, 'data/config_dup_default.yml'))
        else:
            print(
                f'{args.type} is not supported. Only DEL and DUP are viable options!')
            return
        cfg = yaml.load(cfg_stream.open(), Loader=yaml.Loader)
        cfg['SVS']['RAW']['TEST'] = args.variants
        # iterate through cfg entires and set path to package path
        for key_1 in cfg:
            if key_1 == "PRETRAINED_MODEL":
                cfg[key_1] = cfg_stream.parent / cfg[key_1]
            elif key_1 not in ['SVS', 'FEATURES', 'KWARGS', 'CLASSIFIER']:
                for key_2 in cfg[key_1]:
                    cfg[key_1][key_2] = cfg_stream.parent / cfg[key_1][key_2]
    else:
        # read config file
        with pathlib.Path(args.config).open() as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.Loader)

    output = pathlib.Path(args.output)

    predict(cfg, output)

if __name__ == '__main__':
    main()
