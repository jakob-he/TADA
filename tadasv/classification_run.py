"""Classification of pathogenic and non pathogenic variants with variable settings based on input arguments"""
import argparse
import pathlib
import pickle
import numpy as np
import yaml

# own libraries
import tadasv.annotate_tads as annotate_tads
import tadasv.annotate_svs as annotate_svs
import tadasv.lib.plotting as plotting
import tadasv.lib.preprocessing as preprocessing
from tadasv.lib.classifier import Classifier


def argparser():
    # parse inputs
    parser = argparse.ArgumentParser(
        'Full classification between variant sets with variable features. Run classification_run -h for more details')
    parser.add_argument('-c', '--config', help='Config file containing the TAD, SV and annotation locations. The first SV file is considered to be class 0.')
    parser.add_argument('-fs', '--feature_selection', help='Enables the generation of correlation and PCA figures for feature seletion', action='store_true')
    parser.add_argument('-o', '--output', help='Output location.')
    parser.add_argument('-gridcv', '--grid_cv', help='Activate GridSearchCV to find the optimal hyper-parameters.', action='store_true')
    parser.add_argument('-rs', '--random_seed', help='Define a random seed for reproducability.')
    return parser.parse_args()


def run(args):
    # read config file
    with pathlib.Path(args.config).open() as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.Loader)

    # set random seed if available
    if args.random_seed:
        np.random.seed(int(args.random_seed))

    # check if pickled preannotated SVs files are available
    if any(cfg['SVS']['ANNOTATED'].values()):
        sv_dicts = []
        for label, svs in cfg['SVS']['ANNOTATED'].items():
            with pathlib.Path(svs).open('rb') as sv_dict:
                sv_dicts.append(pickle.load(sv_dict))
    else:
        # check if the preannotated TADs are available
        if cfg['TADS']['ANNOTATED']:
            with pathlib.Path(cfg['TADS']['ANNOTATED']).open('rb') as annotated_tads:
                tads = pickle.load(annotated_tads)
        else:
            tads = annotate_tads.annotate(cfg)
        sv_dicts = [svs for label,
                     svs in annotate_svs.annotate(cfg).items()]


    # get labels depending on the feature type
    if cfg['FEATURES'] == 'extended':
        feature_labels = ['PLS Overlap (bp)', 'dELS Overlap (bp)', 'pELS Overlap (bp)', 'Low-DNase Overlap (bp)', 'Stop Codon', 'Start Codon', '5_UTR', '3_UTR', 'Gene Overlap (bp)', 'Enhancer Overlap (bp)', 'DDG2P Overlap (bp)', 'Number of affected Genes', 'Number of affected Enhancers', 'Boundary Distance',
                          'Boundary Stability', 'Gene Distance', 'Enhancer Distance', 'DDG2P Distance', 'Gene LOEUF', 'Enhancer conservation', 'Gene HI', 'CTCF Overlap (bp)', 'HI LogOdds Score', 'Exon Overlap', 'MPOI']
    else:
        feature_labels = [f'{annotation} distance' for annotation in cfg['ANNOTATIONS']]

    # set one class true if a One-Class SVM is specified as classifier
    if cfg['CLASSIFIER'] == 'ossvm':
        one_class = True
    
    if one_class:
        train_set = preprocessing.create_stratified_training_and_test_set(
            sv_dicts[0], feature_type=cfg['FEATURES'], labels=feature_labels, one_class=True)
        # create a second feature df as testset if another set of SVs is present
        if len(sv_dicts) > 1:
            test_set = preprocessing.create_stratified_training_and_test_set(
                sv_dicts[1], feature_type=cfg['FEATURES'], labels=feature_labels, one_class=True)
    else: 
        # create test and training set stratified by the class distribution
        train_set, test_set = preprocessing.create_stratified_training_and_test_set(sv_dicts, feature_type=cfg['FEATURES'], labels=feature_labels)

    # create classifier object with the user-defined keywords
    cls = Classifier(classifier=cfg['CLASSIFIER'], **cfg['KWARGS'], one_class=one_class)

    # # create multiple figures for the purpose of feature slection, if needed
    # if args.feature_selection and one_class:
    #    cls.feature_selection(train_set, cfg['FEATURES'], output_dir=args.output)

    # #train the model
    cls.train(train_set, output_dir=args.output, permut_importance=args.feature_selection, gridcv=args.grid_cv)

    # test on the test split if its not a one-class classifier
    if not one_class:
        cls.test(test_set, save=False, plot=True, output_dir=args.output)
    # if multiple SV sets have been specified in the config file test the one class model on the second set
    elif len(sv_dicts) > 1:
        cls.test_one_class(test_set, train_set, feature_labels, plot=True, output_dir=args.output)


def main():
    # read cli
    args = argparser()

    run(args)


if __name__ == '__main__':
    main()
