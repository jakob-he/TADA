"""Classification of pathogenic and non pathogenic variants with variable settings based on input arguments"""
import argparse
import pathlib
import pickle
import numpy as np
import yaml

# own libraries
import tada.annotate_tads as annotate_tads
import tada.annotate_cnvs as annotate_cnvs
import tada.lib.plotting as plotting
import tada.lib.preprocessing as preprocessing
from tada.lib.classifier import Classifier


def argparser():
    # parse inputs
    parser = argparse.ArgumentParser(
        'Full classification between variant sets with variable features. Run classification_run -h for more details')
    parser.add_argument('-c', '--config', help='Config file containing the TAD, CNV and annotation locations. The first CNV file is considered to be class 0.')
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

    # check if pickled preannotated CNVs files are available
    if any(cfg['CNVS']['ANNOTATED'].values()):
        cnv_dicts = []
        for label, cnvs in cfg['CNVS']['ANNOTATED'].items():
            with pathlib.Path(cnvs).open('rb') as cnv_dict:
                cnv_dicts.append(pickle.load(cnv_dict))
    else:
        # check if the preannotated TADs are available
        if cfg['TADS']['ANNOTATED']:
            with pathlib.Path(cfg['TADS']['ANNOTATED']).open('rb') as annotated_tads:
                tads = pickle.load(annotated_tads)
        else:
            tads = annotate_tads.annotate(cfg)
        cnv_dicts = [cnvs for label,
                     cnvs in annotate_cnvs.annotate(cfg).items()]


    # get labels depending on the feature type
    if cfg['FEATURES'] == 'extended':
        feature_labels = ['Number of affected Genes','Number of affected Enhancers','Boundary Distance', 'Boundary Stability', 'Gene Distance', 'Enhancer Distance', 'DDG2P Distance', 'Gene LOEUF','Enhancer conservation', 'Gene HI', 'CTCF Distance', 'HI LogOdds Score', 'Exon Overlap', 'MPOI']
    else:
        feature_labels = [f'{annotation} distance' for annotation in cfg['ANNOTATIONS']]

    # create test and training set stratified by the class distribution
    train_set, test_set = preprocessing.create_stratified_training_and_test_set(cnv_dicts, feature_type=cfg['FEATURES'], labels=feature_labels)

    # create classifier object with the user-defined keywords
    lr = Classifier(classifier=cfg['CLASSIFIER'], **cfg['KWARGS'])

    # create multiple figures for the purpose of feature slection, if needed
    if args.feature_selection:
       lr.feature_selection(train_set, cfg['FEATURES'], output_dir=args.output)

    #train a model and test in on a seperate test set
    lr.train(train_set, output_dir=args.output, permut_importance=args.feature_selection, gridcv=args.grid_cv)
    lr.test(test_set, save=False, plot=True, output_dir=args.output)

    # plot roc curve
    #plotting.plot_multiple_roc([lr],[test_set],save=True,output=pathlib.Path(args.output) / 'ROC_curve.png')


def main():
    # read cli
    args = argparser()

    run(args)


if __name__ == '__main__':
    main()
