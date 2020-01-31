"""Classification of pathogenic and non pathogenic variants with variable settings based on input arguments"""
import argparse
import pathlib
import pickle
import lib.plotting as plotting
import numpy as np

#own libraries
import lib.preprocessing as preprocessing
from lib.classifier import Classifier


def argparser():
    # parse inputs
    parser = argparse.ArgumentParser('Full classification between pathogenic and non pathogenic variants with variable features. Run classification_run -h for more details')
    parser.add_argument('-c1','--cnvs1', help='Non pathogenic pickeled CNV set.', required=True)
    parser.add_argument('-c2','--cnvs2', help='Pathogenic pickeled CNV set.', required=True)
    parser.add_argument('-f','--feature', help='Feature set. Options: \n basic_binary \n extended binary \n basic continuous \n extended continuous', required=True)
    parser.add_argument('-fs','--feature_selection',help='Enables the generation of correlation and PCA figures for feature seletion',action='store_true')
    parser.add_argument('-cls','--cls', help='Type of classifier. Allowed options are: \n lr = Logistic Regression', required=True)
    parser.add_argument('-kw','--kwargs',default={},help='Keyword arguments for the classifier. They have to be in the following format "keyword_1=arg_1,keyword_2=arg_2". Keywords that should be interpreteted as strings have to be in double quotes.')
    parser.add_argument('-o','--output', help='Output location.')
    parser.add_argument('-gridcv','--grid_cv', help='Activate GridSearchCV to find the optimal hyper-parameters.', action='store_true')
    parser.add_argument('-rs','--random_seed',help='Define a random seed for reproducability.')
    return parser.parse_args()


def run(args):
    # create kwarg dict
    kwargs = args.kwargs
    if kwargs:
        kwarg_split = args.kwargs.split(',')
        kwargs = {}
        for kwarg in kwarg_split:
            key, arg = kwarg.split('=')
            try:
                arg = int(arg)
            except:
                pass
            kwargs[key] = arg

    # set random seed if available
    if args.random_seed:
        np.random.seed(int(args.random_seed))

    # unpickle the annotated CNV files
    with pathlib.Path(args.cnvs1).open('rb') as non_pathogenic_cnvs:
        non_patho_cnvs = pickle.load(non_pathogenic_cnvs)

    with pathlib.Path(args.cnvs2).open('rb') as pathogenic_cnvs:
        patho_cnvs = pickle.load(pathogenic_cnvs)

    # create test and training set
    train_set, test_set = preprocessing.create_stratified_training_and_test_set(non_patho_cnvs,patho_cnvs,feature_type=args.feature,oneHot=False)
    lr = Classifier(classifier=args.cls, **kwargs)
    if args.feature_selection:
        lr.feature_selection(train_set,args.feature,output_dir=args.output)
    lr.train(train_set,output_dir=args.output,permut_importance=True,gridcv=args.grid_cv)
    lr.test(test_set,save=True,plot=True,output_dir=args.output)

    # plot roc curve
    #plotting.plot_multiple_roc([lr],[test_set],save=True,output=pathlib.Path(args.output) / 'ROC_curve.png')


def main():
    args = argparser()

    run(args)

if __name__ == '__main__':
    main()
