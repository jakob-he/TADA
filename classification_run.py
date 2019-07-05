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
    parser = argparse.ArgumentParser('Full classification between pathogenic and non pathogenic variants with variable features.')
    parser.add_argument('-c1','--cnvs1', help='Non pathogenic pickeled CNV set.')
    parser.add_argument('-c2','--cnvs2', help='Pathogenic pickeled CNV set.')
    parser.add_argument('-f','--feature', help='Feature set. Options: \n basic_binary \n extended binary \n basic continuous \n extended continuous')
    parser.add_argument('-cls','--cls', help='Type of classifier. Allowed options are: \n lr = Logistic Regression')
    parser.add_argument('-kw','--kwargs',default={},help='Keyword arguments for the classifier. They have to be in the flowwing format "keyword_1=arg_1,keyword_2=arg_2"')
    parser.add_argument('-o','--output', help='Output location.')
    return parser.parse_args()


def run(args):
    # create kwarg dict
    kwargs = args.kwargs
    if kwargs:
        kwarg_split = args.kwargs.split(',')
        kwargs = {kwarg.split('=')[0]:int(kwarg.split('=')[1]) for kwarg in kwarg_split}

    # unpickle the annotated CNV files
    with pathlib.Path(args.cnvs1).open('rb') as non_pathogenic_cnvs:
        non_patho_cnvs = pickle.load(non_pathogenic_cnvs)

    with pathlib.Path(args.cnvs2).open('rb') as pathogenic_cnvs:
        patho_cnvs = pickle.load(pathogenic_cnvs)

    # create test and training set
    train_set, test_set, scaler, imputer = preprocessing.create_stratified_training_and_test_set(non_patho_cnvs,patho_cnvs,feature_type=args.feature,oneHot=False)
    lr = Classifier(classifier=args.cls, imputer = imputer,scaler=scaler,**kwargs)
    lr.feature_selection(train_set,args.feature,output_dir=args.output)
    lr.train(train_set,output_dir=args.output)
    lr.test(test_set,save=True,plot=True,output_dir=args.output)

    # plot roc curve
    #plotting.plot_multiple_roc([lr],[test_set],save=True,output=pathlib.Path(args.output) / 'ROC_curve.png')


def main():
    # set random seed
    np.random.seed(42)

    args = argparser()

    run(args)

if __name__ == '__main__':
    main()
