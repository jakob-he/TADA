"""Classification of pathogenic and non pathogenic variants with variable settings based on input arguments"""
import argparse
import pathlib
import pickle
import numpy as np

# own libraries
import lib.preprocessing as preprocessing
import lib.utils as utils

# plotting
import seaborn as sns
import matplotlib.pyplot as plt



def argparser():
    # parse inputs
    parser = argparse.ArgumentParser('Full classification between pathogenic and non pathogenic variants with variable features. Run predict_variants -h for more details.')
    parser.add_argument('-c','--cnvs', help='BED or VCF file of CNVs.', required=True)
    parser.add_argument('-t','--tads', default='data/default_annotated_TADs.p',help='Pickeled file with annotated TADs.')
    parser.add_argument('-m','--model', default='data/default_model.p',help='pickled scikit model which is supposed to be used for classification.')
    parser.add_argument('-o','--output', default='./',help='Output location.')
    parser.add_argument('-l','--labeled',action='store_true',help='True if variants are labeled with column "label".')
    parser.add_argument('-f','--feature', default='extended_continuous', help='Feature set. Options: \n basic_binary \n extended binary \n basic continuous \n extended continuous')
    parser.add_argument('-csv','--csv', action='store_false', help='Return CSV file with the pathogencity odds and functional annotation of each CNV.')
    return parser.parse_args()


def run(args):
    #load annotated TAD data
    tads = pathlib.Path(args.tads)
    tads = pickle.load(tads.open('rb'))

    #load CNVs
    cnvs = utils.objects_from_file(args.cnvs,'CNV')

    #create cnv dict
    cnvs = utils.create_chr_dictionary_from_beds(cnvs)

    #annotate CNVS
    annotated_cnvs = utils.annotate_cnvs(tads,cnvs)

    #save raw CNV object as pickle file
    with open(pathlib.Path(args.output) / 'Annotated_CNVs.p', 'wb') as output:
        pickle.dump(annotated_cnvs, output)

    # load model
    model = pickle.load(pathlib.Path(args.model).open('rb'))

    # seperate labels and test data if label is in the cnvs
    y = []
    if args.labeled:
        for chrom in annotated_cnvs:
            for cnv in annotated_cnvs[chrom]:
                if cnv.tads:
                    y.append(int(cnv.data['LABEL']))

    # get feature df
    feature_df = preprocessing.create_feature_df(annotated_cnvs, args.feature)


    test_set = [scaled_feature_df, np.array(y)]
    y_pred_scores = model.pipeline.predict_proba(scaled_feature_df)[:, 1]

    if args.csv:
        feature_df = preprocessing.create_feature_df(annotated_cnvs,args.feature,csv=True)
        feature_df['Predicted Pathogenicity'] = y_pred_scores
        feature_df.to_csv(pathlib.Path(args.output) /  'Annotated_CNVs.csv',sep='\t',header=True,index=False)


def main():
    args = argparser()

    run(args)

if __name__ == '__main__':
    main()
