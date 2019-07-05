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
    parser = argparse.ArgumentParser('Full classification between pathogenic and non pathogenic variants with variable features.')
    parser.add_argument('-v','--variants', help='Bed formatted Variant file.')
    parser.add_argument('-t','--tads', help='Pickeled file with annotated TADs.')
    parser.add_argument('-m','--model', help='pickled scikit model which is supposed to be used for classification.')
    parser.add_argument('-o','--output', help='Output location.')
    parser.add_argument('-l','--labeled',action='store_true',help='True if variants are labeled with column "label".')
    parser.add_argument('-f','--feature', help='Feature set. Options: \n basic_binary \n extended binary \n basic continuous \n extended continuous')
    return parser.parse_args()


def run(args):
    #load annotated TAD data
    tads = pathlib.Path(args.tads)
    tads = pickle.load(tads.open('rb'))

    #load CNVs
    cnvs = utils.objects_from_file(args.variants,'CNV')

    #create cnv dict
    cnvs = utils.create_chr_dictionary_from_beds(cnvs)

    #annotate CNVS
    annotated_cnvs = utils.annotate_cnvs(tads,cnvs)

    #save raw CNV object as pickle file
    with open(pathlib.Path(args.output) / 'annotated_variants.p', 'wb') as output:
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

    # get scaled feature df
    feature_df = preprocessing.create_feature_df(annotated_cnvs, args.feature)

    # plot feature df
    # fig, axes = plt.subplots(ncols=feature_df.shape[1],figsize=(20,5))
    # for ax, col in zip(axes, feature_df.columns):
    #     sns.distplot(feature_df[col],ax=ax)
    # plt.tight_layout()
    # plt.show()

    if 'continuous' in args.feature:
        scaled_feature_df = preprocessing.scale_and_impute_df(feature_df,model.scaler, model.imputer)


    test_set = [scaled_feature_df, np.array(y)]
    y_pred_scores = model.test(test_set,save=True,plot=False,output_dir=args.output)


def main():
    args = argparser()

    run(args)

if __name__ == '__main__':
    main()
