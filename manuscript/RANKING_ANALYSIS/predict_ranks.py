import argparse
import pathlib
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
import pickle

import tada.lib.preprocessing as preprocessing
from tada.lib.classifier import Classifier

def argparser():
    parser = argparse.ArgumentParser("Predict ranks of Variants for groups of 1:99 pathogenic to non-pathogenic.")
    parser.add_argument("-p","--pathogenic",help="Path to the pathogenic Variants.")
    parser.add_argument("-n","--non-pathogenic",help="Path to the non-pathogenic Variants.")
    #parser.add_argument("-m","--model",help="Path to the pickled model.")
    parser.add_argument("-o","--output",help="Path to the output directory.")
    return parser.parse_args()


def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx == len(array):
        return False, array[idx-1]
    elif idx > 0 and (math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return True, array[idx-1]
    else:
        return True, array[idx]

def get_rank(row, bins, bin_dict):
    # test if the variant is larger than the largest
    # non-pathogenic variant and there are 299 or more variants
    # in the same size bin
    in_bins, closest_bin = find_nearest(bins, row['SIZE'])
    if in_bins and len(bin_dict[closest_bin]) >= 99:
        # sample hundred non-pathogenic variants from the bin
        sampled_variants = bin_dict[closest_bin].sample(n=99)

        # merge sampled variants and the pathogenic variant
        merged_sampled_variants = sampled_variants.append(row, ignore_index=True)
        #predict_df = merged_sampled_variants.drop(['CHR', 'START', 'END', 'SIZE', 'LABEL'],axis=1)

        # predict pathogenicity
        #y_pred_scores = model.pipeline.predict_proba(predict_df)[:,1]
        #merged_sampled_variants['PATHOGENICITY'] = y_pred_scores
        # sort merged variants by pathogenicity and report the rank of the pathogenic variants
        #merged_sampled_variants.reset_index(inplace=True)
        merged_sampled_variants.sort_values(['Pathogenicity Score'], ascending=False, inplace=True, ignore_index=True)
        return merged_sampled_variants[merged_sampled_variants["LABEL"]==1].index[0]+1
    return 100

def main():
    # parse cli
    args = argparser()

    output_dir = pathlib.Path(args.output)
    pathogenic_path = pathlib.Path(args.pathogenic)
    non_pathogenic_path = pathlib.Path(args.non_pathogenic)
    #model_path = pathlib.Path(args.model)

    # load pretrained model
    #model = pickle.load(model_path.open('rb'))

    # read variant data as pandas DataFrame
    pathogenic_df = pd.read_csv(pathogenic_path,sep="\t",header=0,index_col=False)
    non_pathogenic_df = pd.read_csv(non_pathogenic_path,sep="\t",header=0,index_col=False)

    # add size to the Dataframes
    for idx, df in enumerate([non_pathogenic_df,pathogenic_df]):
        df['SIZE'] = np.log10(pd.to_numeric(df['END'])-pd.to_numeric(df['START']))

    # add labels to the DataFrame
    pathogenic_df['LABEL'] = 1
    non_pathogenic_df['LABEL'] = 0

    # bin non-pathogenic variants by size
    non_pathogenic_sizes = non_pathogenic_df['SIZE']
    uniques = np.sort(non_pathogenic_sizes.unique())
    bins = np.linspace(start=min(uniques),stop=max(uniques),num=60)
    bin_dict = {}
    for bin in bins:
        # get variants smaller than the bin size
        bin_variants = non_pathogenic_df[non_pathogenic_df['SIZE'] <= bin]
        bin_dict[bin] = bin_variants
        non_pathogenic_df = non_pathogenic_df[non_pathogenic_df['SIZE'] > bin]

    # pathogenic_sizes = [size for size in pathogenic_df['SIZE'] if size <= max(uniques)]
    # plt.Figure(figsize=(12,10))
    # plt.hist(non_pathogenic_sizes,bins=bins,color="b",alpha=0.8)
    # plt.hist(pathogenic_sizes,bins=bins,color="r",alpha=0.8)
    # plt.savefig('size_plot.png')


    # pick a pathogenic variant and get the
    # corresponding bin of non-pathogenic variants
    # subsequently predict the pathogencity
    bin_lists = []
    rank_bins = [10,20,40,60,80,100]
    for seed in np.random.choice(100, 30, replace=False):
        print(seed)
        np.random.seed(seed)
        ranks = [rank for rank in pathogenic_df.apply(get_rank,args=(bins,bin_dict),axis=1) if rank!=100]
        # bin ranks
        digitized = np.digitize(ranks,rank_bins,right=True)
        counts = np.unique(digitized,return_counts=True)
        bin_lists.append(counts[1])

    total_variants = np.sum(bin_lists[0])
    mean_rank_per = [(mean / total_variants)*100 for mean in np.mean(bin_lists,axis=0)]
    std_per = [(std / total_variants)*100 for std in np.std(bin_lists,axis=0)]
    print(total_variants)
    print(mean_rank_per)
    print(std_per)
    # plot ranks
    sns.set_style("white")
    sns.set_palette("muted")
    plt.Figure(figsize=(12,10))
    index = [0,1,2,3,4,5]
    sns.barplot(x=index,y=mean_rank_per, **{'yerr':std_per})
    sns.despine()
    plt.xticks(index, ['1-10','11-20','21-40','41-60','61-80','81-100'])
    plt.ylim(0,100)
    plt.ylabel('Percentage', fontsize=15)
    plt.xlabel('Ranks', fontsize=15)
    plt.title(f'n={total_variants}')
    plt.savefig('Ranking_Plot.png')



if __name__ == "__main__":
    main()
