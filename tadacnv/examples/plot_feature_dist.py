"""Load and annotate a set of CNVs"""
#standard libraries
import argparse
import pickle
import pathlib

#own libraries
import tada.lib.utils as utils
import tada.lib.preprocessing as preprocessing
import tada.lib.plotting as plotting

#third party libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def argparser():
    parser = argparse.ArgumentParser(description="Annotate a set of CNVs. Run annotate_cnvs -h for more details.")
    parser.add_argument('-c1', '--cnvs1',
                        help='Path to the first CNV file.',required=True)
    parser.add_argument('-c2', '--cnvs2', help='Path to the second CNV file.')
    return parser.parse_args()

def main():
    args = argparser()

    feature_df_1 = pd.read_csv(pathlib.Path(args.cnvs1),header=0,sep="\t")
    feature_df_2 = pd.read_csv(pathlib.Path(args.cnvs2),header=0,sep="\t")

    feature_df_1['Pathogenicity'] = np.repeat('Pathogenic',feature_df_1.shape[0])
    feature_df_2['Pathogenicity'] = np.repeat('Non-Pathogenic',feature_df_2.shape[0])

    merged_df = pd.concat([feature_df_1,feature_df_2])

    plotting.plot_feature_dist(merged_df, exclude_features=['CHR','START','END'], by='Pathogenicity')
    plt.savefig('Pathogenic_Non_Pathogenic_extended_cont_feature_dist.png',bbox_inches='tight')





if __name__ == "__main__":
    main()
