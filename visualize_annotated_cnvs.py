"""Visualize the annotated SVs"""

# standard libraries
import argparse
import pathlib
import numpy as np
import pandas as pd

# own libraries
from lib import plotting

# plotting
import matplotlib.pyplot as plt


def argparser():
    parser = argparse.ArgumentParser(description="Visualize annotated CNVs.")
    parser.add_argument('-c1', '--cnvs_1',
                        help='Path to the first CNV csv-file', required=True)
    parser.add_argument('-c2', '--cnvs_2',
                        help='Path to the second CNV csv-file', required=True)
    parser.add_argument('-o', '--output',
                        help='Output path for generated figure.', required=True)
    return parser.parse_args()


def main():
    # parse input arguments
    args = argparser()
    cnv_1_file_path = pathlib.Path(args.cnvs_1)
    cnv_2_file_path = pathlib.Path(args.cnvs_2)

    # read annotated CNV dataframes
    cnvs_1 = pd.read_csv(cnv_1_file_path,header=0,index_col=False,sep='\t')
    cnvs_2 = pd.read_csv(cnv_2_file_path,header=0,index_col=False,sep='\t')


    # add label columns to both dataframes
    cnvs_1['label'] = np.repeat(['Non Pathogenic'],cnvs_1.shape[0])
    cnvs_2['label'] = np.repeat(['Pathogenic'],cnvs_2.shape[0])

    # concatenate the two dataframes
    merged_cnvs = pd.concat([cnvs_1,cnvs_2])

    # plot distrubtion by label exlucding the location columns
    plotting.plot_feature_dist(merged_cnvs,exclude_features=['CHR','START','END'],by='label')

    # save fig to output location
    plt.savefig(pathlib.Path(args.output) / 'feature_dist.png',bbox_inches='tight')


if __name__ == "__main__":
    main()
