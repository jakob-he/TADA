"""Produces plots for multiple results files."""

# standard libraries
import argparse
import re
import pickle
import pathlib
import numpy as np

# own libararies
import lib.utils as utils
import lib.plotting as plotting

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]


def argparser():
    parser = argparse.ArgumentParser(description="Visualize annotated CNVs.")
    parser.add_argument('-f', '--files', nargs="+",
                        help='Path to the result log files.')
    parser.add_argument('-o', '--output', help='Output location.')
    return parser.parse_args()


def main():
    # read cli
    args = argparser()
    args.files.sort(key=alphanum_key)
    # iterate through files list and load the average 10-fold CV
    avg_10_fold_scores = []
    non_pathogenic_support = []
    for file in args.files:
        cv_avg, support  = utils.read_result_file(file)
        avg_10_fold_scores.append(cv_avg)
        non_pathogenic_support.append(support)

    # plot the results with the corresponding allele count values
    plotting.plot_avg_prec_scores(avg_10_fold_scores,[1,2,3,4,5,6,7,8,9,10],support=non_pathogenic_support,save=True,output=args.output)


if __name__ == "__main__":
    main()
