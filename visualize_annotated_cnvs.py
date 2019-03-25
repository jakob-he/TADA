"""Visualize the annotated SVs"""

# standard libraries
import argparse
import pickle
import pathlib

# plotting
import matplotlib.pyplot as plt

# own libraries
import lib.plotting as plotting


def argparser():
    parser = argparse.ArgumentParser(description="Visualize annotated CNVs.")
    parser.add_argument('-c', '--cnvs', default='annotated_CNVs.p',
                        help='Path to the pickeled TAD file')
    return parser.parse_args()


def main():
    # parse input arguments
    args = argparser()

    # load annotated CNV data
    cnvs = pathlib.Path(args.cnvs)
    cnvs = pickle.load(cnvs.open('rb'))

    #see how many cnvs are boundary_spanning
    for chrom in cnvs:
        for cnv in cnvs[chrom]:
            if len(cnv.tads)>1:
                print(cnv)


if __name__ == "__main__":
    main()
