"""Annotate tads with the overlapping genes and enhancers"""
# standard libraries
import argparse
import pickle
import pathlib

# own libraries
from lib.tad import Tad
from lib.gene import Gene
from lib.enhancer import Enhancer
import lib.utils as utils


def argparser():
    parser = argparse.ArgumentParser(description="Annotate TADs.")
    parser.add_argument('-t', '--tads', default='data/TAD/hg38/H1-ESC_Dixon_2015-raw_TADs.txt',
                        help='Path to the TAD boundary BED file.')
    parser.add_argument(
        '-an','--annotation_names',nargs='*',help='Names of the annotations (e.g. genes).')
    parser.add_argument('-af','--annotation_files',nargs='*',help='Paths to the annotation files.')
    parser.add_argument('-o','--output', default='annotated_TADs.p',help='Output file.')
    return parser.parse_args()


def run(args):
    # create bed objects from TADS
    tad_beds = utils.objects_from_file(args.tads, 'TAD')

    # create empty enhancer and gene dicts
    annotation_dicts = {}
    for idx, annotation_name in enumerate(args.annotation_names):
        annotation_dicts[annotation_name] = utils.create_chr_dictionary_from_beds(utils.objects_from_file(
            args.annotation_files[idx], 'Bed'))

    # create dict with chromsomes as keys
    tad_dict = utils.create_chr_dictionary_from_beds(tad_beds)

    # Annotate TADs with overlapping enhancer and genes
    annotated_tads = utils.create_annotated_tad_dict(
        tad_dict, annotation_dicts)

    return annotated_tads


def main():
    # parse input arguments
    args = argparser()

    annotated_tads = run(args)

    #save object as pickle file
    with open(args.output, "wb") as output:
        pickle.dump(annotated_tads, output)


if __name__ == '__main__':
    main()
