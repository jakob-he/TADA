"""Load and annotate a set of CNVs"""
# standard libraries
import argparse
import pickle
import pathlib

#own libraries
import lib.utils as utils


def argparser():
    parser = argparse.ArgumentParser(description="Annotate a set of CNVs.")
    parser.add_argument('-t', '--tads', default='annotated_TADs.p',
                        help='Path to the pickeled TAD file.')
    parser.add_argument('-c', '--cnvs', help='Path to the CNV file.')
    return parser.parse_args()


def main():
    # parse input arguments
    args = argparser()

    # load annotated TAD data
    tads = pathlib.Path(args.tads)
    tads = pickle.load(tads.open('rb'))

    # load CNVs
    cnvs = utils.objects_from_file(args.cnvs,'CNV')

    #create cnv dict
    cnvs = utils.create_chr_dictionary_from_beds(cnvs)

    #annotate CNVS
    annotated_cnvs = utils.annotate_cnvs(tads,cnvs)

    for cnv in annotated_cnvs['chr1']:
        print(cnv)


if __name__ == "__main__":
    main()
