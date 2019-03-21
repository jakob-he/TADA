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


def run(args):
    # load annotated TAD data
    tads = pathlib.Path(args.tads)
    tads = pickle.load(tads.open('rb'))


    # load CNVs
    cnvs = utils.objects_from_file(args.cnvs,'CNV')

    #create cnv dict
    cnvs = utils.create_chr_dictionary_from_beds(cnvs)

    #annotate CNVS
    annotated_cnvs = utils.annotate_cnvs(tads,cnvs)

    #save object as pickle file
    with open("annotated_Eichler_CNVs.p", "wb") as output:
        pickle.dump(annotated_cnvs, output)

def main():
    # parse input arguments
    args = argparser()

    run(args)


if __name__ == "__main__":
    main()
