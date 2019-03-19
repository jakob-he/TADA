# standart libraries
import argparse
import json

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
        '-g', '--genes', default='data/GENE/genes.bed', help='Path to the Gene BED file.')
    parser.add_argument('-e', '--enhancer', default='data/ENHANCER/FANTOM.bed',
                        help='Path to the Enhancer BED file.')
    return parser.parse_args()

def run(args):
    # create bed objects from TADS
    tad_beds = utils.objects_from_file(args.tads,Tad)
    enhancer_beds = utils.objects_from_file(args.enhancer,Enhancer,column_names=['RANDOM'])
    gene_beds = utils.objects_from_file(args.genes,Gene,column_names=['ENSEMBLE_ID'])

    #create dict with chromsomes as keys
    gene_dict = utils.create_chr_dictionary_from_beds(gene_beds)
    enhancer_dict = utils.create_chr_dictionary_from_beds(enhancer_beds)
    tad_dict = utils.create_chr_dictionary_from_beds(tad_beds)

    #Annotate TADs with overlapping enhancer and genes
    annotated_tads = utils.create_annotated_dict(tad_dict,gene_dict,enhancer_dict)



    return annotated_tad



def main():
    # parse input arguments
    args = argparser()

    run(args)

if __name__ == '__main__':
    main()
