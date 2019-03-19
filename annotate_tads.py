# standart libraries
import argparse

# own libraries
from lib.bed import Bed
from lib.bed_parser import bed_objects_from_file


def argparser():
    parser = argparse.ArgumentParser(description="Annotate TADs.")
    parser.add_argument('-t', '--tads', default='data/TAD/hg38/H1-ESC_Dixon_2015-raw_TADs.txt',
                        help='Path to the TAD boundary BED file.')
    parser.add_argument(
        '-g''--genes', default='data/GENE/genes.bed', help='Path to the Gene BED file.')
    parser.add_argument('-e''--enhancer', default='data/ENHANCER/FANTOM.bed',
                        help='Path to the Enhancer BED file.')
    return parser.parse_args()


def main():
    # parse input arguments
    args = argparser()

    # create bed objects from TADS
    tad_beds = bed_objects_from_file(args.tads)
    enhancer_beds = bed_objects_from_file(args.enhancer)


if __name__ == '__main__':
    main()
