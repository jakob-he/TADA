"""Annotate tads with the overlapping genes and enhancers"""
# standard libraries
import argparse
import pickle
import pathlib

# own libraries
from lib.tad import Tad
import lib.utils as utils


def argparser():
    parser = argparse.ArgumentParser(description="Annotate TADs. Run annotate_tads -h for more details")
    parser.add_argument('-t', '--tads',
                        help='Path to the TAD boundary BED file.', required=True)
    parser.add_argument(
        '-an','--annotation_names',nargs='*',help='Names of the annotations. ')
    parser.add_argument('-af','--annotation_files',nargs='*',help='Paths to the annotation files.')
    parser.add_argument('-ga','--gene_annotation',action='store_true',help='Annotate Genes with Exons. The gene (-g) and exon (-e) bed files have to be included.')
    parser.add_argument('-g','--genes',help='Path to a BED file containing gene annotations.')
    parser.add_argument('-e','--exons',help='Path to a BED file containing exon annotations with Gene IDs in the 4th column.')
    parser.add_argument('-o','--output', default='annotated_TADs.p', help='Output file.')
    return parser.parse_args()


def run(args):
    # create bed objects from TADS
    tad_beds = utils.objects_from_file(args.tads, 'TAD')

    # create annotation dicts
    annotation_dicts = {}
    for idx, annotation_name in enumerate(args.annotation_names):
        annotation_dicts[annotation_name] = utils.create_chr_dictionary_from_beds(utils.objects_from_file(
            args.annotation_files[idx], 'Bed'))

    # if needed annotate genes with enhancers
    if args.gene_annotation:
        exon_dict = {'exons':utils.create_chr_dictionary_from_beds(utils.objects_from_file(args.exons,'Bed'))}
        gene_dict = utils.create_chr_dictionary_from_beds(utils.objects_from_file(args.genes,'Gene'))
        annotated_genes = utils.create_annotated_bed_dict(gene_dict,exon_dict,filer_exons=True)
        annotation_dicts['genes'] = annotated_genes

    # create dict with chromsomes as keys
    tad_dict = utils.create_chr_dictionary_from_beds(tad_beds)

    # Annotate TADs with overlapping annotations
    annotated_tads = utils.create_annotated_bed_dict(
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
