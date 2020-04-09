"""Annotate TADs with a set of functional annotations."""
# standard libraries
import argparse
import pickle
import pathlib
import yaml

# own libraries
from lib.tad import Tad
import lib.utils as utils


def argparser():
    parser = argparse.ArgumentParser(
        description="Annotate TADs. Run annotate_tads -h for more details")
    parser.add_argument('-c', '--config', default='config.yml', help='Path')
    parser.add_argument('-o', '--output', default='./',
                        help='Path to the output directory.')
    return parser.parse_args()


def annotate(cfg):
    print('Annotating TADs...')

    # create bed objects from TADS
    tad_beds = utils.objects_from_file(cfg['TADS']['RAW'], 'TAD')

    # create dict with chromsomes as keys
    tad_dict = utils.create_chr_dictionary_from_beds(tad_beds)

    # create annotation dicts

    # if the extended feature set is to be used:
    # annotate genes with exons and P-O-interactions
    # and add all other annotations to the annotation dict
    annotation_dicts = {}
    if cfg['FEATURES'] == 'extended':
        for annotation_name in cfg['ANNOTATIONS']:
            if annotation_name != 'GENES':
                annotation_dicts[annotation_name] = utils.create_chr_dictionary_from_beds(utils.objects_from_file(
                    cfg['ANNOTATIONS'][annotation_name], 'BED'))
            else:
                gene_dict = utils.create_chr_dictionary_from_beds(
                    utils.objects_from_file(cfg['ANNOTATIONS']['GENES'], 'GENE'))

        annotated_genes = utils.create_annotated_bed_dict(gene_dict, {
                                                          **{'EXONS': annotation_dicts['EXONS']}, **{'POINT': annotation_dicts['POINT']}}, gene_annotation=True)
        annotation_dicts['GENES'] = annotated_genes

        # remove exons and P-O-interactions from the dict since they are no longer needed
        for key in ['EXONS', 'POINT']:
            del annotation_dicts[key]

    # otherwise just add the annotations to the annotation dict
    else:
        for annotation_name in cfg['ANNOTATIONS']:
            annotation_dicts[annotation_name] = utils.create_chr_dictionary_from_beds(utils.objects_from_file(
                cfg['ANNOTATIONS'][annotation_name], 'BED'))

    # Annotate TADs with overlapping annotations
    annotated_tads = utils.create_annotated_bed_dict(
        tad_dict, annotation_dicts)

    return annotated_tads


def main():
    # parse input arguments
    args = argparser()

    # get config data
    with pathlib.Path(args.config).open() as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.Loader)

    annotated_tads = annotate(cfg)

    # save object as pickle file
    output_dir = pathlib.Path(args.output)
    with open(output_dir / 'Annotated_TADs.p', 'wb') as output:
        pickle.dump(annotated_tads, output)

    print(f'Annotated TADs saved in {output_dir.absolute() / "Annotated_TADs.p"}')


if __name__ == '__main__':
    main()
