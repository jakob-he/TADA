"""Annotate CNVs with either a user-defined or predefined feature set."""
# standard libraries
import argparse
import pickle
import pathlib
import yaml

# own libraries
import tada.annotate_tads as annotate_tads
import tada.lib.utils as utils
import tada.lib.preprocessing as preprocessing
import tada.lib.plotting as plotting

# third party libraries
import pandas as pd
import matplotlib.pyplot as plt

def argparser():
    parser = argparse.ArgumentParser(description="Annotate CNVs. Run annotate_cnvs -h for more details.")
    parser.add_argument('-c', '--config', default = 'config.yml', help='Path to the config file containing TAD,CNV and annotation locations.')
    parser.add_argument('-p', '--pickle', action='store_true', help='Save annotated CNV objects as pickled file. Default is false.')
    parser.add_argument('-o', '--output', default='./', help='Directory for the output files.')
    return parser.parse_args()


def annotate(cfg):
    # check if the TAD file is a pickle file
    # if not the TADs need to be annoted first

    if cfg['TADS']['ANNOTATED']:
        with pathlib.Path(cfg['TADS']['ANNOTATED']).open('rb') as annotated_tads:
            tads = pickle.load(annotated_tads)
    else:
        tads = annotate_tads.annotate(cfg)
    # for each set of cnvs in the config file
    # read them into a dict with chromosomes as keys
    # and annotate each CNVs according to the TAD environment
    annotated_cnvs = {}
    for cnvs in cfg['CNVS']['RAW']:
        cnv_bed = utils.objects_from_file(cfg['CNVS']['RAW'][cnvs],'CNV')
        cnv_dict = utils.create_chr_dictionary_from_beds(cnv_bed)
        # annotate CNVS
        print(f'Annotate {cnvs} CNVs..')
        annotated_cnvs[cnvs] = utils.annotate_cnvs(tads, cnv_dict)
    return annotated_cnvs


def main():
    # parse input arguments
    args = argparser()

    # read config file
    with pathlib.Path(args.config).open() as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.Loader)


    annotated_cnvs = annotate(cfg)

    # get labels depending on the feature type
    if cfg['FEATURES'] == 'extended':
        feature_labels = ['Number of affected Genes','Number of affected Enhancers','Boundary Distance', 'Boundary Stability', 'Gene Distance', 'Enhancer Distance', 'DDG2P Distance', 'Gene LOEUF','Enhancer conservation', 'Gene HI', 'CTCF Distance', 'HI LogOdds Score', 'Exon Overlap', 'MPOI']
    else:
        feature_labels = [f'{annotation} distance' for annotation in cfg['ANNOTATIONS']]

    # save CNVs as csv-file and, if necessary, as pickle file,
    output_path = pathlib.Path(args.output)
    for label, cnvs in annotated_cnvs.items():
        if args.pickle:
            with open(output_path / f'Annotated_{label}.p', "wb") as output:
                pickle.dump(cnvs, output)
        feature_df = preprocessing.create_feature_df(cnvs,cfg['FEATURES'],feature_labels,csv=True)
        feature_df.to_csv(output_path / f'Annotated_{label}.csv',sep='\t',header=True,index=False)
        print(f'{label} CNVs saved in {output_path / f"Annotated_{label}.csv"}')


if __name__ == "__main__":
    main()
