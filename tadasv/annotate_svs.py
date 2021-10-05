"""Annotate SVs with either a user-defined or predefined feature set."""
# standard libraries
import argparse
import pickle
import pathlib
import yaml
import pkg_resources

# own libraries
import tadasv.annotate_tads as annotate_tads
import tadasv.lib.utils as utils
import tadasv.lib.preprocessing as preprocessing
import tadasv.lib.plotting as plotting

# third party libraries
import pandas as pd
import matplotlib.pyplot as plt

def argparser():
    parser = argparse.ArgumentParser(description="Annotate SVs. Run annotate_svs -h for more details.")
    parser.add_argument('-c', '--config', default = 'config.yml', help='Path to the config file containing TAD,SV and annotation locations.')
    parser.add_argument('-d', '--default', action='store_true', help='Use default settings and annptations. This requires the variant path to be set via -v!')
    parser.add_argument('-v', '--variants', help='Path to the SV bed/vcf-file. Only usable in combination with the -d flag!')
    parser.add_argument('-t', '--type', help='Type of variant (Either DEL or DUP). Only usable in combination with the -d flag!')
    parser.add_argument('-p', '--pickle', action='store_true', help='Save annotated SV objects as pickled file. Default is false.')
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
    # for each set of svs in the config file
    # read them into a dict with chromosomes as keys
    # and annotate each SVs according to the TAD environment
    annotated_svs = {}
    for svs in cfg['SVS']['RAW']:
        sv_bed = utils.objects_from_file(cfg['SVS']['RAW'][svs],'SV')
        sv_dict = utils.create_chr_dictionary_from_beds(sv_bed)
        # annotate SVS
        print(f'Annotate {svs} SVs..')
        annotated_svs[svs] = utils.annotate_svs(tads, sv_dict)
    return annotated_svs


def main():
    # parse input arguments
    args = argparser()

    if args.default:
        # depending on the variant type get config file from package
        if args.type.upper() == 'DEL':
            cfg_stream = pathlib.Path(pkg_resources.resource_filename(__name__, 'data/config_del_default.yml'))
        elif args.type.upper() == 'DUP':
            cfg_stream = pathlib.Path(pkg_resources.resource_filename(
                __name__, 'data/config_dup_default.yml'))
        else:
            print(f'{args.type} is not supported. Only DEL and DUP are viable options!')
            return
        cfg = yaml.load(cfg_stream.open(), Loader=yaml.Loader)
        cfg['SVS']['RAW']['TEST'] = args.variants
        # iterate through cfg entires and set path to package path
        for key_1 in cfg:
            if key_1 == "PRETRAINED_MODEL":
                cfg[key_1] = cfg_stream.parent / cfg[key_1]
            elif key_1 not in ['SVS', 'FEATURES', 'KWARGS', 'CLASSIFIER']:
                for key_2 in cfg[key_1]:
                    cfg[key_1][key_2] =  cfg_stream.parent / cfg[key_1][key_2]
    else:
        # read config file
        with pathlib.Path(args.config).open() as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.Loader)
            
    annotated_svs = annotate(cfg)

    # get labels depending on the feature type
    if cfg['FEATURES'] == 'extended':
        feature_labels = ['PLS Overlap (bp)', 'dELS Overlap (bp)', 'pELS Overlap (bp)', 'Low-DNase Overlap (bp)', 'Stop Codon', 'Start Codon', '5_UTR', '3_UTR', 'Gene Overlap (bp)', 'Enhancer Overlap (bp)', 'DDG2P Overlap (bp)', 'Number of affected Genes', 'Number of affected Enhancers', 'Boundary Distance',
                          'Boundary Stability', 'Gene Distance', 'Enhancer Distance', 'DDG2P Distance', 'Gene LOEUF', 'Enhancer conservation', 'Gene HI', 'CTCF Overlap (bp)', 'HI LogOdds Score', 'Exon Overlap', 'MPOI']
    else:
        feature_labels = [f'{annotation} distance' for annotation in cfg['ANNOTATIONS']]

    # save SVs as csv-file and, if necessary, as pickle file,
    output_path = pathlib.Path(args.output)
    for label, svs in annotated_svs.items():
        if args.pickle:
            with open(output_path / f'Annotated_{label}.p', "wb") as output:
                pickle.dump(svs, output, protocol=4)
        feature_df = preprocessing.create_feature_df(svs,cfg['FEATURES'],feature_labels,csv=True)
        feature_df.to_csv(output_path / f'Annotated_{label}.csv',sep='\t',header=True,index=False)
        print(f'{label} SVs saved in {output_path / f"Annotated_{label}.csv"}')


if __name__ == "__main__":
    main()