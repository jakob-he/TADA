"""Load and annotate a set of CNVs"""
#standard libraries
import argparse
import pickle
import pathlib

#own libraries
import lib.utils as utils
import lib.preprocessing as preprocessing
import lib.plotting as plotting

#third party libraries
import pandas as pd
import matplotlib.pyplot as plt

def argparser():
    parser = argparse.ArgumentParser(description="Annotate a set of CNVs. Run annotate_cnvs -h for more details.")
    parser.add_argument('-t', '--tads',default='data/default_annotated_TADs.p',
                        help='Path to the pickeled TAD file.')
    parser.add_argument('-c', '--cnvs', help='Path to the CNV file.', required=True)
    parser.add_argument('-vcf', '--vcf', action='store_true', help='Needs to be set if the CNV file is a VCF with the location at the second position.')
    parser.add_argument('-p', '--pickle', action='store_false', help='Save annotated CNV objects as pickled file. Default is True.')
    parser.add_argument('-csv', '--csv', action='store_false', help='Save a CSV file in additon to the pickled object. Specific features sets can be definied with -f.')
    parser.add_argument('-pl', '--plot', action='store_true', help='Plot Distribution of features for CNV set. Default is false.')
    parser.add_argument('-f','--features',default='extended_continuous',help='Features for the CSV file. TADs need to be annotated with the corresponding features.')
    parser.add_argument('-o', '--output', default='./', help='Output File.')
    return parser.parse_args()


def run(args):
    print('Annotating CNVs...')
    #load annotated TAD data
    tads = pathlib.Path(args.tads)
    tads = pickle.load(tads.open('rb'))
    output_path = pathlib.Path(args.output)

    #load CNVs
    cnvs = utils.objects_from_file(args.cnvs,'CNV',vcf=args.vcf)

    #create cnv dict
    cnvs = utils.create_chr_dictionary_from_beds(cnvs)

    #annotate CNVS
    annotated_cnvs = utils.annotate_cnvs(tads,cnvs)
    
    #save raw CNV object as pickle file
    with open(output_path / 'Annotated_CNVs.p', "wb") as output:
        pickle.dump(annotated_cnvs, output)

    if args.csv:
        feature_df = preprocessing.create_feature_df(annotated_cnvs,args.features,csv=True)
        feature_df.to_csv(output_path / 'Annotated_CNVs.csv',sep='\t',header=True,index=False)

        #if args.plot:
        #    plotting.plot_feature_dist(feature_df, exclude_features=['CHR','START','END'])
        #    plt.savefig(output_path / 'Annotated_CNVs.png',bbox_inches='tight')



def main():
    #parse input arguments
    args = argparser()

    run(args)


if __name__ == "__main__":
    main()
