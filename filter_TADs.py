"""Filter annotated TADs and save them as bed files. The filters are currently hardcoded."""

import argparse
import pathlib
import pickle
import numpy as np
import lib.utils as utils


def argparser():
    parser = argparse.ArgumentParser(description="Filter TADs and save them to bed files.")
    parser.add_argument('-t', '--tads', default='annotated_TADs.p',help='Path to the pickeled TAD file')
    return parser.parse_args()




def main():
    #read cli
    args = argparser()

    # load annotated TAD data
    tad_path = pathlib.Path(args.tads)
    tad_dict = pickle.load(tad_path.open('rb'))

    tads_with_pli_genes = {'1':[],'0.9':[],'0.5':[],'0.1':[],'0':[]}
    tads_with_conserved_enhancer = {'1':[],'0.9':[],'0.5':[],'0.1':[],'0':[]}
    for chrom in tad_dict:
        for tad in tad_dict[chrom]:
            if tad.genes:
                pLIs = [float(gene.data['pLI']) if gene.data['pLI']!='NA' else 0 for gene in tad.genes]
                for threshold in tads_with_pli_genes:
                    if any(pLI >= float(threshold) for pLI in pLIs):
                        tads_with_pli_genes[threshold].append(tad)
            if tad.enhancers:
                conservation_scores =  [float(enhancer.data['Phastcon']) if enhancer.data['Phastcon']!='None' else 0 for enhancer in tad.enhancers]
                for threshold in tads_with_conserved_enhancer:
                    if any(phastcon >= float(threshold) for phastcon in conservation_scores):
                        tads_with_conserved_enhancer[threshold].append(tad)




    for threshold, tads in tads_with_pli_genes.items():
        utils.to_bed(tads,pathlib.Path(f'tads_with_pli_{threshold}_genes.bed'),label=f'TADs_with_pLI_{threshold}_genes')

    for threshold, tads in tads_with_conserved_enhancer.items():
        utils.to_bed(tads,pathlib.Path(f'tads_with_phastcon_{threshold}_enhancer.bed'),label=f'TADs_with_phastcon_{threshold}_enhancer')









if __name__ == '__main__':
    main()
