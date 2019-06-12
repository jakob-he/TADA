"""Load and annotate a set of CNVs"""
#standard libraries
import argparse
import pickle
import pathlib

#own libraries
import lib.utils as utils

#third party libraries
import pandas as pd

def argparser():
    parser = argparse.ArgumentParser(description="Annotate a set of CNVs.")
    parser.add_argument('-t', '--tads', default='annotated_TADs.p',
                        help='Path to the pickeled TAD file.')
    parser.add_argument('-c', '--cnvs', help='Path to the CNV file.')
    parser.add_argument('-vcf', '--vcf', action='store_true', help='Needs to be set if the CNV file is a VCF with the location at the second position.')
    parser.add_argument('-csv', '--csv', action='store_true', help='Save a CSV file in additon to the pickled object.')
    parser.add_argument('-o', '--output', default='annotated_CNVS.p', help='Output File.')
    return parser.parse_args()


def run(args):
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
    with open(output_path, "wb") as output:
        pickle.dump(annotated_cnvs, output)

    if args.csv:
        output_df = {'chr':[],'start':[],'stop':[],'gene':[],'hgnc_id':[],'type':[],'mode':[],'mech':[],'syndrome':[],'TAD_boundary_overlap':[],'ctcf_overlap':[]}
        for chrom in annotated_cnvs:
            for cnv in annotated_cnvs[chrom]:
                genes = []
                hgnc_id = []
                type = []
                mode = []
                mech = []
                syndrome = []
                ctcf_overlap = []
                boundary = []
                if cnv.tads:
                    for tad in cnv.tads:
                        if tad.annotations['genes']:
                            for gene in tad.annotations['genes']:
                                genes.append(gene.data['gene'])
                                hgnc_id.append(gene.data['hgnc_id'])
                                type.append(gene.data['type'])
                                mode.append(gene.data['mode'])
                                mech.append(gene.data['mech'])
                                syndrome.append(gene.data['syndrome'])

                    if any(ctcf_distance for ctcf_distance in cnv.annotation_distances['ctcf']):
                        ctcf_overlap.append('YES')
                    else:
                        ctcf_overlap.append('NO')
                else:
                    ctcf_overlap.append('NO')
                if cnv.boundary_spanning:
                    boundary.append('YES')
                else:
                    boundary.append('NO')

                output_df['chr'].append(cnv.chr)
                output_df['start'].append(cnv.start)
                output_df['stop'].append(cnv.end)
                output_df['gene'].append(', '.join(genes))
                output_df['hgnc_id'].append(', '.join(hgnc_id))
                output_df['type'].append(', '.join(type))
                output_df['mode'].append(', '.join(mode))
                output_df['mech'].append(', '.join(mech))
                output_df['syndrome'].append(', '.join(syndrome))
                output_df['ctcf_overlap'].append(', '.join(ctcf_overlap))
                output_df['TAD_boundary_overlap'].append(', '.join(boundary))


        # save annotated cnvs as CSV
        output_df = pd.DataFrame(output_df)
        output_df.drop_duplicates(inplace=True)
        output_df.to_csv(output_path.stem + '.csv',sep='\t',index=False)


def main():
    #parse input arguments
    args = argparser()

    run(args)


if __name__ == "__main__":
    main()
