"""Returns a single bed element for each cluster of elements with the specified reciprocal overlap"""

import argparse
import pathlib
import numpy as np
from collections import defaultdict


def argparser():
    parser = argparse.ArgumentParser('Selects the representitive with the smallest size from a cluster with the defined reciprocal overlap.')
    parser.add_argument('-bed','--bedfile',help='Path to the bedfile.')
    parser.add_argument('-r','--reciprocal',default=0.9,help='Reciprocal overlap defined as a fraction (0-1).')
    parser.add_argument('-o','--output',help='Output file location.')
    parser.add_argument('-p','--pathogenic',action='store_true',help='Select variant without regarding the fourth column.')
    return parser.parse_args()

def create_chrom_dict(file_path, pathogenic):
    chrom_dict = defaultdict(list)
    if pathogenic:
        for line in file_path.open():
            chr, start, end, patho = line.rstrip('\n').split('\t')[0:4]
            chrom_dict[chr].append([int(start),int(end),patho])
    else:
        for line in file_path.open():
            chr, start, end, data_set = line.rstrip('\n').split('\t')[0:4]
            chrom_dict[chr].append([int(start),int(end),data_set])
    return chrom_dict

def get_smallest_variant(variant_set, pathogenic):
    sizes = [variant[1]-variant[0] for variant in variant_set]
    if pathogenic:
        return variant_set[np.argmin(sizes)]
    else:
        data_sets = [variant[2] for variant in variant_set]
        if any(data_set == 'Eichler' for data_set in data_sets):
            eichler_variants = [variant for i, variant in enumerate(variant_set) if data_sets[i] == 'Eichler']
            eichler_sizes = [variant[1]-variant[0] for variant in eichler_variants]
            return eichler_variants[np.argmin(eichler_sizes)]
        elif any(data_set == 'GnomAD' for data_set in data_sets):
            gnomad_variants = [variant for i, variant in enumerate(variant_set) if data_sets[i] == 'GnomAD']
            gnomad_sizes = [variant[1]-variant[0] for variant in gnomad_variants]
            return gnomad_variants[np.argmin(gnomad_sizes)]
        else:
            return variant_set[np.argmin(sizes)]

def test_overlap(interval_a,interval_b,overlap):
    overlap = max(0,min(interval_a[1],interval_b[1])-max(interval_a[0],interval_b[0]))
    if overlap != 0:
        rec_overlap_1 = overlap / (interval_a[1]-interval_a[0])
        rec_overlap_2 = overlap / (interval_b[1]-interval_b[0])
        if rec_overlap_1 >= 0.9 and rec_overlap_2 >= 0.9:
            return True

def run(args):
    bed_path = pathlib.Path(args.bedfile)
    output_path = pathlib.Path(args.output)
    rec_overlap = float(args.reciprocal)

    # create chromosome dictionary from bed file
    chrom_dict = create_chrom_dict(bed_path,args.pathogenic)

    # get overlapping intervals and cluster them
    # for each cluster select the smallest variant
    # and save it in a new chromosome dict
    new_chrom_dict = defaultdict(list)
    for chrom in chrom_dict.keys():
        cluster_dict = defaultdict(list)
        for i, variant_1 in enumerate(chrom_dict[chrom]):
            for j, variant_2 in enumerate(chrom_dict[chrom]):
                    if test_overlap(variant_1,variant_2,rec_overlap):
                        cluster_dict[i].append(j)


        # iterate through each cluster i.e. variant index
        # get all overlapping variants if not present in an earlier cluster
        # select the smallest variant of the cluster for the new variant set
        for cluster_id in list(cluster_dict.keys()):
            overlapping_idx = []
            if cluster_dict[cluster_id]:
                for id in cluster_dict[cluster_id]:
                    overlapping_idx.extend(cluster_dict[id])
                unique_idx = list(set([id for id in overlapping_idx if cluster_dict[id]]))
                for idx in unique_idx:
                    del cluster_dict[idx]
                if len(unique_idx) > 1:
                    smallest_variant = get_smallest_variant([chrom_dict[chrom][idx] for idx in unique_idx], args.pathogenic)
                else:
                    smallest_variant = chrom_dict[chrom][unique_idx[0]]
                new_chrom_dict[chrom].append(smallest_variant)

    with output_path.open('w+') as output:
        for chrom in new_chrom_dict:
            for variant in new_chrom_dict[chrom]:
                output.write(f'{chrom}\t{variant[0]}\t{variant[1]}\t{variant[2]}\n')


def main():
    args = argparser()

    run(args)

if __name__ == '__main__':
    main()
