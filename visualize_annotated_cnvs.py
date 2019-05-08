"""Visualize the annotated SVs"""

# standard libraries
import argparse
import pickle
import pathlib
import matplotlib.pyplot as plt

# plotting
import matplotlib.pyplot as plt

# own libraries
import lib.plotting as plotting


def argparser():
    parser = argparse.ArgumentParser(description="Visualize annotated CNVs.")
    parser.add_argument('-c', '--cnvs', nargs="+", default='annotated_CNVs.p',
                        help='Path to the pickeled TAD file')
    return parser.parse_args()


def main():
    # parse input arguments
    args = argparser()


    for cnv_set in args.cnvs:
            # load annotated CNV data
            cnv_set = pathlib.Path(cnv_set)
            cnv_set = pickle.load(cnv_set.open('rb'))

            #Iterate through every cnv and get the distance to the closest gene/enhancer
            gene_distances = []
            enhancer_distances = []
            for chrom in cnv_set:
                for cnv in cnv_set[chrom]:
                    if cnv.tads and cnv.gene_distances:
                        gene_distances.append(cnv.gene_distances[0])

            #plot the gene_distances
            plt.figure()
            plt.hist(gene_distances,bins=200, histtype = 'step')
            plt.show()



if __name__ == "__main__":
    main()
