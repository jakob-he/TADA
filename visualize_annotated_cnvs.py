"""Visualize the annotated SVs"""

# standard libraries
import argparse
import pickle
import pathlib
import numpy as np

# plotting
import matplotlib.pyplot as plt
import seaborn as sns


def argparser():
    parser = argparse.ArgumentParser(description="Visualize annotated CNVs.")
    parser.add_argument('-c', '--cnvs', nargs="+", default='annotated_CNVs.p',
                        help='Path to the pickeled CNV file')
    return parser.parse_args()


def main():
    # parse input arguments
    args = argparser()

    colors = ['red','green']

    sns.set_style("ticks")
    sns.set_context("paper")
    sns.set_palette("Paired")
    plt.figure()
    for idx, cnv_set in enumerate(args.cnvs):
            cnv_counter = 0
            # load annotated CNV data
            cnv_path = pathlib.Path(cnv_set)
            cnv_set = pickle.load(cnv_path.open('rb'))

            #Iterate through every cnv and get the distance to the closest gene/enhancer
            gene_distances = []
            enhancer_distances = []
            boundary_breaking = []
            for chrom in cnv_set:
                for cnv in cnv_set[chrom]:
                    if cnv.tads and cnv.gene_distances and cnv.enhancer_distances:
                        cnv_counter+=1
                        gene_distances.append(cnv.gene_distances[0])
                        enhancer_distances.append(cnv.enhancer_distances[0])
                        boundary_breaking.append(cnv.boundary_spanning)

            #plot the gene_distances
            gene_distances = [np.log10(value) if value!=0 else value for value in gene_distances]
            sns.kdeplot(gene_distances,label=f'Gene distance {cnv_path.stem} (N={cnv_counter})').set(xlim=(0))

            #plot the enhancer_distances
            gene_distances = [np.log10(value) if value!=0 else value for value in enhancer_distances]
            sns.kdeplot(gene_distances,label=f'Enhancer distance {cnv_path.stem}').set(xlim=(0))




    plt.legend()
    plt.xlabel('Distance log10(bp)')
    sns.despine(top=True,right=True)
    plt.show()




if __name__ == "__main__":
    main()
