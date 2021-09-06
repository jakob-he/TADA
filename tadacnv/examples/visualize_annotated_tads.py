"""Visualize the annotated TADs"""

# standard libraries
import argparse
import pickle
import pathlib

# plotting
import matplotlib.pyplot as plt

# own libraries
import lib.plotting as plotting


def argparser():
    parser = argparse.ArgumentParser(description="Visualize annotated TADs.")
    parser.add_argument('-t', '--tads', default='annotated_TADs.p',
                        help='Path to the pickeled TAD file')
    return parser.parse_args()


def main():
    # parse input arguments
    args = argparser()

    # load annotated TAD data
    tads = pathlib.Path(args.tads)
    tads = pickle.load(tads.open('rb'))

    #plot the distribution of pLI values

    # plot number of genes and enhancers for chromosome 1
    #plt.figure()
    #ax = plotting.plot_tad_element_dist(tads['chr1'])
    #ax.set_title('Number of overlapping Genes / Enhancer')
    #ax.set_xticklabels(['Genes', 'Enhancer'])
    #plt.show()

    #TODO plot number of genes for all chromsomes
    # plt.figure()
    # for chrom in tads:
    #     plotting.plot_tad_element_dist(tads[chrom])
    # ax.set_title('Number of overlapping Genes per Chromsomes')
    # ax.set_xticklabels(tads.keys())
    # plt.show()

if __name__ == "__main__":
    main()
