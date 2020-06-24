"""Plots the distribution of variants across the genome, stained by data source."""

import argparse
import pathlib
import pandas as pd

# plotting
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib.lines import Line2D


def chromosome_collections(df, y_positions, height,  **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


def plot_chromsome_distribution(ideo,variants,ax):
    # Height of each ideogram
    chrom_height = 0.5

    # Spacing between consecutive ideograms
    chrom_spacing = 1

    # Height of the variant track. Should be smaller than `chrom_spacing` in order to
    # fit correctly
    variant_height = 0.8

    # Padding between the top of a gene track and its corresponding ideogram
    variant_padding = 0.1

    # Decide which chromosomes to use
    chromosome_list = ['chr%s' % i for i in range(1, 23)]

    # Keep track of the y positions for ideograms and genes for each chromosome,
    # and the center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    variant_ybase = {}
    chrom_centers = {}

    # Iterate in reverse so that items in the beginning of `chromosome_list` will
    # appear at the top of the plot
    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        variant_ybase[chrom] = ybase - variant_height - variant_padding
        ybase += chrom_height + chrom_spacing

    # Filter out chromosomes not in our list
    ideo = ideo[ideo['chrom'].apply(lambda x: x in chromosome_list)]

    # Add a new column for width
    ideo['width'] = ideo.end - ideo.start

    # Colors for different chromosome stains and variant sources
    color_lookup_ideogram = {
        'gneg': (1., 1., 1.),
        'gpos25': (.6, .6, .6),
        'gpos50': (.4, .4, .4),
        'gpos75': (.2, .2, .2),
        'gpos100': (0., 0., 0.),
        'acen': (.8, .4, .4),
        'gvar': (.8, .8, .8),
        'stalk': (.9, .9, .9),
    }

    color_lookup_variants = {
        'GnomAD': '#e01b22',
        'Eichler': '#22e01b',
        'Biobank': '#1b28e0',
        'DGV': '#e07a1b'
    }

    # Add a new column for colors
    ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup_ideogram[x])

    # Same thing for the variants
    variants = variants[variants['chrom'].apply(lambda x: x in chromosome_list)]
    variants['width'] = variants.end - variants.start
    variants['colors'] = variants['origin'].apply(
        lambda x: color_lookup_variants[x])

    # Now all we have to do is call our function for the ideogram data...
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, linewidths=1, edgecolors='black'):
        ax.add_collection(collection)

    # ...and the gene data
    for collection in chromosome_collections(
        variants, variant_ybase, variant_height, alpha=0.5, linewidths=0
    ):
        ax.add_collection(collection)

    # add custom legend
    custom_lines = [Line2D([0], [0], color=color_lookup_variants['GnomAD'], lw = 3),
                    Line2D([0], [0], color=color_lookup_variants['Eichler'], lw = 3),
                    Line2D([0], [0], color=color_lookup_variants['Biobank'], lw = 3),
                    Line2D([0], [0], color=color_lookup_variants['DGV'], lw = 3)]

    ax.legend(custom_lines, ['GnomAD', 'Eichler', 'Biobank', 'DGV'],loc='lower right')

    # Axes tweaking
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    ax.axis('tight')
    return ax

if __name__ == "__main__":
    main()
