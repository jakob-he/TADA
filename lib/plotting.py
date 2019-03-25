"""Functions that can be used to visualize bed files and annotate TADs"""
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns


def plot_size_dist(beds, output='./', plot_type='h', save=False):
    """Plots the size distribution of elements in a bed file.
    Args:
        beds: List of bed objects (e.g. genes).
        output: Path were the figure is saved.
        plot_type: 'h' for histogram, 'b' for boxplot, 'v' for violinplot
        save: True if figure should be saved rather than shown directly.
    """
    sizes = [bed.end - bed.start for bed in beds]
    plt.figure()
    sns.set_style('whitegrid')
    if plot_type == 'h':
        ax = sns.distplot(sizes, kde=False)
    elif plot_type == 'b':
        ax = sns.boxplot(y=sizes)
    elif plot_type == 'v':
        ax = sns.violinplot(y=sizes)

    ax.set_title('Size Distribution')
    ax.set_ylabel('Amount')
    ax.set_xlabel('Size')
    if save:
        plt.savefig(output / 'size_histogram.png')
    else:
        plt.show()

def plot_annotation_dist(beds, annotation, output='./',plot_type='h', save = False):
    """Plots the distribution of an annotation in the bed file in a box plot.
    Args:
        beds: List of bed objects (e.g. genes).
        annotation: Name of the column which is going to be plotted. (e.g. conservation)
        output: Path were the figure is saved.
        plot_type: 'h' for histogram, 'b' for boxplot, 'v' for violinplot
        save: True if figure should be saved rather than shown directly.
    """
    plt.figure()
    sns.set_style('whitegrid')
    annotations = [bed.data[annotation] for bed in beds]
    if plot_type == 'h':
        ax = sns.distplot(annotations, kde=False)
    elif plot_type == 'b':
        ax = sns.boxplot(y=annotations)
    elif plot_type == 'v':
        ax = sns.violinplot(y=annotations)
    ax.set_title(f'Distribution of {annotation}')
    ax.set_ylabel('Amount')
    ax.set_xlabel(f'{annotation}')
    if save:
        plt.savefig(output / f'{annotation}_distribution.png')
    else:
        plt.show()



def plot_tad_element_dist(tads, output='./', save=False, genes=True, enhancer=True):
    """Plots the distribution of the number of genes and enhancers in the list of Tads"""
    length = 1
    if genes and enhancer:
        genes = []
        enhancer = []
        for tad in tads:
            genes.append(tad.count_genes())
            enhancer.append(tad.count_enhancer())
        elements = [genes, enhancer]
        length = 2
    elif genes and not enhancer:
        elements = [tad.count_genes() for tad in tads]
    elif enhancer and not genes:
        elements = [tad.count_enhancer() for tad in tads]

    sns.set_style('whitegrid')
    ax = sns.boxplot(x=np.arange(0,length), y=elements)
    return ax
