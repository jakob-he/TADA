"""Functions that can be used to visualize bed files and annotate TADs"""
# own libaries
from .classifier import Classifier

# third party libraries
import numpy as np
import pathlib
import pathlib
import pandas as pd
from sklearn.metrics import confusion_matrix
from itertools import permutations

# plotting
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from matplotlib.ticker import StrMethodFormatter


def plot_confusion_matrix(y_true, y_pred, classes,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax


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
    ax.set_xlabel('Size')
    if save:
        plt.savefig(output / 'size_histogram.png')
    else:
        plt.show()


def plot_annotation_dist(beds, annotation, output='./', plot_type='h', save=False):
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
    excpetions = ['None', 'NA']
    annotations = [float(bed.data[annotation])
                   for bed in beds if not bed.data[annotation] in excpetions]
    # calculate exact quartiles and other statistics
    df = pd.DataFrame({annotation: annotations})
    if plot_type == 'h':
        ax = sns.distplot(annotations, kde=True)
    elif plot_type == 'b':
        ax = sns.boxplot(y=annotations)
    elif plot_type == 'v':
        ax = sns.violinplot(y=annotations)
    ax.set_title(f'Distribution of {annotation} values')
    ax.set_xlabel(f'{annotation}')
    if save:
        output = pathlib.Path(output)
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
    ax = sns.boxplot(x=np.arange(0, length), y=elements)
    return ax


def plot_corr(df):
    '''Function plots a graphical correlation matrix for each pair of columns in the dataframe.

    Input:
        df: pandas DataFrame
        size: vertical and horizontal size of the plot'''

    corr = df.corr()
    sns.heatmap(corr,
                xticklabels=corr.columns.values,
                yticklabels=corr.columns.values)
    plt.tight_layout()
    plt.show()


def plot_multiple_roc(classifiers: [Classifier], test_sets, save=False, output=''):
    """Plots roc curve for multiple classifier."""
    plt.figure(figsize=(12, 10))
    plt.plot([0, 1], [0, 1], linestyle='--', label='random classification')
    plt.ylabel('TPR')
    plt.xlabel('FPR')
    plt.title(f'ROC curve')
    for idx, classifier in enumerate(classifiers):
        fpr, tpr = classifier.test(test_sets[idx])
        plt.plot(fpr, tpr, marker='.', label=classifier.name)
    plt.legend()
    if save:
        plt.savefig(pathlib.Path(output))
    else:
        plt.show()


def plot_avg_prec_scores(scores, k, support=[], k_name='Allele Count Threshold', save=False, output=''):
    """Plots a line plot with average precision means for each fiven K"""
    plt.figure(figsize=(12, 10))

    fig, ax1 = plt.subplots()

    color = 'tab:blue'
    ax1.plot(k, scores, color=color)
    ax1.set_ylabel('10-fold average precision mean', color=color)
    ax1.set_xlabel(k_name)
    ax1.set_xticks(k)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()

    color = 'tab:red'
    ax2.plot(k, support, color=color)
    ax2.set_ylabel('Number of test variants', color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()
    if save:
        plt.savefig(pathlib.Path(output) / 'Avg_Prec_per_AF.png')
    else:
        plt.show()


def plot_correlation_node_graph(cor):
    """Plots a node graph where nodes correspond to features and edges to correlations
    Args:
        cor = pandasDataFrame with pairwise correlation values.
    Output:
        fig = matplotlib figure of the node graph.
    """
    # Create edge dictionary
    orig, dest, cors = zip(*[(cor.columns[row], cor.columns[col], cor.iloc[row][cor.columns[col]])
                             for row, col in permutations(range(0, cor.shape[1]), 2) if cor.iloc[row][cor.columns[col]] >= 0.1])

    # replace  spaces with newline in node descriptors
    orig = [ori.replace(' ', '\n') for ori in orig]
    dest = [des.replace(' ', '\n') for des in dest]

    # Build a dataframe for the edges
    edge_df = pd.DataFrame({'from': orig, 'to': dest})

    # Build the node graph
    G = nx.from_pandas_edgelist(edge_df, 'from', 'to')

    # Graph with Custom nodes:
    fig = plt.figure(figsize=(12, 10))
    pos = nx.spring_layout(G, scale=0.5)
    nx.draw(G, pos=pos, with_labels=True, node_size=10000, node_color="#98db84", node_shape="o",
            alpha=0.7, linewidths=2, font_size=20, font_weight='bold', font_color='#000000')

    edge_labels = {}
    for idx, cor in enumerate(cors):
        key = (orig[idx], dest[idx])
        edge_labels[key] = round(cor, 4)
    nx.draw_networkx_edge_labels(
        G, pos, edge_labels=edge_labels, font_color='black', font_size=20, label_pos=0.5)
    return fig


def plot_feature_dist(df: pd.DataFrame, exclude_features=[]):
    """Plot the distribution of all columns in a pandas dataframe expect those sepcified in exlude_features
    Args:
        df = pandas DataFrame
        exclude_features = list of feature names i.e. strings
    Output:
        ax = matplotlib axis
    """
    df.drop(columns=exclude_features,inplace=True)
    ax = df.astype('float64').hist(bins=25,grid=False,color='#659666',zorder=2,rwidth=0.9,figsize=(12,10))
    for axis in ax:
        for x in axis:
            # Despine
            x.spines['right'].set_visible(False)
            x.spines['top'].set_visible(False)

            # Switch off ticks and change label size
            x.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="on", left="off", right="off", labelleft="on", labelsize=6)

            # Draw horizontal axis lines
            vals = x.get_yticks()
            for tick in vals:
                x.axhline(y=tick, linestyle='solid', alpha=0.4, color='#eeeeee', zorder=1)

            # Modify title
            x.title.set_fontsize(10)
            x.title.set_fontweight('bold')

            # Format y-axis label
            x.yaxis.set_major_formatter(StrMethodFormatter('{x:,g}'))
    return ax
