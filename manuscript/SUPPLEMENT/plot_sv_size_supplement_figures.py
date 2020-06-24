"""
Plots Supplemental Figure 1 and 2: Size Distribution of Pathogenic and Non-Pathogenic Variants.
"""

import argparse
import pathlib
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

from scipy import stats
import numpy as np
from decimal import Decimal

import chromosome_plot

def argparser():
    parser = argparse.ArgumentParser('Plot the size distribution of Pathogenic and Non-Pathogenic Variants.')
    parser.add_argument('-b', '--bedfiles', nargs='+', help = 'Input bedfiles.')
    parser.add_argument('-t', '--title', help='Title of the figure')
    parser.add_argument('-i', '--ideogram', help='Ideogram of hg19.')
    return parser.parse_args()

def create_dfs(args):
    df_dict = {}
    # read variant data as pandas Dataframes
    # non pathogenic
    # before size matching
    non_patho_df_before_size_match_del = pd.read_csv(pathlib.Path(args.bedfiles[0]), header = None, index_col = False, sep = "\t", names = ['chrom','start','end','origin'])
    non_patho_df_before_size_match_dup = pd.read_csv(pathlib.Path(args.bedfiles[1]), header = None, index_col = False, sep = "\t", names = ['chrom','start','end','origin'])

    # after size matching
    non_patho_df_after_size_match_del = pd.read_csv(pathlib.Path(args.bedfiles[2]), header = None, index_col = False, sep = "\t", names = ['chrom','start','end','origin'])
    non_patho_df_after_size_match_dup = pd.read_csv(pathlib.Path(args.bedfiles[3]), header = None, index_col = False, sep = "\t", names = ['chrom','start','end','origin'])

    # pathogenic
    # before size matching
    patho_df_before_size_match_del = pd.read_csv(pathlib.Path(args.bedfiles[4]), header = None, index_col = False, sep = "\t", names = ['chrom','start','end','source','pathogenicity'])
    patho_df_before_size_match_dup = pd.read_csv(pathlib.Path(args.bedfiles[5]), header = None, index_col = False, sep = "\t", names = ['chrom','start','end','source','pathogenicity'])

    # split pathogenic by pathogenicity annotation
    patho_df_before_size_match_del_p = patho_df_before_size_match_del[patho_df_before_size_match_del['pathogenicity'] == 'Pathogenic']
    patho_df_before_size_match_del_l = patho_df_before_size_match_del[patho_df_before_size_match_del['pathogenicity'] == 'Likely']
    patho_df_before_size_match_del_u = patho_df_before_size_match_del[patho_df_before_size_match_del['pathogenicity'] == 'Unknown']

    patho_df_before_size_match_dup_p = patho_df_before_size_match_dup[patho_df_before_size_match_dup['pathogenicity'] == 'Pathogenic']
    patho_df_before_size_match_dup_l = patho_df_before_size_match_dup[patho_df_before_size_match_dup['pathogenicity'] == 'Likely']
    patho_df_before_size_match_dup_u = patho_df_before_size_match_dup[patho_df_before_size_match_dup['pathogenicity'] == 'Unknown']

    # after size matching
    patho_df_after_size_match_del = pd.read_csv(pathlib.Path(args.bedfiles[6]), header = None, index_col = False, sep = "\t", names = ['chrom','start','end','source','pathogenicity'])
    patho_df_after_size_match_dup = pd.read_csv(pathlib.Path(args.bedfiles[7]), header = None, index_col = False, sep = "\t", names = ['chrom','start','end','source','pathogenicity'])

    # split pathogenic by pathogenicity annotation
    patho_df_after_size_match_del_p = patho_df_after_size_match_del[patho_df_after_size_match_del['pathogenicity'] == 'Pathogenic']
    patho_df_after_size_match_del_l = patho_df_after_size_match_del[patho_df_after_size_match_del['pathogenicity'] == 'Likely']
    patho_df_after_size_match_del_u = patho_df_after_size_match_del[patho_df_after_size_match_del['pathogenicity'] == 'Unknown']

    patho_df_after_size_match_dup_p = patho_df_after_size_match_dup[patho_df_after_size_match_dup['pathogenicity'] == 'Pathogenic']
    patho_df_after_size_match_dup_l = patho_df_after_size_match_dup[patho_df_after_size_match_dup['pathogenicity'] == 'Likely']
    patho_df_after_size_match_dup_u = patho_df_after_size_match_dup[patho_df_after_size_match_dup['pathogenicity'] == 'Unknown']


    # add all dfs to the df_dict
    df_dict['non_patho_df_before_size_match_del'] = non_patho_df_before_size_match_del
    df_dict['non_patho_df_before_size_match_dup'] = non_patho_df_before_size_match_dup

    df_dict['non_patho_df_after_size_match_del'] = non_patho_df_after_size_match_del
    df_dict['non_patho_df_after_size_match_dup'] = non_patho_df_after_size_match_dup

    df_dict['patho_df_before_size_match_del'] = patho_df_before_size_match_del
    df_dict['patho_df_before_size_match_dup'] = patho_df_before_size_match_dup

    df_dict['patho_df_before_size_match_del_p'] = patho_df_before_size_match_del_p
    df_dict['patho_df_before_size_match_del_l'] = patho_df_before_size_match_del_l
    df_dict['patho_df_before_size_match_del_u'] = patho_df_before_size_match_del_u

    df_dict['patho_df_before_size_match_dup_p'] = patho_df_before_size_match_dup_p
    df_dict['patho_df_before_size_match_dup_l'] = patho_df_before_size_match_dup_l
    df_dict['patho_df_before_size_match_dup_u'] = patho_df_before_size_match_dup_u

    df_dict['patho_df_after_size_match_del'] = patho_df_after_size_match_del
    df_dict['patho_df_after_size_match_dup'] = patho_df_after_size_match_dup

    df_dict['patho_df_after_size_match_del_p'] = patho_df_after_size_match_del_p
    df_dict['patho_df_after_size_match_del_l'] = patho_df_after_size_match_del_l
    df_dict['patho_df_after_size_match_del_u'] = patho_df_after_size_match_del_u

    df_dict['patho_df_after_size_match_dup_p'] = patho_df_after_size_match_dup_p
    df_dict['patho_df_after_size_match_dup_l'] = patho_df_after_size_match_dup_l
    df_dict['patho_df_after_size_match_dup_u'] = patho_df_after_size_match_dup_u

    return df_dict

def plot_size_dist_box(df_1,df_2,xticks,ax):
    for df in [df_1,df_2]:
        df['size'] = df['end'] - df['start']

    size_1 = [np.log10(size) for size in df_1['size']]
    size_2 = [np.log10(size) for size in df_2['size']]

    # compute K-S p-values
    D, pvalue = stats.ks_2samp(df_1['size'],df_2['size'])

    ax.boxplot([size_1,size_2], showfliers=False)

    # include p-values into the plot
    x1, x2 = 1, 2
    y, h, col = max(np.percentile(size_1,75),np.percentile(size_2,75)), 1, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c=col)
    ax.text(1.5,y+h, f'K-S Statistic-{Decimal(D):.2e}', ha='center',fontsize='8', va='bottom', color=col)
    ax.set_xticklabels(['Non-Pathogenic','Pathogenic'])

    return ax

def ecdf(data: pd.Series):
    #create sorted array of unique elements
    uniques = np.sort(data.unique())
    #compute envenly spaced elements as x elements for the ecdf
    ecdf_x = np.linspace(start=min(uniques),stop=max(uniques),num=60)
    #get size of the raw values
    size = data.shape[0]
    #compute y values for the ecdf
    ecdf_y = []
    for x in ecdf_x:
        #count raw values below or equal to x
        cum_amount = data[data <= x].shape[0]
        data = data[data >= x]
        #save the cummulative amount of values as y values
        ecdf_y.append(cum_amount)
    return ecdf_x.tolist(),ecdf_y

def plot_size_histogram(df_1,df_2,ax):
    df_1['size'] = df_1['size'].apply(np.log10)
    df_2['size'] = df_2['size'].apply(np.log10)
    ecdf_x,ecdf_y = ecdf(df_1['size'])
    ax.hist(x='size',density=False,histtype='stepfilled',bins=ecdf_x,data=df_1[df_1['origin'] == 'GnomAD'],label='GnomAD', alpha=0.8)
    ax.hist(x='size',density=False,histtype='stepfilled',bins=ecdf_x,data=df_1[df_1['origin'] == 'Eichler'],label='Eichler', alpha=0.8)
    ax.hist(x='size',density=False,histtype='stepfilled',bins=ecdf_x,data=df_1[df_1['origin'] == 'Biobank'],label='Biobank', alpha=0.25)
    ax.hist(x='size',density=False,histtype='stepfilled',bins=ecdf_x,data=df_1[df_1['origin'] == 'DGV'],label='DGV', alpha=0.25)
    ax.hist(x='size',density=False,histtype='stepfilled',bins=ecdf_x,data=df_2,label='Decipher', alpha=0.8)
    return ax

def plot_pie_chart(df_1,column,ax):
    counts = df_1[column].value_counts()
    size = 0.2
    counts_percentages = [count/sum(counts) for count in counts]
    labels = [label for label in counts.index]
    ax.pie(counts_percentages,radius=0.8,labels=labels,autopct='%1.1f%%',wedgeprops=dict(width=size, edgecolor='w'))
    #ax.axis('equal')
    return ax

def plot_figure_1(df_dict):
    # Figure 1: Size Distributions and Proportions by Origin for Deletions
    fig = plt.figure(constrained_layout=False, figsize=(15,15))
    spec = gridspec.GridSpec(ncols=2, nrows=3, figure=fig)
    f_ax1 = fig.add_subplot(spec[0, 0])
    f_ax2 = fig.add_subplot(spec[0, 1])
    f_ax3 = fig.add_subplot(spec[1, 0])
    f_ax4 = fig.add_subplot(spec[1, 1])
    f_ax5 = fig.add_subplot(spec[2,:])

    # Plot A: Size Distribution of Deletions before Size-Matching
    f_ax1 = plot_size_dist_box(df_dict['non_patho_df_before_size_match_del'],df_dict['patho_df_before_size_match_del'], ['Non-Pathogenic','Pathogenic'], f_ax1)
    f_ax1.set_title('A')
    # Plot B: Size Distrubition of Duplications before Size-Matching
    f_ax2 = plot_size_dist_box(df_dict['non_patho_df_before_size_match_dup'],df_dict['patho_df_before_size_match_dup'], ['Non-Pathogenic','Pathogenic'], f_ax2)
    f_ax2.set_title('B')
    # Plot C: Size Distribution of Deletions before Size-Matching
    f_ax3 = plot_size_dist_box(df_dict['non_patho_df_after_size_match_del'],df_dict['patho_df_after_size_match_del'], ['Non-Pathogenic','Pathogenic'], f_ax3)
    f_ax3.set_title('C')
    # Plot D: Size Distrubition of Duplications before Size-Matching
    f_ax4 = plot_size_dist_box(df_dict['non_patho_df_after_size_match_dup'],df_dict['patho_df_after_size_match_dup'], ['Non-Pathogenic','Pathogenic'], f_ax4)
    f_ax4.set_title('D')
    # Plot E: Size Distribution Histogram
    f_ax5 = plot_size_histogram(df_dict['non_patho_df_before_size_match_del'],df_dict['patho_df_before_size_match_del'],f_ax5)
    f_ax5.set_ylabel('Variant Count')
    f_ax5.set_xlabel('Size log10(bp)')
    f_ax5.legend()
    f_ax5.set_title('E')

    #fig.text(0.5, 0.04, 'common X', ha='center')
    plt.gcf().text(0.065, 0.66, 'CNV Length log10(bp)', va='center', rotation='vertical')
    plt.subplots_adjust(left=0.1)
    return fig

def plot_figure_2(df_dict,ideo):
    # Figure 2: Proportion of Variants by source before and after size_matching
    fig = plt.figure(constrained_layout=False, figsize=(15,15))
    spec = gridspec.GridSpec(ncols=4, nrows=3, figure=fig)
    f_ax1 = fig.add_subplot(spec[0, 0])
    f_ax2 = fig.add_subplot(spec[0, 1])
    f_ax3 = fig.add_subplot(spec[0, 2])
    f_ax4 = fig.add_subplot(spec[0, 3])
    f_ax5 = fig.add_subplot(spec[1, 0])
    f_ax6 = fig.add_subplot(spec[1, 1])
    f_ax7 = fig.add_subplot(spec[1, 2])
    f_ax8 = fig.add_subplot(spec[1, 3])
    f_ax9 = fig.add_subplot(spec[2, :])


    # Plot A: Proportion of Non-Pathogenic Variants by Source after and before size-matching
    f_ax1 = plot_pie_chart(df_dict['non_patho_df_before_size_match_del'],'origin', f_ax1)
    f_ax1.set_title('A')
    f_ax2 = plot_pie_chart(df_dict['non_patho_df_after_size_match_del'],'origin', f_ax2)
    f_ax2.set_title('B')
    f_ax3 = plot_pie_chart(df_dict['non_patho_df_before_size_match_dup'],'origin', f_ax3)
    f_ax3.set_title('C')
    f_ax4 = plot_pie_chart(df_dict['non_patho_df_after_size_match_dup'],'origin', f_ax4)
    f_ax4.set_title('D')

    # Plot B: Proportion of Pathogenic Variants by Pathogenicity after and before size-matching
    f_ax5 = plot_pie_chart(df_dict['patho_df_before_size_match_del'],'pathogenicity', f_ax5)
    f_ax5.set_title('E')
    f_ax6 = plot_pie_chart(df_dict['patho_df_after_size_match_del'],'pathogenicity', f_ax6)
    f_ax6.set_title('F')
    f_ax7 = plot_pie_chart(df_dict['patho_df_before_size_match_dup'],'pathogenicity', f_ax7)
    f_ax7.set_title('G')
    f_ax8 = plot_pie_chart(df_dict['patho_df_after_size_match_dup'],'pathogenicity', f_ax8)
    f_ax8.set_title('H')

    f_ax9 = chromosome_plot.plot_chromsome_distribution(ideo,df_dict['non_patho_df_before_size_match_del'], f_ax9)
    f_ax9.set_title('I')
    #fig.text(0.5, 0.04, 'common X', ha='center')
    # plt.gcf().text(0.065, 0.75, 'A', va='center', rotation='horizontal')
    # plt.gcf().text(0.065, 0.25, 'B', va='center', rotation='horizontal')
    # plt.subplots_adjust(left=0.1)

    return fig






def main():
    # read cli
    args = argparser()

    # get variant dfs
    df_dict = create_dfs(args)

    # plot Figure 1
    fig_1 = plot_figure_1(df_dict)
    plt.savefig('SV_sizes.png',bbox_inches='tight')

    # Read in ideogram.txt, downloaded from UCSC Table Browser
    ideo = pd.read_csv(
        pathlib.Path(args.ideogram),
        skiprows=1, sep='\t',
        names=['chrom', 'start', 'end', 'name', 'gieStain']
    )
    # plot Figure 2
    fig_2 = plot_figure_2(df_dict, ideo)
    plt.savefig('SV_proportions.png',bbox_inches='tight')

if __name__ == '__main__':
    main()
