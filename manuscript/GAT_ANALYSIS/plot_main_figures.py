"""Plots ther anlysis results of GAT"""
import argparse
import pathlib
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from collections import OrderedDict


def argparser():
    parser = argparse.ArgumentParser('Plots l2fold of a GAT analysis with the qvalue.')
    parser.add_argument('-a','--pathogenic',help='Path to the pathogenic anlysis file.')
    parser.add_argument('-b','--nonpathogenic',help='Path to the non-pathogenic anlysis file.')
    parser.add_argument('-c','--compare',help='Path to the count comparision of pathogenic and non-pathogenic.')
    parser.add_argument('-o','--output',help='Figure title.')
    parser.add_argument('-q','--qvalue',type=float,help='Q-value threshold.')
    args = parser.parse_args()
    return args

def main():
    #get input arguments
    args = argparser()
    pathogenic_ana = pathlib.Path(args.pathogenic)
    nonpathogenic_ana = pathlib.Path(args.nonpathogenic)
    compare_ana = pathlib.Path(args.compare)

    #define dict for y axis ticks
    tick_dict = {'HMM_Active_Enhancer':'Activer Enhancer',
                 'HMM_Insulator':'Insulator',
                 'HMM_Poised_Promotor':'Poised Promotor',
                 'HMM_Strong_Enhancer':'Strong Enhancer',
                 'HMM_Transcriptional_elongation':'Transcriptional Elongation',
                 'HMM_Transcriptional_transition':'Transcriptional Transition',
                 'HMM_Weak_Enhancer':'Weak Enhancer',
                 'HMM_Weak_Promotor':'Weak Promotor',
                 'HMM_Weak_Repetitive':'Repetitiv',
                 'HMM_Weak_Transcribed':'Weak Transcribed',
                 'Crossover_Paternal_Cov_75':'Crossover Peaks Paternal 75%',
                 'Crossover_Paternal_Cov_90':'Crossover Peaks Paternal 90%',
                 'Crossover_Maternal_Cov_75':'Crossover Peaks Maternal 75%',
                 'Crossover_Maternal_Cov_90':'Crossover Peaks Maternal 90%',
                 'all_enhancer_75_conservation':r'Conserved Enhancer',
                 'all_enhancer_90_conservation':r'Highly Conserved Enhancer',
                 'fantom_5_all_enhancer':'Enhancer',
                 'fantom_5_brain_enhancer':'Brain',
                 'Phylop_Negative':'Negative PhyloP Enhancer',
                 'Phylop_Positive':'Positive PhyloP Enhancer',
                 'Phylop_Positive_75':r'Conserved Enhancer',
                 'Phylop_Positive_90':r'Enhancer with PhyloP $\geq$ 90 Percentile',
                 'TADS_high_conservation_enhancer':'TADS containing an enhancer with Phastcon > 0.9',
                 'TADs_high_mean_conservation':'TADS with mean Phastcon > 0.2',
                 'TADs_high_pLI_genes':'TADS containing a gene with pLI=1',
                 'TADs_high_mean_pLI':'TADS with mean pLI > 0.8',
                 'gnomAD_genes':'Genes',
                 'telomeres':'Telomeric Regions',
                 'TADs_with_phastcon_1_enhancer':'TADs with Phastcon = 1',
                 'TADs_with_phastcon_0.9_enhancer':'TADs with Phastcon >= 0.9',
                 'TADs_with_phastcon_0.5_enhancer':'TADs with Phastcon >= 0.5',
                 'TADs_with_phastcon_0.1_enhancer':'TADs with Phastcon >= 0.1',
                 'TADs_with_phastcon_0_enhancer':'TADs with Phastcon >= 0',
                 'TADs_with_pLI_1_genes':'TADs with pLI = 1',
                 'TADs_with_pLI_0.9_genes':'TADs with pLI >= 0.9',
                 'TADs_with_pLI_0.5_genes':'TADS with pLI >= 0.5',
                 'TADs_with_pLI_0.51_genes':'TADS with pLI >= 0.1',
                 'TADs_with_pLI_0_genes':'TADS with pLI >= 0',
                 'tads_without_functional_elements':'TADs without genes and enhancers',
                 'CTCF': 'CTCF Binding Sites',
                 'DDG2P_genes': 'DDG2P Genes',
                 'high_HI':r'HI Genes',
                 'low_HI':r'HS Genes',
                 'high_loeuf':r'pLoF Tolerant Genes',
                 'low_loeuf':r'pLoF Intolerant Genes',
                 'Segmental_duplications':'Segmental Duplications',
                 'TAD_boundaries':'TAD Boundaries',
                 'high_pLI':r'GnomAD Genes with pLI $\geq$ 0.9',
                 'low_pLI':r'GnomAD Genes with pLI $\leq$ 0.1',
                 }


    # load pathogenic data into a df
    cols = ['Telomeric Regions','TAD Boundaries','CTCF Binding Sites','Highly Conserved Enhancer','Conserved Enhancer','Enhancer','HI Genes','HS Genes','pLoF Intolerant Genes','pLoF Tolerant Genes','DDG2P Genes','Genes']
    annotation_dict = {annotation:idx for idx,annotation in enumerate(cols)}

    patho_df = pd.read_csv(pathogenic_ana,header=0,sep='\t')
    patho_df['annotation'] = [tick_dict[annotation] for annotation in patho_df['annotation']]
    patho_df['annotation'] = [annotation_dict[annotation] for annotation in patho_df['annotation']]
    patho_df.sort_values(by='annotation',inplace=True)
    patho_df.reset_index(inplace=True)
    patho_df['label'] = 'patho'

    # load pathogenic data into a df
    nonpatho_df = pd.read_csv(nonpathogenic_ana,header=0,sep='\t')
    nonpatho_df['annotation'] = [tick_dict[annotation] if annotation in tick_dict else annotation for annotation in nonpatho_df['annotation']]
    nonpatho_df['annotation'] = [annotation_dict[annotation] for annotation in nonpatho_df['annotation']]
    nonpatho_df.sort_values(by='annotation',inplace=True)
    nonpatho_df.reset_index(inplace=True)
    nonpatho_df['label'] = 'non_patho'



    merged_df = pd.concat([patho_df,nonpatho_df])
    merged_df.reset_index(inplace=True)


    # create figure
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, gridspec_kw = {'width_ratios':[5, 1]}, figsize=(15,12))

    # set barheight
    barheight = 0.40
    bar_padding_1 = 0.05
    bar_padding_2 = 0.1
    bars = np.arange(patho_df.shape[0])
    patho_bars = [bar + barheight for bar in bars]
    nonpatho_bars = [bar + barheight + bar_padding_1 for bar in patho_bars]
    bar_ticks = [bar + barheight/2 + bar_padding_1/2 for bar in patho_bars]


    f.subplots_adjust(hspace=0.025, wspace=0.05)

    # ax1.grid(b=True,axis='both',which='major', color='lightgrey', linewidth=1.0)
    ax1.barh(patho_bars, patho_df['l2fold'], height=barheight, color = '#913a1d', linewidth=1.5,edgecolor='#781e00')
    ax1.barh(nonpatho_bars, nonpatho_df['l2fold'], height=barheight, color = '#1d9191', linewidth=1.5,edgecolor='#007878')


    #sns.barplot(x="l2fold", y="annotation", data=df, color='#007878',ax=ax1)

    ax1.set_xlim(-1.5,1.5)
    ax1.set_ylim(0,13)
    ax1.set_ylabel('')
    ax1.set_xlabel('log2 FoldChange',fontsize=15)
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.set_ticks_position('none')
    ax1.axvline(0,linewidth=1, color='black')
    ax1.set_axisbelow(True)
    #ax1.axhline(nonpatho_bars[0] + barheight/2 + 0.05, xmax = 0.8, linewidth=1,color='#b9bab7')
    for bar in nonpatho_bars:
        ax1.axhline(bar + barheight/2 + 0.05,linewidth=1,color='#b9bab7')


    # read in file with compared counts
    compare_df = pd.read_csv(compare_ana,header=0,sep='\t',comment='#')
    compare_df['annotation'] = [tick_dict[annotation] if annotation in tick_dict else annotation for annotation in compare_df['annotation']]
    compare_df['annotation'] = [annotation_dict[annotation] for annotation in compare_df['annotation']]
    compare_df.sort_values(by='annotation',inplace=True)
    compare_df.reset_index(inplace=True)


    sizes = [min(200,abs(300*fold)) for fold in compare_df['observed']]

    for idx,size in enumerate(sizes):
        facecolors = 'black'
        alpha = 0.5
        if compare_df['qvalue'][idx] <= args.qvalue:
            alpha = 1.0
        ax2.scatter(x=0,y=bar_ticks[idx],marker='s',s=size,color='black',linewidths=2,alpha=alpha,facecolors=facecolors)

    ax2.set_frame_on(False) #Remove both axes
    ax2.get_yaxis().set_visible(False)
    ax2.get_xaxis().set_visible(False)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax2.set_title('Overlap Difference',fontsize=15)

    childrenLS=ax1.get_children()
    barlist=list(filter(lambda x: isinstance(x, matplotlib.patches.Rectangle), childrenLS))
    for idx, qvalue in enumerate(merged_df['qvalue']):
        if qvalue > args.qvalue:
            barlist[idx].set_color('#bdbebf')

    patho_patch = mpatches.Patch(color='#913a1d', label='Pathogenic')
    non_patho_patch = mpatches.Patch(color='#1d9191', label='Non-Pathogenic')
    legend = ax2.legend(by_label.values(), by_label.keys(),loc='center left',bbox_to_anchor=(1,0.5),framealpha=0.05,labelspacing=1)
    plt.setp(legend.get_title(),fontsize=15)
    ax1.set_yticks(bar_ticks)
    ax1.set_yticklabels(labels=cols)
    ax1.tick_params(labelsize=15)
    ax1.legend(handles=[patho_patch,non_patho_patch],loc='lower right',labelspacing=1,fontsize=15,bbox_to_anchor=(1,0.15),borderaxespad=0.7)
    f.savefig(f'{args.output}.png',bbox_inches='tight')


if __name__=='__main__':
    main()
