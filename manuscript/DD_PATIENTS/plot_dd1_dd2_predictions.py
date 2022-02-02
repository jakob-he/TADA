"""Plot the pathogenicity scores for DD2 duplications with the spiked in pathogenic DD1 duplication."""
import argparse
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def argparser():
    parser = argparse.ArgumentParser(
        description="Annotate a set of CNVs. Run annotate_cnvs -h for more details.")
    parser.add_argument('-v', '--variant_preds',
                        help='Path to the CNV file containing variants and predictions.', required=True)
    parser.add_argument('-o', '--output', help='Output location.', required=True)
    return parser.parse_args()


def main():
    args = argparser()

    cnv_df = pd.read_csv(pathlib.Path(args.variant_preds), header=0, sep="\t")

    cnv_df['pathogenic'] = np.select([(cnv_df['START'] == 67442001) | (cnv_df['START'] == 67958880)], [True], default=False)
    
    print(cnv_df)

    non_pathogenic = cnv_df[cnv_df['pathogenic'] == False]
    pathogenic = cnv_df[cnv_df['pathogenic'] == True]
    print(non_pathogenic)
    print(pathogenic)

    colors = ['#A82B0B', '#0B88A8']
    patients = ['DD1', 'DD2']
    fig, ax = plt.subplots(figsize=(12,10))
    ax.scatter(non_pathogenic.index,non_pathogenic['Pathogenicity Score'], color=colors[1], s=20)
    ax.scatter(pathogenic.index,
            pathogenic['Pathogenicity Score'], color=colors[0], s=20)
    ax.plot(cnv_df.index, np.repeat(0.5, cnv_df.shape[0]), color='lightgrey', linestyle='dashed', label='Pathogenicity Threshold (0.5)')
    for idx, row in pathogenic.iterrows():
        print(idx)
        ax.text(idx,row['Pathogenicity Score'],patients.pop(), fontsize=15)
    ax.set_ylim(0,1)
    ax.set_yticks(np.linspace(0,1,11))
    ax.set_yticklabels([str(round(bin,1)) for bin in np.linspace(0, 1, 11)],fontsize=15)
    ax.set_ylabel('Pathogenicity Score', fontsize=15)
    ax.set_xlabel('Variants', fontsize=15)
    ax.set_xticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.legend(loc='upper right', fontsize=15)
    fig.savefig(pathlib.Path(args.output) / 'Pathogenicity_DD2_DUPs.png', dpi=300)
    


if __name__ == "__main__":
    main()
