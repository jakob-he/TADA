#!/usr/bin/python3

import argparse
import pathlib
import pandas as pd

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score


def argparser():
    parser = argparse.ArgumentParser('Generate ROC-Curves for one or multiple tab seperated txt files containing actual labels and predictions.')
    parser.add_argument('-f','--files', nargs='+', help='Paths to the prediction txt files.')
    parser.add_argument('-l','--labels', nargs='+', help='Labels for the data of the corresponding txt file. Should be in the same order as the files argument.')
    parser.add_argument('-o','--output',help='Output path for the generated figure.')
    return parser.parse_args()


def main():
    # read cli
    args = argparser()

    plt.figure(figsize=(12, 10))
    for idx,file in enumerate(args.files):
        file_path = pathlib.Path(file)
        file_df = pd.read_csv(file_path,header=None,index_col=False,sep='\t',names=['label','prediction'])

        # calculate roc curve
        fpr, tpr, thresholds = roc_curve(file_df['label'], file_df['prediction'])
        auc = roc_auc_score(file_df['label'], file_df['prediction'])

        plt.plot(fpr, tpr, label=f'{args.labels[idx].replace("_"," ")} AUC = {round(auc,4)}')

    plt.plot([0, 1], [0, 1], color='grey', lw=2, linestyle='--', label='Random Classifier AUC = 0.5')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.tight_layout()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend(fontsize=15)
    plt.savefig(pathlib.Path(args.output))









if __name__ == "__main__":
    main()
