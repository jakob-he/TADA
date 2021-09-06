#!/usr/bin/python3

import argparse
import pathlib
import pandas as pd

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.calibration import calibration_curve


def argparser():
    parser = argparse.ArgumentParser('Generate ROC-Curves for one or multiple tab seperated txt files containing actual labels and predictions.')
    parser.add_argument('-f','--files', nargs='+', help='Paths to the prediction txt files.')
    parser.add_argument('-l','--labels', nargs='+', help='Labels for the data of the corresponding txt file. Should be in the same order as the files argument.')
    parser.add_argument('-o','--output_dir',help='Output dir for the generated figures.')
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
    plt.xlabel('False Positive Rate',fontsize=15)
    plt.ylabel('True Positive Rate',fontsize=15)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend(fontsize=15)
    plt.savefig(pathlib.Path(args.output_dir) / 'ROC_Curves.png')

    plt.rcdefaults()
    plt.figure(figsize=(12, 10))
    ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
    ax2 = plt.subplot2grid((3, 1), (2, 0))
    ax1.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")
    for idx,file in enumerate(args.files):
        file_path = pathlib.Path(file)
        file_df = pd.read_csv(file_path,header=None,index_col=False,sep='\t',names=['label','prediction'])

        # calculate roc curve
        fraction_of_positives, mean_predicted_value = calibration_curve(file_df['label'], file_df['prediction'], n_bins=10)

        ax1.plot(mean_predicted_value, fraction_of_positives, "s-",label=f'')
        ax2.hist(file_df['prediction'], range=(0, 1), bins=10, label=args.labels[idx].replace("_"," "),histtype="step", lw=2)

    ax1.set_ylabel("Fraction of positives")
    ax1.set_ylim([-0.05, 1.05])
    ax1.legend(loc="lower right")
    ax1.set_title('Calibration plots  (reliability curve)')

    ax2.set_xlabel("Mean predicted value")
    ax2.set_ylabel("Count")
    ax2.legend(loc="upper center", ncol=2)

    plt.tight_layout()
    plt.savefig(pathlib.Path(args.output_dir) / f'rf_calibration.png')









if __name__ == "__main__":
    main()
