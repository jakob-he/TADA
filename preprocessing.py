"""Preprocessing of annotated CNVs for ML models."""
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from lib.lr import LR

#testing libararies
import argparse
import pathlib
import pickle
import lib.plotting as plotting

def create_feature_df(cnv_dict):
    """Creates a pandas Dataframe containing cnvs as rows and features as columns"""
    #get features for each CNV
    cnv_features = []
    for chrom in cnv_dict:
        for cnv in cnv_dict[chrom]:
            cnv_features.append(cnv.get_features())

    feature_df = pd.DataFrame(data=cnv_features,columns=['Gene Overlap','Enhancer Overlap','TAD with high pLI','TAD with high Phastcon'])
    return feature_df


def create_stratified_training_and_test_set(cnv_dict_1,cnv_dict_2):
    """Splits the merged feature dataframe of two CNV sets into a training set stratified by label."""
    # create feature dataframes for both cnv dicts
    df_0 = create_feature_df(cnv_dict_1)
    df_1 = create_feature_df(cnv_dict_2)

    # add labels to the dataframe. The first cnv dict is considered to be class 0.
    df_0['label'] = np.repeat(0,df_0.shape[0])
    df_1['label'] = np.repeat(1,df_1.shape[0])

    # concatenate dataframe
    df_merged = pd.concat([df_0,df_1],ignore_index=True)

    # sort values
    df_merged.sort_values(['Gene Overlap','Enhancer Overlap','TAD with high pLI','TAD with high Phastcon'], inplace=True,ascending=False)
    # plotting.plot_corr(df_merged)

    # define X and y
    X = df_merged.loc[: , df_merged.columns != 'label'].values
    Y = df_merged['label'].values.reshape(-1,1)

    # One hot encode all features and labels
    oneHot = OneHotEncoder(categories='auto')

    oneHot.fit(X)
    x = oneHot.transform(X).toarray()

    oneHot.fit(Y)
    y = oneHot.transform(Y).toarray()

    # create training and test set stratified by class labels
    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42, stratify=y)

    # create validation and training set
    X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.2, random_state=42, stratify=y_train)

    return ([X_train,y_train],[X_val,y_val],[X_test,y_test])



def main():
    # parse inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('-c1','--cnvs1',help='Non pathogenic pickeled CNV set.')
    parser.add_argument('-c2','--cnvs2',help='Pathogenic pickeled CNV set.')
    args = parser.parse_args()

    # unpickle the annotated CNV files
    with pathlib.Path(args.cnvs1).open('rb') as non_pathogenic_cnvs:
        non_patho_cnvs = pickle.load(non_pathogenic_cnvs)

    with pathlib.Path(args.cnvs2).open('rb') as pathogenic_cnvs:
        patho_cnvs = pickle.load(pathogenic_cnvs)


    # create test and training set
    train_set, val_set, test_set = create_stratified_training_and_test_set(non_patho_cnvs,patho_cnvs)

    lr = LR(train_set,val_set,test_set,0.0005,200)
    lr.train()







if __name__ == '__main__':
    main()
