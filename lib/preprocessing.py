"""Preprocessing of annotated CNVs for ML models."""
#third part libaries
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder, RobustScaler, MaxAbsScaler

# tensorflow classifier
#from lib.lr import LR

#testing libararies
import argparse
import pathlib
import pickle
import lib.plotting as plotting

def create_feature_df(cnv_dict,feature_type):
    """Creates a pandas Dataframe containing cnvs as rows and features as columns"""
    #get features for each CNV
    if feature_type=='basic_binary':
        features = ['Boundary Breaking','Gene Overlap','Enhancer Overlap']
    if feature_type=='extended_binary':
        features = ['Boundary Breaking','Gene Overlap','Enhancer Overlap','DDG2P Genes Overlap','CTCF','TADs with high pLI','TADs with high Phastcon']
    if feature_type=='basic_continous':
        features = ['Boundary Distance','Gene Distance','Enhancer Distance']
    if feature_type=='extended_continous':
        features = ['Boundary Distance','Gene Distance','Enhancer Distance','DDG2P distance','gene pLI','enhancer conservation','HI score gene','CTCF distance']

    cnv_features = []
    for chrom in cnv_dict:
        for cnv in cnv_dict[chrom]:
            if cnv.tads:
                cnv_features.append(cnv.annotate(feature_type))

    feature_df = pd.DataFrame(data=cnv_features,columns=features)
    return feature_df


def create_stratified_training_and_test_set(cnv_dict_1,cnv_dict_2,feature_type,oneHot=False,validation=False,exclude_features=[]):
    """Splits the merged feature dataframe of two CNV sets into a training set stratified by label."""
    # create feature dataframes for both cnv dicts
    df_0 = create_feature_df(cnv_dict_1,feature_type)
    df_1 = create_feature_df(cnv_dict_2,feature_type)

    # exclude features
    for feature in exclude_features:
        df_0.drop([feature],axis=1,inplace=True)
        df_1.drop([feature],axis=1,inplace=True)

    # add labels to the dataframe. The first cnv dict is considered to be class 0.
    df_0['label'] = np.repeat(0,df_0.shape[0])
    df_1['label'] = np.repeat(1,df_1.shape[0])

    # concatenate dataframe
    df_merged = pd.concat([df_0,df_1],ignore_index=True)

    # sort by all values
    # df_merged.sort_values(by=df_merged.columns, inplace=True,ascending=False)
    # plotting.plot_corr(df_merged)


    # define X and y
    X = df_merged.loc[: , df_merged.columns != 'label']

    if oneHot:
        X = df_merged.loc[: , df_merged.columns != 'label'].values
        Y = df_merged['label'].values.reshape(-1,1)
        # One hot encode all features and labels
        oneHot = OneHotEncoder(categories='auto')

        oneHot.fit(X)
        X = oneHot.transform(X).toarray()

        oneHot.fit(Y)
        Y = oneHot.transform(Y).toarray()
    else:
        Y = df_merged['label']

    if 'continous' in feature_type:
        # replace NA value with max of the column
        X.replace({column:{-1:X[column].mean()} for column in X.columns},inplace=True)

    # create training and test set stratified by class labels
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=42, stratify=Y)

    # scale continous features.
    #TODO: Adapt this for validation sets.
    if 'continous' in feature_type:
        scaler = MinMaxScaler(feature_range=(-1,1))
        transformer = scaler.fit(X_train)
        X_train = pd.DataFrame(transformer.transform(X_train),columns=X_train.columns)
        X_test = pd.DataFrame(transformer.transform(X_test),columns=X_test.columns)

    if validation:
        # create validation and training set
        X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.2, random_state=42, stratify=y_train)

        return ([X_train,y_train],[X_val,y_val],[X_test,y_test])
    else:
        return ([X_train,y_train],[X_test,y_test])
