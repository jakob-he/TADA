"""Preprocessing of annotated SVs for ML models."""
# third part libaries
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from sklearn.impute import MissingIndicator

# testing libararies
import argparse
import pathlib
import pickle
from . import plotting


def create_feature_df(sv_dict, feature_type, labels, csv=False):
    """Creates a pandas Dataframe containing svs as rows and features as columns"""
    # get features for each SV
    sv_features = []
    if csv:
        for chrom in sv_dict:
            for sv in sv_dict[chrom]:
                if sv.tads:
                    sv_features.append(
                        np.append([sv.chr, sv.start, sv.end], sv.annotate(feature_type)))
        feature_df = pd.DataFrame(data=sv_features, columns=[
                                  'CHR', 'START', 'END'] + labels)
    else:
        for chrom in sv_dict:
            for sv in sv_dict[chrom]:
                if sv.tads:
                    sv_features.append(sv.annotate(feature_type))

        feature_df = pd.DataFrame(data=sv_features, columns=labels)
    return feature_df


def create_stratified_training_and_test_set(sv_dicts, feature_type, labels, oneHot=False, feature_dfs = [], exclude_features=[], one_class=False):
    """Splits the merged feature dataframe of two SV sets into a training set stratified by label."""
    if one_class:
        df = create_feature_df(sv_dicts, feature_type, labels)
        return df
    else:
        # create feature dataframes for both sv dicts
        if feature_dfs:
            df_0 = feature_dfs[0]
            df_1 = feature_dfs[1]
        else:
            df_0 = create_feature_df(sv_dicts[0], feature_type, labels)
            df_1 = create_feature_df(sv_dicts[1], feature_type, labels)

        # exclude features
        for feature in exclude_features:
            df_0.drop([feature], axis=1, inplace=True)
            df_1.drop([feature], axis=1, inplace=True)

        # add labels to the dataframe. The first sv dict is considered to be class 0.
        df_0['label'] = np.repeat(0, df_0.shape[0])
        df_1['label'] = np.repeat(1, df_1.shape[0])

        # concatenate dataframe
        df_merged = pd.concat([df_0, df_1], ignore_index=True)

        # sort by all values
        # df_merged.sort_values(by=df_merged.columns, inplace=True,ascending=False)
        # plotting.plot_corr(df_merged)

        # define X and y
        X = df_merged.loc[:, df_merged.columns != 'label']
        
        if oneHot:
            X = df_merged.loc[:, df_merged.columns != 'label'].values
            Y = df_merged['label'].values.reshape(-1, 1)
            # One hot encode all features and labels
            oneHot = OneHotEncoder(categories='auto')

            oneHot.fit(X)
            X = oneHot.transform(X).toarray()

            oneHot.fit(Y)
            Y = oneHot.transform(Y).toarray()
        else:
            Y = df_merged['label']

        # Add Missing indicators to the feature dataframe
        # TODO: This is only valid for the etended feature set
        # if feature_type == 'extended_continuous':
        #     indicator = MissingIndicator(features='all')
        #     indicator = pd.DataFrame(indicator.fit_transform(X), columns=['Boundary Distance','Gene Distance Missing','Enhancer Distance Missing','DDG2P Distance Missing','Gene LOEUF','Enhancer conservation','Gene HI','CTCF Distance Missing','HI LogOdds Score'])
        #     X = pd.concat([X,indicator[['Gene Distance Missing','Enhancer Distance Missing','DDG2P Distance Missing','CTCF Distance Missing']]],axis=1)

        # create training and test set stratified by class labels
        X_train, X_test, y_train, y_test = train_test_split(
            X, Y, test_size=0.3, stratify=Y)

        return ([X_train, y_train], [X_test, y_test])


def create_stratified_training_and_test_set_one_class(sv_dicts, feature_type, labels, exclude_features=[], ):
    """Splits the merged feature dataframe of two SV sets into a training set stratified by label."""
    # create feature dataframes for the sv dict
    df_0 = create_feature_df(sv_dicts[0], feature_type, labels)
    df_1 = create_feature_df(sv_dicts[1], feature_type, labels)

    # exclude features
    for feature in exclude_features:
        df_0.drop([feature], axis=1, inplace=True)
        df_1.drop([feature], axis=1, inplace=True)

    # sort by all values
    # df_merged.sort_values(by=df_merged.columns, inplace=True,ascending=False)
    # plotting.plot_corr(df_merged)
    return (df_0, df_1)
