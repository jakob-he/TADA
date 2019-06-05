"""Preprocessing of annotated CNVs for ML models."""
#third part libaries
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder

# scikit classifier
from lib.classifier import Classifier

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
    cnv_features = []
    if feature_type=='binary':
        features = ['Boundary Breaking','Gene Overlap','Enhancer Overlap','TAD with high pLI','TAD with high Phastcon']
        for chrom in cnv_dict:
            for cnv in cnv_dict[chrom]:
                cnv_features.append(cnv.get_binary_features())

    feature_df = pd.DataFrame(data=cnv_features,columns=features)
    return feature_df


def create_stratified_training_and_test_set(cnv_dict_1,cnv_dict_2,feature_type='binary',oneHot=False,validation=False,exclude_features=[]):
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

    # create training and test set stratified by class labels
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=42, stratify=Y)

    if validation:
        # create validation and training set
        X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.2, random_state=42, stratify=y_train)

        return ([X_train,y_train],[X_val,y_val],[X_test,y_test])
    else:
        return ([X_train,y_train],[X_test,y_test])


def main():
    # set random seed
    np.random.seed(42)

    # parse inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('-c1','--cnvs1', help='Non pathogenic pickeled CNV set.')
    parser.add_argument('-c2','--cnvs2', help='Pathogenic pickeled CNV set.')
    args = parser.parse_args()

    # unpickle the annotated CNV files
    with pathlib.Path(args.cnvs1).open('rb') as non_pathogenic_cnvs:
        non_patho_cnvs = pickle.load(non_pathogenic_cnvs)

    with pathlib.Path(args.cnvs2).open('rb') as pathogenic_cnvs:
        patho_cnvs = pickle.load(pathogenic_cnvs)


    # create test and training set
    train_set_1, test_set_1 = create_stratified_training_and_test_set(non_patho_cnvs,patho_cnvs,oneHot=False)
    lr_1 = Classifier(solver='liblinear',name='All Binary Features',class_weight='balanced')
    lr_1.feature_selection(train_set_1)
    lr_1.train(train_set_1)

    # train second classifier only on overlap
    train_set_2, test_set_2 = create_stratified_training_and_test_set(non_patho_cnvs,patho_cnvs,oneHot=False,exclude_features=['Boundary Breaking','TAD with high pLI','TAD with high Phastcon'])
    lr_2= Classifier(solver='liblinear',name='Gene + Enhancer Overlap',class_weight='balanced')
    lr_2.feature_selection(train_set_2)
    lr_2.train(train_set_2)

    # train second classifier only TAD features
    train_set_3, test_set_3 = create_stratified_training_and_test_set(non_patho_cnvs,patho_cnvs,oneHot=False,exclude_features=['Boundary Breaking','Enhancer Overlap','Gene Overlap'])
    lr_3= Classifier(solver='liblinear',name='TAD Features',class_weight='balanced')
    lr_3.feature_selection(train_set_3)
    lr_3.train(train_set_3)

    # plot roc curve
    plotting.plot_multiple_roc([lr_1,lr_2,lr_3],[test_set_1,test_set_2,test_set_3],save=True,output='binary_features_ROC_cuves.png')

    # tensorflow logistic regression
    #lr = LR(train_set,val_set,test_set,0.0005,200)
    #lr.train()



if __name__ == '__main__':
    main()
