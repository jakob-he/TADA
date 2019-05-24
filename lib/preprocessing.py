"""Preprocessing of annotated CNVs for ML models."""
import pandas as pd

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
    """Splits the merged feature dataframe of two CNV sets into a training set stratified by label"""
