"""implementation of a basic classifier class which allows to use different algorithms in the same framework."""
# classifiers
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

# preprocessing
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score

# Feature/model selection
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from scipy.cluster import hierarchy
#from sklearn.feature_selection import RFECV

# visualization
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
import seaborn as sns

# metrics
from sklearn.metrics import f1_score
from sklearn.metrics import r2_score
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score

# third party libararies
import numpy as np
import pathlib
import pandas as pd
import pickle
from collections import defaultdict

#own libraries
from . import utils
from . import plotting

classifier_dict = {'lr':LogisticRegression,'rf':RandomForestClassifier}

gridcv_parameters_dict = {'lr':
                   {'cls__solver': 'liblinear',
                    'cls__C': np.arange(0,1,0.1),
                    'cls__penalty': ["l1","l2"]},
                   'rf':
                   {'cls__max_depth': [50, 100, None],
                    'cls__max_features': ['auto', 'log2'],
                    'cls__min_samples_leaf': [1, 3, 5],
                    'cls__min_samples_split': [2, 4, 6, 10],
                    'cls__n_estimators': [100, 300, 500],
                    'cls__oob_score' : True}}

tuned_parameters_dict = {'lr':
                   {'cls__solver': 'liblinear',
                    'cls__C': 0.1,
                    'cls__penalty': "l1"},
                   'rf':
                   {'cls__max_depth': 50,
                    'cls__max_features': 'auto',
                    'cls__min_samples_leaf': 1,
                    'cls__min_samples_split': 2,
                    'cls__n_estimators': 500,
                    'cls__oob_score' : True}}

class Classifier():
    def __init__(self, classifier='lr', name='',**kwargs):
        """Initialization of a new Classifier object. This requires the name of classifier."""
        if name:
            self.name = name
        else:
            self.name = classifier
        self.classes = ['non-pathogenic','pathogenic']
        self.pipeline = Pipeline(steps=([('imputer', SimpleImputer(missing_values=np.nan, strategy='mean')), ('scaler', MinMaxScaler()), ('cls', classifier_dict[classifier](**kwargs))]))
        self.trained = False

    def feature_selection(self,train_set,feature_type,save=False,output_dir='./'):
        """Visualize the correlation in between the features."""
        labels = train_set[0].columns

        # scale and impute the training data
        scaled_train_df = self.pipeline[0:2].fit_transform(train_set[0])

        output_path = pathlib.Path(output_dir)

        # compute Pearson correlation between features
        # TODO: update phi_coeff for numpy array
        # if 'binary' in feature_type:
        #     cor = scaled_train_df.corr(utils.phi_coeff)
        # else:
        cor = np.corrcoef(scaled_train_df,rowvar=False)

        fig = plotting.plot_correlation_heatmap(cor,labels)
        plt.savefig(output_path / f'{self.name}_feature_correlations.png')

        #compute the partial correlation for all features
        P_corr = utils.partial_corr(scaled_train_df,feature_type)

        # plot partial correlation heatmap
        fig = plotting.plot_correlation_heatmap(P_corr, labels)
        plt.savefig(output_path / f'{self.name}_feature_partial_correlations.png')

        #plot partial correlation as node graph
        P_corr = pd.DataFrame(P_corr,columns=labels)
        P_corr.rename({idx:label for idx,label in enumerate(labels)},inplace=True)
        fig = plotting.plot_correlation_node_graph(P_corr)
        fig.savefig(output_path / f'{self.name}_correlation_node_graph.png',bbox_inches='tight')

        # plot correlation dendrogram based on hierachical clustering
        plt.figure(figsize=(10,12))
        cor_linkage = hierarchy.ward(cor)
        dendro = hierarchy.dendrogram(cor_linkage, labels=train_set[0].columns,
                                      leaf_rotation=90)
        plt.tight_layout()
        plt.savefig(output_path / f'{self.name}_correlation_dendogram.png')

        # PCA
        # plt.figure(figsize=(12,10))
        # pca = PCA(n_components=2)
        # projected = pca.fit_transform(scaled_train_df)
        # colors = ['red' if label else 'green' for label in merged_labels]
        # plt.scatter(projected[:, 0], projected[:, 1],
        #             c=colors, edgecolor='none', alpha=0.05)
        # plt.xlabel('PC 1')
        # plt.ylabel('PC 2')
        # plt.yticks(np.arange(train_set[0].shape[1])+0.5, va="center")
        # plt.tight_layout()
        # plt.savefig(pathlib.Path(output_dir) / f'{self.name}_PCA.png')


    def train(self,train_set,output_dir='./',permut_importance=True,gridcv=False):
        skf = StratifiedKFold(n_splits=5)

        # extract feature names
        names = train_set[0].columns

        # transform panda dataframes into np.arrays
        if any([type(set)==pd.DataFrame for set in train_set]):
            train_set = [set.values for set in train_set]


        if gridcv:
            # Initialize GridSearch object
            gscv = GridSearchCV(self.pipeline, parameters_dict[self.name], cv = skf, verbose = 1, scoring = "roc_auc")
            # Fit gscv
            print(f"Tuning {self.name}...")
            gscv.fit(train_set[0], train_set[1])

            # Get best parameters and score
            print(f'Best parameters: {gscv.best_params_}')
            print(f'Best ROC-AUC score {gscv.best_score_}')
        else:
            # set classifier params
            self.pipeline.set_params(**tuned_parameters_dict[self.name])
            # compute cross validation score
            roc_auc_scores = cross_val_score(self.pipeline, train_set[0], train_set[1], cv=skf, scoring='roc_auc')
            print(f'5-fold CV ROC-AUC Mean: {np.mean(roc_auc_scores)}')
            self.pipeline.fit(train_set[0],train_set[1])

        self.trained = True

        # save model to pickle file
        with open(pathlib.Path(output_dir) / f'{self.name}_model.p','wb') as model_output:
            pickle.dump(self, model_output)

        # print feature importance by decreased accuracy
        if permut_importance:
            # scale and impute training data
            scaled_train_df = self.pipeline[0:2].transform(train_set[0])

            # compute pearson correlation and cluster feature with a distance <= 1
            cor = np.corrcoef(scaled_train_df,rowvar=False)
            corr_linkage = hierarchy.ward(cor)
            cluster_ids = hierarchy.fcluster(corr_linkage, 1, criterion='distance')

            # save clustered feature indices
            cluster_id_to_feature_ids = defaultdict(list)
            for idx, cluster_id in enumerate(cluster_ids):
                cluster_id_to_feature_ids[cluster_id].append(idx)

            # iterate through clusters, permutating the corresponding features
            scores = defaultdict()
            acc = utils.oob_classifier_accuracy(self.pipeline, train_set[0], train_set[1])
            for cluster, feature_idx in cluster_id_to_feature_ids.items():
                X_t = train_set[0].copy()
                X_t[:,feature_idx] = np.random.permutation(X_t[:,feature_idx])
                shuff_acc = utils.oob_classifier_accuracy(self.pipeline, X_t, train_set[1])
                scores[cluster] = acc-shuff_acc

            # sort the scores and plot them in decending order
            sorted_scores, sorted_clusters = list(zip(*sorted([(round(score, 4), cluster) for cluster, score in scores.items()], reverse=False)))
            feature_idx = [cluster_id_to_feature_ids[cluster_id] for cluster_id in sorted_clusters]
            feature_names = [[names[idx] for idx in feature_list] for feature_list in feature_idx]
            feature_ticks = ['\n'.join(features) for features in feature_names]
            heights = [len(features) for features in feature_names]
            coordinates = np.arange(len(feature_ticks))
            # # get plotting coordinates according to the size of each cluster
            # coordinate = 0
            # coordinates = []
            # for idx, features in enumerate(feature_names):
            #     if idx != 0:
            #         coordinate += len(features)
            #         coordinates.append(coordinate)
            #     else:
            #         coordinates.append(coordinate)
            plt.figure(figsize=(12,10))
            plt.barh(coordinates,sorted_scores,color='#007878')
            plt.yticks(coordinates,feature_ticks)
            plt.tight_layout()
            plt.savefig(pathlib.Path(output_dir) / 'feature_importance_acc_decrease.png')

        # # plot roc auc as a function of added features by importance
        # features = []
        # roc_auc_scores = []
        # for scores, feat in sorted_scores:
        #     features.append(np.argwhere(names == feat).flat[0])
        #     X_t = train_set[0][:,features]
        #     roc_auc_scores.append(cross_val_score(self.pipeline, X_t, train_set[1], cv=skf, scoring='roc_auc'))
        # plt.figure(figsize=(12,10))
        # plt.plot(coordinates,[np.mean(roc_auc_score) for roc_auc_score in roc_auc_scores],color='#007878')
        # plt.xticks(coordinates,[feat for scores, feat in sorted_scores],rotation=90)
        # plt.tight_layout()
        # plt.savefig(pathlib.Path(output_dir) / 'roc_auc_over_features.png')

    def test(self,test_set,save=False,output_dir='./',plot=False):
        print('testing...')
        if self.trained:
            # predict class labels and probabilities
            y_pred = self.pipeline.predict(test_set[0])
            y_pred_scores = self.pipeline.predict_proba(test_set[0])[:, 1]

            #report classification metrics
            print(f'ROC AUC on test-set: {roc_auc_score(test_set[1],y_pred_scores)}')
            print(classification_report(test_set[1],y_pred,target_names=self.classes))

            # calculate roc curve
            precision, recall, thresholds = precision_recall_curve(test_set[1], y_pred_scores)

            if plot:
                # plot confusion matrix
                plotting.plot_confusion_matrix(test_set[1],y_pred,classes=['Non-pathogenic','Pathogenic'],cmap=plt.cm.BuGn)
                if save:
                    plt.tight_layout()
                    plt.savefig(pathlib.Path(output_dir) / f'{self.name}_Confusion_Matrix.png')
                # plot feature importance
                # IMPOPRTANT! to be able to compare the coefficients all the features have to be on the same scale.
                if self.name == 'rf':
                    coef = self.pipeline['cls'].feature_importances_
                    ylabel = 'Feature Importance'
                else:
                    coef = self.pipeline['cls'].coef_[0]
                    ylabel = 'Model Coefficients'
                plt.figure(figsize=(12,10))
                feature_index = np.arange(len(coef))
                plt.bar(feature_index,coef,color='#007878')
                plt.ylabel(ylabel,fontsize=25)
                plt.tick_params(labelsize=25)
                if len(coef) < 5:
                    plt.xticks(ticks=feature_index,labels=test_set[0].columns)
                else:
                    plt.xticks(ticks=feature_index,labels=test_set[0].columns, rotation='vertical')
                if save:
                    plt.tight_layout()
                    plt.savefig(pathlib.Path(output_dir) / f'{self.name}_Coefficients.png')

                # plot the roc curve
                plt.figure(figsize=(12,10))
                plt.step(recall, precision, color='#2c7056', alpha=0.2, where='post')
                plt.fill_between(recall, precision, alpha=0.2, color='#2c7056', step='post')
                plt.xlabel('Recall',fontsize=25)
                plt.ylabel('Precision',fontsize=25)
                plt.ylim([0.0, 1.05])
                plt.xlim([0.0, 1.0])
                plt.tick_params(labelsize=25)
                if save:
                    plt.tight_layout()
                    plt.savefig(pathlib.Path(output_dir) / f'{self.name}_PRC.png')
            return precision, recall
        else:
            return y_pred_scores
