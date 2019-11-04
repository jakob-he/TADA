"""implementation of a basic classifier class which allows to use different algorithms in the same framework."""
# classifiers
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from imblearn.ensemble import BalancedRandomForestClassifier

# preprocessing
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score

# Feature selection
from sklearn.decomposition import PCA

# visualization
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
import seaborn as sns

# metrics
from sklearn.metrics import f1_score
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score

# third party libararies
import numpy as np
import pathlib
import pandas as pd
import pickle

#own libraries
from . import utils
from . import plotting

classifier_dict = {'lr':LogisticRegression,'rf':RandomForestClassifier,'brf':BalancedRandomForestClassifier}

class Classifier():
    def __init__(self, classifier='lr', name='',**kwargs):
        """Initialization of a new Classifier object. This requires the name of classifier."""
        if name:
            self.name = name
        else:
            self.name = classifier
        self.classes = ['non-pathogenic','pathogenic']
        self.pipeline = Pipeline(steps=([('imputer', SimpleImputer()), ('scaler', MinMaxScaler(feature_range=(-1,1))), ('cls', classifier_dict[classifier](**kwargs))]))
        self.trained = False

    def feature_selection(self,train_set,test_set,feature_type,save=False,output_dir='./'):
        """Visulaizes the correlation in between the independent feateure and the class labels."""
        # merge training and test data
        merged_df = pd.concat([train_set[0],test_set[0]])

        merged_labels = pd.concat([train_set[1],test_set[1]])

        # compute Pearson correlation between features
        if 'binary' in feature_type:
            cor = merged_df.corr(utils.phi_coeff)
        else:
            cor = merged_df.corr()

        plt.figure(figsize=(12,12))
        sns.heatmap(cor, annot=True, cbar=False, cmap=plt.cm.Greens,annot_kws={"size": 20})
        plt.tick_params(labelsize=15)
        plt.yticks(np.arange(train_set[0].shape[1])+0.5, va="center")
        plt.tight_layout()
        plt.savefig(pathlib.Path(output_dir) / f'{self.name}_feature_correlations.png')

        #compute the partial correlation for all features
        P_corr = utils.partial_corr(merged_df,feature_type)
        plt.figure(figsize=(12,12))
        sns.heatmap(P_corr, annot=True, cbar=False,cmap=plt.cm.Greens, annot_kws={"size": 20})
        plt.tick_params(labelsize=15)
        plt.yticks(np.arange(train_set[0].shape[1])+0.5, va="center")
        plt.tight_layout()
        plt.savefig(pathlib.Path(output_dir) / f'{self.name}_feature_partial_correlations.png')

        #plot partial correlation as node graph
        fig = plotting.plot_correlation_node_graph(P_corr)
        fig.savefig(pathlib.Path(output_dir) / f'{self.name}_correlation_node_graph.png',bbox_inches='tight')

        # PCA
        plt.figure(figsize=(12,10))
        pca = PCA(n_components=2)
        projected = pca.fit_transform(merged_df)
        colors = ['red' if label else 'green' for label in merged_labels]
        plt.scatter(projected[:, 0], projected[:, 1],
                    c=colors, edgecolor='none', alpha=0.05)
        plt.xlabel('PC 1')
        plt.ylabel('PC 2')
        plt.yticks(np.arange(train_set[0].shape[1])+0.5, va="center")
        plt.tight_layout()
        plt.savefig(pathlib.Path(output_dir) / f'{self.name}_PCA.png')


    def train(self,train_set,output_dir='./'):
        skf = StratifiedKFold(n_splits=10)
        print('training with 10-fold CV...')

        # transform panda dataframes into np.arrays
        if any([type(set)==pd.DataFrame for set in train_set]):
            train_set = [set.values for set in train_set]

        average_precision = []

        average_precision = cross_val_score(self.pipeline, train_set[0], train_set[1], cv=skf, scoring='average_precision')
        print(f'10-fold CV Average Precision Mean: {np.mean(average_precision)}')

        self.pipeline.fit(train_set[0],train_set[1])
        self.trained = True
        # save model to pickle file
        with open(pathlib.Path(output_dir) / f'{self.name}_model.p','wb') as model_output:
            pickle.dump(self, model_output)

    def test(self,test_set,save=False,output_dir='./',plot=False):
        print('testing...')
        if self.trained:
            # predict class labels and probabilities
            y_pred = self.pipeline.predict(test_set[0])
            y_pred_scores = self.pipeline.predict_proba(test_set[0])[:, 1]

            #report classification metrics
            print(f'Average Precision on test-set: {average_precision_score(test_set[1],y_pred_scores)}')
            print(classification_report(test_set[1],y_pred,target_names=self.classes))

            # calculate roc curve
            precision, recall, thresholds = precision_recall_curve(test_set[1], y_pred_scores)

            if plot:
                # plot confusion matrix
                plotting.plot_confusion_matrix(test_set[1],y_pred,classes=['Non-pathogenic','Pathogenic'],cmap=plt.cm.Greens)
                if save:
                    plt.tight_layout()
                    plt.savefig(pathlib.Path(output_dir) / f'{self.name}_Confusion_Matrix.png')
                # plot feature importance
                # IMPOPRTANT! to be able to compare the coefficients all the features have to be on the same scale.
                if self.name == 'rf' or self.name == 'brf':
                    coef = self.pipeline['cls'].feature_importances_
                    ylabel = 'Feature Importance'
                else:
                    coef = self.pipeline['cls'].coef_[0]
                    ylabel = 'Model Coefficients'
                plt.figure(figsize=(12,10))
                feature_index = np.arange(len(coef))
                plt.bar(feature_index,coef,color='#2c7056')
                plt.ylabel(ylabel,fontsize=15)
                plt.tick_params(labelsize=15)
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
                plt.xlabel('Recall',fontsize=15)
                plt.ylabel('Precision',fontsize=15)
                plt.ylim([0.0, 1.05])
                plt.xlim([0.0, 1.0])
                plt.tick_params(labelsize=15)
                if save:
                    plt.tight_layout()
                    plt.savefig(pathlib.Path(output_dir) / f'{self.name}_ROC_Curve.png')
            return precision, recall
        else:
            return y_pred_scores
