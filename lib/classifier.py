"""implementation of a basic classifier class which allows to use different algorithms in the same framework."""
# classifiers
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from imblearn.ensemble import BalancedRandomForestClassifier

# CV
from sklearn.model_selection import StratifiedKFold

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
    def __init__(self, scaler, imputer, classifier='lr', name='',**kwargs):
        """Initialization of a new Classifier object. This requires the name of classifier."""
        if name:
            self.name = name
        else:
            self.name = classifier

        self.clf = classifier_dict[classifier](**kwargs)
        self.classes = ['non-pathogenic','pathogenic']
        self.trained = False
        self.scaler = scaler
        self.imputer = imputer

    def feature_selection(self,train_set,feature_type,save=False,output_dir='./'):
        """Visulaizes the correlation in between the independent feateure and the class labels."""
        # merge independent and dependent variables
        merged_df = train_set[0].copy()
        merged_df['label'] = train_set[1].copy()

        # compute Pearson correlation between features
        plt.figure(figsize=(12,10))

        if 'binary' in feature_type:
            cor = merged_df.corr(utils.phi_coeff)
        else:
            cor = merged_df.corr()

        sns.heatmap(cor, annot=True, cmap=plt.cm.Reds)
        plt.title('Pairwise Pearson Correlation Values')
        plt.savefig(pathlib.Path(output_dir) / f'{self.name}_feature_correlations.png')

        # compute correlation with the output variable
        cor_target = abs(cor['label'])

        # report highly correlated features
        relevant_features = cor_target[cor_target>0.5]
        print(f'Features with correlation larger than 0.5: {relevant_features}')

        # PCA
        plt.figure(figsize=(12,10))
        pca = PCA(n_components=2)
        projected = pca.fit_transform(merged_df.drop(['label'],axis=1))
        colors = ['red' if label else 'green' for label in merged_df['label']]
        plt.scatter(projected[:, 0], projected[:, 1],
                    c=colors, edgecolor='none', alpha=0.05)
        plt.xlabel('PC 1')
        plt.ylabel('PC 2')
        plt.savefig(pathlib.Path(output_dir) / f'{self.name}_PCA.png')


    def train(self,train_set,output_dir='./'):
        skf = StratifiedKFold(n_splits=10)
        print('training with 10-fold CV...')

        # transform panda dataframes into np.arrays
        if any([type(set)==pd.DataFrame for set in train_set]):
            train_set = [set.values for set in train_set]

        average_precision = []
        for train, test in skf.split(train_set[0],train_set[1]):
            self.clf = self.clf.fit(train_set[0][train],train_set[1][train])
            y_pred = self.clf.predict_proba(train_set[0][test])
            avg_prec_sc = average_precision_score(train_set[1][test],y_pred[:, 1])
            average_precision.append(avg_prec_sc)
            print(f'Average Precision: {avg_prec_sc}')

        print(f'10-fold CV Average Precision Mean: {np.mean(average_precision)}')
        self.trained = True
        # save model to pickle file
        with open(pathlib.Path(output_dir) / f'{self.name}_model.p','wb') as model_output:
            pickle.dump(self, model_output)

    def test(self,test_set,save=False,output_dir='./',plot=False):
        print('testing...')
        if self.trained:
            # predict class labels and probabilities
            y_pred = self.clf.predict(test_set[0])
            y_pred_scores = self.clf.predict_proba(test_set[0])[:, 1]

            #report classification metrics
            print(classification_report(test_set[1],y_pred,target_names=self.classes))

            # plot confusion matrix
            plotting.plot_confusion_matrix(test_set[1],y_pred,classes=['Non-pathogenic','Pathogenic'])
            if save:
                plt.tight_layout()
                plt.savefig(pathlib.Path(output_dir) / f'{self.name}_Confusion_Matrix.png')

            # calculate roc curve
            precision, recall, thresholds = precision_recall_curve(test_set[1], y_pred_scores)

            if plot:
                # plot feature importance
                # IMPOPRTANT! to be able to compare the coefficients all the features have to be on the same scale.
                if self.name == 'rf' or self.name == 'brf':
                    coef = self.clf.feature_importances_
                else:
                    coef = self.clf.coef_[0]
                plt.figure(figsize=(12,10))
                feature_index = np.arange(len(coef))
                plt.bar(feature_index,coef)
                plt.title('Coefficients')
                plt.xlabel('Features')
                plt.xticks(ticks=feature_index,labels=test_set[0].columns)
                if save:
                    plt.tight_layout()
                    plt.savefig(pathlib.Path(output_dir) / f'{self.name}_Coefficients.png')

                # plot the roc curve
                plt.figure(figsize=(12,10))
                plt.step(recall, precision, color='b', alpha=0.2, where='post')
                plt.fill_between(recall, precision, alpha=0.2, color='b', step='post')
                plt.xlabel('Recall')
                plt.ylabel('Precision')
                plt.ylim([0.0, 1.05])
                plt.xlim([0.0, 1.0])
                plt.title('2-class Precision-Recall curve')
                if save:
                    plt.tight_layout()
                    plt.savefig(pathlib.Path(output_dir) / f'{self.name}_ROC_Curve.png')
            return precision, recall
        else:
            return y_pred_scores
