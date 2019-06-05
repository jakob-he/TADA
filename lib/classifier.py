"""implementation of a basic classifier class which allows to use different algorithms in the same framework."""
# classifiers
from sklearn.linear_model import LogisticRegression

#CV
from sklearn.model_selection import StratifiedKFold

# visualization
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
import seaborn as sns

# metrics
from sklearn.metrics import f1_score
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

import pathlib
import pandas as pd

#own libraries
from . import utils

classifier_dict = {'lr':LogisticRegression}

class Classifier():
    def __init__(self, classifier='lr',name='',**kwargs):
        """Initialization of a new Classifier object. This requires the name of classifier."""
        if name:
            self.name = name
        else:
            self.name = classifier

        self.clf = classifier_dict[classifier](**kwargs)
        self.classes = ['non-pathogenic','pathogenic']
        self.trained = False

    def feature_selection(self,train_set,save=False,output_dir='./'):
        """Visulaizes the correlation in between the independent feateure and the class labels."""
        # merge independent and dependent variables
        merged_df = train_set[0].copy()
        merged_df['label'] = train_set[1].copy()

        # compute Pearson correlation between features
        plt.figure(figsize=(12,10))
        cor = merged_df.corr(utils.phi_coeff)
        sns.heatmap(cor, annot=True, cmap=plt.cm.Reds)
        plt.title('Pariwise Pearson Correlation Values')
        plt.savefig(pathlib.Path(output_dir) / f'{self.name}_feature_correlations.png')

        # compute correlation with the output variable
        cor_target = abs(cor['label'])

        # report highly correlated features
        relevant_features = cor_target[cor_target>0.5]
        print(f'Features with correlation larger than 0.5: {relevant_features}')

    def train(self,train_set):
        skf = StratifiedKFold(n_splits=10)
        print('training with 10-fold CV...')

        # transform panda dataframes into np.arrays
        if any([type(set)==pd.DataFrame for set in train_set]):
            train_set = [set.values for set in train_set]

        for train, test in skf.split(train_set[0],train_set[1]):
            self.clf = self.clf.fit(train_set[0][train],train_set[1][train])
            y_pred = self.clf.predict_proba(train_set[0][test])
            print(f'AUC: {roc_auc_score(train_set[1][test],y_pred[:, 1])}')

        self.trained = True

    def test(self,test_set,save=False,output_dir='./',plotting=False):
        print('testing...')
        if self.trained:
            # predict class labels and probabilities
            y_pred = self.clf.predict(test_set[0])
            y_pred_scores = self.clf.predict_proba(test_set[0])[:, 1]
            print(y_pred_scores)

            #report classification metrics
            print(classification_report(test_set[1],y_pred,target_names=self.classes))

            # calculate roc curve
            fpr, tpr, thresholds = roc_curve(test_set[1], y_pred_scores)

            if plotting:
                # plot the roc curve
                plt.figure(figsize=(12,10))
                plt.plot([0, 1], [0, 1], linestyle='--',label='random classification')
                plt.plot(fpr, tpr, marker='.')
                plt.ylabel('TPR')
                plt.xlabel('FPR')
                plt.title(f'{self.name} ROC curve')
                plt.legend()
                if save:
                    plt.savefig(pathlib.Path(output_dir) / f'{self.classifier_name}_ROC_Curve.png')
            return fpr, tpr
