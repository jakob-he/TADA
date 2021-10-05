"""implementation of a basic classifier class which allows to use different algorithms in the same framework."""
# classifiers
from sklearn.svm import OneClassSVM
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

# preprocessing
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score

# Feature/model selection
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from scipy.cluster import hierarchy
import shap
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
from sklearn.calibration import calibration_curve

# third party libararies
import numpy as np
import pathlib
import pandas as pd
import pickle
from collections import defaultdict

# own libraries
from . import utils
from . import plotting

classifier_dict = {'lr': LogisticRegression, 'rf': RandomForestClassifier, 'ossvm':OneClassSVM}

gridcv_parameters_dict = {'lr':
                          {'cls__solver': 'liblinear',
                           'cls__C': np.arange(0, 1, 0.1),
                           'cls__penalty': ["l1", "l2"]},
                          'rf':
                          {'cls__max_depth': [50, 100, None],
                              'cls__max_features': ['auto', 'log2'],
                              'cls__min_samples_leaf': [1, 3, 5],
                              'cls__min_samples_split': [2, 4, 6, 10],
                              'cls__n_estimators': [300, 500, 1000],
                              'cls__oob_score': [True]}}

tuned_parameters_dict = {'lr':
                         {'cls__solver': 'liblinear',
                          'cls__C': 0.1,
                          'cls__penalty': "l1"},
                         'rf':
                         {'cls__max_depth': 50,
                             'cls__max_features': 'auto',
                             'cls__min_samples_leaf': 5,
                             'cls__min_samples_split': 4,
                             'cls__n_estimators': 500,
                             'cls__oob_score': True}}


class Classifier():
    def __init__(self, classifier, name='', one_class=False, **kwargs):
        """Initialization of a new Classifier object. This requires the name of classifier."""
        self.one_class = one_class
        if name:
            self.name = name
        else:
            self.name = classifier
        self.classes = ['non-pathogenic', 'pathogenic']
        self.pipeline = Pipeline(steps=([('imputer', SimpleImputer(missing_values=np.nan, strategy='mean')), (
            'scaler', StandardScaler()), ('cls', classifier_dict[classifier](**kwargs))]))
        self.trained = False

    def feature_selection(self, train_set, feature_type, save=False, output_dir='./'):
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
        cor = np.corrcoef(scaled_train_df, rowvar=False)

        fig = plotting.plot_correlation_heatmap(cor, labels)
        plt.savefig(output_path / f'{self.name}_feature_correlations.png')

        # compute the partial correlation for all features
        P_corr = utils.partial_corr(scaled_train_df, feature_type)

        # plot partial correlation heatmap
        fig = plotting.plot_correlation_heatmap(P_corr, labels)
        plt.savefig(
            output_path / f'{self.name}_feature_partial_correlations.png')

        # plot partial correlation as node graph
        P_corr = pd.DataFrame(P_corr, columns=labels)
        P_corr.rename(
            {idx: label for idx, label in enumerate(labels)}, inplace=True)
        fig = plotting.plot_correlation_node_graph(P_corr)
        fig.savefig(
            output_path / f'{self.name}_correlation_node_graph.png', bbox_inches='tight')

        # plot correlation dendrogram based on hierachical clustering
        plt.figure(figsize=(10, 12))
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

    def train(self, train_set, output_dir='./', permut_importance=True, gridcv=False):
        if self.one_class:
            self.train_one_class(train_set, output_dir, gridcv)
        else:
            skf = StratifiedKFold(n_splits=5)

            # extract feature names
            names = train_set[0].columns

            # transform panda dataframes into np.arrays
            if any([type(set) == pd.DataFrame for set in train_set]):
                train_set = [set.values for set in train_set]

            if gridcv:
                # Initialize GridSearch object
                gscv = GridSearchCV(
                    self.pipeline, gridcv_parameters_dict[self.name], n_jobs=3, cv=skf, verbose=1, scoring="roc_auc")
                # Fit gscv
                print(f"Tuning {self.name}...")
                gscv.fit(train_set[0], train_set[1])

                # Get best parameters and score
                print(f'Best parameters: {gscv.best_params_}')
                print(f'Best ROC-AUC score {gscv.best_score_}')

                # train pipeline on all data with the best parameters
                self.pipeline.set_params(**gscv.best_params_)
                self.pipeline.fit(train_set[0], train_set[1])
            else:
                # set classifier params
                self.pipeline.set_params(**tuned_parameters_dict[self.name])
                # compute cross validation score
                roc_auc_scores = cross_val_score(
                    self.pipeline, train_set[0], train_set[1], cv=skf, scoring='roc_auc')
                print(f'5-fold CV ROC-AUC Mean: {np.mean(roc_auc_scores)}')
                self.pipeline.fit(train_set[0], train_set[1])

            # print feature importance by decreased accuracy
            if permut_importance:
                # scale and impute training data
                scaled_train_df = self.pipeline[0:2].transform(train_set[0])

                # compute pearson correlation and cluster feature with a distance <= 1
                cor = np.corrcoef(scaled_train_df, rowvar=False)
                corr_linkage = hierarchy.ward(cor)
                cluster_ids = hierarchy.fcluster(
                    corr_linkage, 1, criterion='distance')

                # save clustered feature indices
                cluster_id_to_feature_ids = defaultdict(list)
                for idx, cluster_id in enumerate(cluster_ids):
                    cluster_id_to_feature_ids[cluster_id].append(idx)

                # iterate through clusters, permutating the corresponding features
                acc = utils.oob_classifier_accuracy(
                    self.pipeline, train_set[0], train_set[1])
                scores = defaultdict(list)
                for i in np.random.choice(100, 30, replace=False):
                    np.random.seed(i)
                    for cluster, feature_idx in cluster_id_to_feature_ids.items():
                        X_t = train_set[0].copy()
                        X_t[:, feature_idx] = np.random.permutation(
                            X_t[:, feature_idx])
                        shuff_acc = utils.oob_classifier_accuracy(
                            self.pipeline, X_t, train_set[1])
                        scores[cluster].append(acc - shuff_acc)
                    mean_scores = {cluster:[np.mean(score),np.std(score)] for cluster, score in scores.items()}

                # sort the scores and plot them in decending order
                sorted_scores, sorted_clusters, sorted_std = list(zip(
                    *sorted([(round(score[0], 4), cluster, round(score[1], 4)) for cluster, score in mean_scores.items()], reverse=False)))
                feature_idx = [cluster_id_to_feature_ids[cluster_id]
                            for cluster_id in sorted_clusters]
                feature_names = [[names[idx] for idx in feature_list]
                                for feature_list in feature_idx]
                feature_ticks = ['\n'.join(features) for features in feature_names]
                heights = [len(features) for features in feature_names]
                coordinates = np.arange(len(feature_ticks))
                fig = plt.figure(figsize=(12, 10))
                plt.barh(coordinates, sorted_scores, xerr=sorted_std, color='#007878')
                plt.yticks(coordinates, feature_ticks)
                plt.xlabel('Mean-Decrease in Accuracy')
                plt.gca().spines['top'].set_visible(False)
                plt.gca().spines['right'].set_visible(False)
                plt.tight_layout()
                plt.savefig(pathlib.Path(output_dir) /
                            'feature_importance_acc_decrease.png')
        
    def train_one_class(self, train_set, output_dir='./', gridcv=False):
        # transform panda dataframes into np.arrays
        # if any([type(set) == pd.DataFrame for set in train_set]):
        #     train_set = [set.values for set in train_set]
        
        if gridcv:
            # Initialize GridSearch object
            gscv = GridSearchCV(
                self.pipeline, gridcv_parameters_dict[self.name], n_jobs=3, cv=5, verbose=2, scoring="neg_mean_squared_error")

            # Fit gscv
            print(f"Tuning {self.name}...")
            gscv.fit(train_set)

            # Get best parameters and score
            print(f'Best parameters: {gscv.best_params_}')
            print(f'Best Score (MSE) {gscv.best_score_}')

            # train pipeline on all data with the best parameters
            self.pipeline.set_params(**gscv.best_params_)
            self.pipeline.fit(train_set)
        else:
            self.pipeline.fit(train_set)

        self.trained = True

        # save model to pickle file
        with open(pathlib.Path(output_dir) / f'{self.name}_model.p', 'wb') as model_output:
            pickle.dump(self, model_output, protocol=4)

    def test(self, test_set, save=False, output_dir='./', plot=False):
        print('testing...')
        if self.trained:
            # predict class labels and probabilities
            y_pred = self.pipeline.predict(test_set[0])
            y_pred_scores = self.pipeline.predict_proba(test_set[0])[:, 1]

            if save:
                with open(pathlib.Path(output_dir) /'ROC_scores.txt','w') as test_output:
                    roc_scores = [f"{test_set[1].values[idx]}\t{pred}" for idx, pred in enumerate(y_pred_scores)]
                    test_output.write("\n".join(roc_scores))

            # report classification metrics
            print(
                f'ROC AUC on test-set: {roc_auc_score(test_set[1],y_pred_scores)}')
            print(classification_report(
                test_set[1], y_pred, target_names=self.classes))

            # calculate roc curve
            fpr, tpr, thresholds = roc_curve(test_set[1], y_pred_scores)
            auc = roc_auc_score(test_set[1], y_pred_scores)

            if plot:
                # plot confusion matrix
                plotting.plot_confusion_matrix(test_set[1], y_pred, classes=[
                                               'Non-pathogenic', 'Pathogenic'], cmap=plt.cm.BuGn)
                if save:
                    plt.tight_layout()
                    plt.savefig(pathlib.Path(output_dir) /
                                f'{self.name}_Confusion_Matrix.png')
                # plot feature importance
                if self.name == 'rf':
                    coef = self.pipeline['cls'].feature_importances_
                    ylabel = 'Feature Importance'
                else:
                    coef = self.pipeline['cls'].coef_[0]
                    ylabel = 'Model Coefficients'
                plt.figure(figsize=(12, 10))
                feature_index = np.arange(len(coef))
                plt.bar(feature_index, coef, color='#007878')
                plt.ylabel(ylabel, fontsize=25)
                plt.tick_params(labelsize=25)
                if len(coef) < 5:
                    plt.xticks(ticks=feature_index, labels=test_set[0].columns)
                else:
                    plt.xticks(ticks=feature_index,
                               labels=test_set[0].columns, rotation='vertical')
                if save:
                    plt.tight_layout()
                    plt.savefig(pathlib.Path(output_dir) /
                                f'{self.name}_Coefficients.png')

                # plot the roc curve
                plt.figure(figsize=(12, 10))
                plt.plot(fpr, tpr, color='#007878', label='Random Forest')
                plt.plot([0, 1], [0, 1], color='grey', lw=2, linestyle='--')
                # axis labels
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                # show the legend
                plt.legend()
                if save:
                    plt.tight_layout()
                    plt.savefig(pathlib.Path(output_dir) /
                                f'{self.name}_ROC.png')

                # plot calibration curve
                plt.rcdefaults()
                plt.figure(figsize=(12, 10))
                ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
                ax2 = plt.subplot2grid((3, 1), (2, 0))
                ax1.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")
                fraction_of_positives, mean_predicted_value = calibration_curve(test_set[1], y_pred_scores, n_bins=10)

                ax1.plot(mean_predicted_value, fraction_of_positives, "s-",label=f'{self.name}')
                ax2.hist(y_pred_scores, range=(0, 1), bins=10, label=f'{self.name}',histtype="step", lw=2)

                ax1.set_ylabel("Fraction of positives")
                ax1.set_ylim([-0.05, 1.05])
                ax1.legend(loc="lower right")
                ax1.set_title('Calibration plots  (reliability curve)')

                ax2.set_xlabel("Mean predicted value")
                ax2.set_ylabel("Count")
                ax2.legend(loc="upper center", ncol=2)

                plt.tight_layout()
                plt.savefig(pathlib.Path(output_dir) / f'{self.name}_calibration.png')

    def test_one_class(self, test_set, train_set, labels, plot=False, output_dir='.'):
        print('testing...')
        if self.trained:
            # predict class labels and probabilities
            y_pred_train = self.pipeline.predict(train_set)
            y_pred_test = self.pipeline.predict(test_set)

            train_accuracy = 100*list(y_pred_train
                                ).count(-1)/(len(y_pred_train))
            print("Outliers in Training-Set:", train_accuracy)

            test_accuracy = 100*list(y_pred_test
                                ).count(-1)/(len(y_pred_test))
            print("Outliers in Test-Set:", test_accuracy)
       # if plot:
            # tsne_3 = TSNE(n_components=3)
            # tsne_2 = TSNE(n_components=2)
            # train_scaled = self.pipeline[:2].transform(train_set)
            # test_scaled = self.pipeline[:2].transform(test_set)
            
            # # for idx, feature in enumerate(labels):
            # #     plt.figure(figsize=(12,10))
            # #     shap.dependence_plot(
            # #         idx, shap_values, shap_test, feature_names=labels, show=False)
            # #     plt.savefig(pathlib.Path(output_dir) /
            # #             f'SHAP_dependency_plot_{feature}.png', bbox_inches='tight', dpi=300)

            # test_set['anomaly_score'] = self.pipeline[2].decision_function(test_scaled)
            # test_set['anomaly_score_decision'] = y_pred_test
            # test_set.to_csv(pathlib.Path(output_dir) /
            #                 'predicted_test_decipher.csv', header=True, sep='\t')

            # merged = np.concatenate([train_scaled, test_scaled])
            # reduced_2 = tsne_2.fit_transform(merged)
            # reduced_train_2 = reduced_2[:len(train_scaled)]
            # reduced_test_2 = reduced_2[len(train_scaled):]

            # reduced_3 = tsne_3.fit_transform(merged)
            # reduced_train_3 = reduced_3[:len(train_scaled)]
            # reduced_test_3 = reduced_3[len(train_scaled):]

            # # plot dimensionality 3D reduced data 
            # fig = plt.figure(figsize=(12, 10))
            # ax = fig.add_subplot(111, projection='3d')
            # ax.set_xlabel("t-SNE dimension 1")
            # ax.set_ylabel("t-SNE dimension 2")
            # ax.set_zlabel("t-SNE dimension 3")

            # # Plot the compressed data points
            # ax.scatter(reduced_train_3[:, 0], reduced_train_3[:, 1], zs=reduced_train_3[:, 2], s=4, lw=1, label="train", c="green")

            # # Plot the compressed data points
            # ax.scatter(reduced_test_3[:, 0], reduced_test_3[:, 1], zs=reduced_test_3[:, 2],  lw=2, s=60, marker="x", label="test", c="red")
            # ax.legend()
            # fig.savefig(pathlib.Path(output_dir) / 'TSNE_train_test_3D.png', dpi=300)

            # # plot dimensionality 2D reduced data and
            # fig = plt.figure(figsize=(12, 10))
            # ax = fig.add_subplot(111)
            # ax.set_xlabel("t-SNE dimension 1")
            # ax.set_ylabel("t-SNE dimension 2")

            # # Plot the compressed data points
            # ax.scatter(reduced_train_2[:, 0], reduced_train_2[:, 1], s=4, lw=1, label="train", c="green")

            # # Plot the compressed data points
            # ax.scatter(reduced_test_2[:, 0], reduced_test_2[:, 1], lw=2, s=60, marker="x", label="test", c="red")
            # ax.legend()
            # fig.savefig(pathlib.Path(output_dir) /
            #             'TSNE_train_test_2D.png', dpi=300)

            # # plot SVM decision in 2D
            # # define the meshgrid
            # fig = plt.figure(figsize=(12,10))
            # x_min, x_max = reduced_train_2[:, 0].min() - 5, reduced_train_2[:, 0].max() + 5
            # y_min, y_max = reduced_train_2[:, 1].min() - 5, reduced_train_2[:, 1].max() + 5

            # x_ = np.linspace(x_min, x_max, 500)
            # y_ = np.linspace(y_min, y_max, 500)

            # xx, yy = np.meshgrid(x_, y_)

            # # evaluate the decision function on the meshgrid
            # svm = OneClassSVM(kernel='rbf', gamma='scale', nu=0.15)
            # svm.fit(reduced_train_2)
            # z = svm.decision_function(np.c_[xx.ravel(), yy.ravel()])
            # z = z.reshape(xx.shape)

            # # plot the decision function and the reduced data
            # plt.contourf(xx, yy, z, cmap=plt.cm.PuBu)
            # a = plt.contour(xx, yy, z, levels=[0], linewidths=2, colors='darkred')
            # b = plt.scatter(reduced_train_2[:,0], reduced_train_2[:,1], c='white', s=3)
            # c = plt.scatter(
            #     reduced_test_2[:, 0], reduced_test_2[:, 1], c='red', marker='X', s=4)
            # plt.legend([a.collections[0], b, c], ['Decision Boundary',
            #         'Common', 'Pathogenic'], bbox_to_anchor=(1.05, 1))
            # plt.axis('tight')
            # fig.savefig(pathlib.Path(output_dir) /
            #             'TSNE_SVM_decision_boundary.png', dpi=300)

            # shap_train = shap.kmeans(train_scaled, 50)
            # shap_test = shap.sample(test_scaled, 50)

            # explainer = shap.KernelExplainer(
            #     self.pipeline[2].predict, shap_train)
            # shap_values = explainer.shap_values(shap_test)
            # with open(pathlib.Path(output_dir) / 'shap_values.p', "wb") as shap_file:
            #     pickle.dump(shap_values, shap_file, protocol=4)
            # plt.figure(0, figsize=(12,12))
            # shap.summary_plot(shap_values, shap_test,
            #                 show=False, feature_names=labels)
            # plt.savefig(pathlib.Path(output_dir) /
            #             'SHAP_summary_plot.png', bbox_inches='tight', dpi=300)

            # # plt.figure(1, figsize=(12, 12))
            # shap.plots.bar(shap_values, shap_test, max_display=20, show=False)
            # plt.savefig(pathlib.Path(output_dir) /
            #             'SHAP_summary_plot_bar.png', bbox_inches='tight', dpi=300)

            # for idx in range(1,20):
            #     plt.figure(idx, figsize=(20, 15))
            #     shap.force_plot(explainer.expected_value, shap_values[idx],
            #                     show=False, matplotlib=True, feature_names=labels)
            #     plt.savefig(pathlib.Path(output_dir) /
            #                 f'SHAP_force_plot_sample_{idx}.png', bbox_inches='tight', dpi=300)
            #     plt.close()