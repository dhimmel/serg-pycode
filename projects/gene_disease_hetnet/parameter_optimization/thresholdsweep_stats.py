import os
import gzip
import pprint
import csv
import collections

import rpy2
import rpy2.robjects
import rpy2.robjects.numpy2ri
import sklearn
import sklearn.metrics
import sklearn.linear_model
import numpy
import pandas
import pandas.io.parsers

#network_dir = '/home/dhimmels/Documents/serg/ashg13/140310-parsweep'
project_dir = '/home/dhimmels/Documents/serg/gene-disease-hetnet'
network_dir = os.path.join(project_dir, 'networks', '140518-thresholdsweep')



part_path = os.path.join(project_dir, 'partitions.txt.gz')
part_file = gzip.open(part_path)
part_rows = list(csv.DictReader(part_file, delimiter='\t'))
part_file.close()
#part_df = pandas.io.parsers.read_table(part_path, compression='gzip')
training_tuples = {(row['disease_code'], row['gene_symbol'])
                   for row in part_rows
                   if row['part'] == 'train'}

gene_counter = collections.Counter(row['gene_symbol'] for row in part_rows if row['status_int'] == '1')
disease_counter = collections.Counter(row['disease_code'] for row in part_rows if row['status_int'] == '1')


feature_path = os.path.join(network_dir, 'features-exp0.4.txt.gz')
feature_df = pandas.io.parsers.read_table(feature_path, compression='gzip')
is_training = feature_df.apply(lambda x: (x['target'], x['source']) in training_tuples, axis=1)
feature_df = feature_df.loc[is_training]

column_names = feature_df.columns.values.tolist()
feature_names = column_names[4:]

feature_df['PCs'] = feature_df.apply(lambda x: gene_counter[x['source']] - x['status'], axis=1)
feature_df['PCt'] = feature_df.apply(lambda x: disease_counter[x['target']] - x['status'], axis=1)

rpy2.robjects.numpy2ri.activate()
rpy2.robjects.r(
'''
library(glmnet)
RidgePredictions <- function(X, y, X_predict) {
  y <- as.vector(y)
  cv.ridge = glmnet::cv.glmnet(X, y, family='binomial', alpha=0, standardize=TRUE)
  lambda <- cv.ridge$lambda.1se
  predictions <- predict(cv.ridge, s=lambda, newx=X_predict, type='response')
  return(predictions)
}
''')

"""
zdf_features = list()
for name, group in feature_df.groupby('PCs'):
    group_feature_df = group[feature_names]
    group_feature_zdf = (group_feature_df - group_feature_df.mean()) / group_feature_df.std()
    group_feature_zdf['status'] = group['status']
    group_feature_zdf['target'] = group['target']
    zdf_features.append(group_feature_zdf)

feature_zdf = pandas.concat(zdf_features)
feature_zdf = feature_zdf.fillna(0)
feature_df = feature_zdf
"""

grouped = list(feature_df.groupby('target'))

features = list()
metaedges = set()

y_true = feature_df['status'].values
for feature_name in feature_names:
    print feature_name
    # parse feature
    split_feature = feature_name.split('_')
    metapath = split_feature[0]
    thresholds = [tuple(threshold.split('=')) for threshold in split_feature[1:]]
    feature = {'name': feature_name, 'metapath': metapath}
    features.append(feature)
    for metaedge, value in thresholds:
        feature[metaedge] = value
        metaedges.add(metaedge)
    # calculate AUC
    y_score = feature_df[feature_name]
    feature['auc_global'] = sklearn.metrics.roc_auc_score(y_true, y_score)
    # Logreg to include PCs and PCt for GaD
    X = feature_df[[feature_name, 'PCs', 'PCt']].values
    """
    logreg = sklearn.linear_model.LogisticRegression(C=100, random_state=0)
    logreg.fit(X, y_true)
    y_logreg_score = logreg.decision_function(X)
    feature['auc_logreg'] = sklearn.metrics.roc_auc_score(y_true, y_logreg_score)
    """
    y_glmnet_score = rpy2.robjects.r['RidgePredictions'](X, y_true, X)
    feature['auc_logreg'] = sklearn.metrics.roc_auc_score(y_true, y_glmnet_score)

    """# Ridge Classifier with CV to include PCs and PCt for GaD
    alphas = [0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.5, 1.0]
    ridge_cv = sklearn.linear_model.RidgeClassifierCV(alphas, normalize=True)
    ridge_cv.fit(X, y_true)
    print 'CV Alpha', ridge_cv.alpha_
    y_ridge_score = ridge_cv.decision_function(X)
    feature['auc_ridge'] = sklearn.metrics.roc_auc_score(y_true, y_ridge_score)
    # Group by disease
    group_aucs = list()
    group_positives = list()
    for name, group in grouped:
        y_true = group['status']
        y_score = group[feature_name]
        group_aucs.append(sklearn.metrics.roc_auc_score(y_true, y_score))
        group_positives.append(sum(y_true))
    feature['auc_grouped'] = numpy.average(group_aucs)
    feature['auc_grouped_weighted'] = numpy.average(group_aucs, weights=group_positives)
    """

fieldnames = ['name', 'metapath', 'auc_global', 'auc_grouped', 'auc_grouped_weighted', 'auc_logreg'] + list(metaedges)
auc_path = os.path.join(network_dir, 'feature-aucs.txt')
auc_file = open(auc_path, 'w')
writer = csv.DictWriter(auc_file, delimiter='\t', fieldnames=fieldnames, extrasaction='ignore')
writer.writeheader()
writer.writerows(features)
auc_file.close()
