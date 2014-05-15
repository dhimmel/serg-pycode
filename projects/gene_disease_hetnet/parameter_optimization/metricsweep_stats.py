import os
import gzip
import pprint
import csv

import sklearn
import sklearn.metrics
import sklearn.linear_model
import sklearn.preprocessing
import numpy
import pandas
import pandas.io.parsers

project_dir = '/home/dhimmels/Documents/serg/gene-disease-hetnet'
network_dir = os.path.join(project_dir, 'networks', '140313-metricsweep')

part_path = os.path.join(project_dir, 'partitions.txt.gz')
with gzip.open(part_path) as part_file:
    part_reader = csv.DictReader(part_file, delimiter='\t')
    training_tuples = {(row['doid_code'], row['gene'])
                       for row in part_reader if row['part'] == 'train'}


feature_path = os.path.join(network_dir, 'features.txt.gz')
feature_df = pandas.io.parsers.read_table(feature_path, compression='gzip')
is_training = feature_df.apply(lambda x: (x['target'], x['source']) in training_tuples, axis=1)
training_feature_df = feature_df.loc[is_training]
feature_df = training_feature_df


column_names = feature_df.columns.values.tolist()
feature_names = column_names[4:]
stat_df = pandas.DataFrame.from_records(
    data=(x.split(':') for x in feature_names),
    columns=['metric', 'metapath'])
stat_df['feature'] = feature_names
metrics = sorted(set(stat_df['metric']))

metric = metrics[0]
X = feature_df[stat_df.loc[stat_df['metric'] == metric]['feature']].values
y_true = feature_df['status'].values

alphas = [0.01, 0.02, 0.05, 0.01, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75,
          1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0]
# Alpha corresponds to (2*C)^-1 in LogisticRegression
#alpha between 0--1, occasionally in 30s
#
#sklearn.linear_model.LogisticRegression(penalty='l2', C=c, random_state=None)

ridge = sklearn.linear_model.RidgeClassifierCV(alphas)
ridge.fit(X, y_true)
ridge.decision_function(X)
"""
column_names = feature_df.columns.values.tolist()
feature_names = column_names[4:]
stat_df = pandas.DataFrame.from_records(
    data=(x.split(':') for x in feature_names),
    columns=['metric', 'metapath'])
metrics = sorted(set(stat_df['metric']))

y_true = feature_df['status'].values

def get_auc(y_score):
    return sklearn.metrics.roc_auc_score(y_true, y_score)

stat_df['auc'] = feature_df[feature_names].fillna(0).apply(get_auc, axis=0).values

stat_df_path = os.path.join(network_dir, 'feature-aucs.txt')
stat_df.to_csv(stat_df_path, sep='\t', index=False)
"""