import os
import gzip
import pprint
import csv

import sklearn
import sklearn.metrics
import numpy
import pandas
import pandas.io.parsers

#network_dir = '/home/dhimmels/Documents/serg/ashg13/140310-parsweep'
project_dir = '/home/dhimmels/Documents/serg/gene-disease-hetnet'
network_dir = os.path.join(project_dir, 'networks', '140313-thresholdsweep')



part_path = os.path.join(project_dir, 'partitions.txt.gz')
part_file = gzip.open(part_path)
part_reader = csv.DictReader(part_file, delimiter='\t')
#part_df = pandas.io.parsers.read_table(part_path, compression='gzip')
training_tuples = {(row['doid_code'], row['gene'])
                   for row in part_reader
                   if row['part'] == 'train'}
part_file.close()

feature_path = os.path.join(network_dir, 'features-exp0.4.txt.gz')
feature_df = pandas.io.parsers.read_table(feature_path, compression='gzip')
is_training = feature_df.apply(lambda x: (x['target'], x['source']) in training_tuples, axis=1)
training_feature_df = feature_df.loc[is_training]
feature_df = training_feature_df

column_names = feature_df.columns.values.tolist()
feature_names = column_names[4:]


grouped = list(feature_df.groupby('target'))

features = list()
metaedges = set()
for feature_name in feature_names:
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
    y_true = feature_df['status']
    y_score = feature_df[feature_name]
    feature['auc_global'] = sklearn.metrics.roc_auc_score(y_true, y_score)
    group_aucs = list()
    group_positives = list()
    for name, group in grouped:
        y_true = group['status']
        y_score = group[feature_name]
        group_aucs.append(sklearn.metrics.roc_auc_score(y_true, y_score))
        group_positives.append(sum(y_true))
    feature['auc_grouped'] = numpy.average(group_aucs)
    feature['auc_grouped_weighted'] = numpy.average(group_aucs, weights=group_positives)

fieldnames = ['name', 'metapath', 'auc_global', 'auc_grouped', 'auc_grouped_weighted'] + list(metaedges)
auc_path = os.path.join(network_dir, 'feature-aucs-training.txt')
auc_file = open(auc_path, 'w')
writer = csv.DictWriter(auc_file, delimiter='\t', fieldnames=fieldnames, extrasaction='ignore')
writer.writeheader()
writer.writerows(features)
auc_file.close()
