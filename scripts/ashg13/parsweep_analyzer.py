import os
import gzip
import pprint
import csv

import sklearn
import sklearn.metrics
import numpy
import pandas
import pandas.io.parsers

network_dir = '/home/dhimmels/Documents/serg/ashg13/140302-parsweep'


feature_path = os.path.join(network_dir, 'features-subset-fine.txt.gz')
feature_df = pandas.io.parsers.read_table(feature_path, compression='gzip')

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
    #metaedges |= {threshold[0] for threshold in thresholds}
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
auc_path = os.path.join(network_dir, 'feature-aucs.txt')
auc_file = open(auc_path, 'w')
writer = csv.DictWriter(auc_file, delimiter='\t', fieldnames=fieldnames, extrasaction='ignore')
writer.writeheader()
writer.writerows(features)
auc_file.close()
