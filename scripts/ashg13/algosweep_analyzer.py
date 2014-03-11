import os
import gzip
import pprint
import csv

import sklearn
import sklearn.metrics
import sklearn.linear_model
import numpy
import pandas
import pandas.io.parsers

network_dir = '/home/dhimmels/Documents/serg/ashg13/140305-algosweep'

feature_path = os.path.join(network_dir, 'features-part.txt.gz')
feature_df = pandas.io.parsers.read_table(feature_path, compression='gzip', nrows=20000)

column_names = feature_df.columns.values.tolist()
feature_names = column_names[5:]
metrics = sorted({x.split(':')[0] for x in feature_names})

y = feature_df['status'].values
for metric in metrics:
    pass
metric = metrics[0]
features_for_metric = [f for f in feature_names if f.startswith(metric)]
X = feature_df[features_for_metric].values


numpy.core.defchararray.startswith(feature_names, metrics[0])
feature_df[['PC:G-i-G-a-D', 'PC:G-m-C5.B-m-G-a-D']]

for damping_exponent in [0.1 * x for x in range(11)]:
    pass
    sklearn.linear_model.ElasticNetCV()



feature_array = numpy.array(feature_df[5:])

column_names[]
feature_dfvalues