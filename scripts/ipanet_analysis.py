import collections
import csv
import os
import random
import shutil
import string
import sys

import networkx

import networks.networkx_extensions as nxext
import networks.features

################################################################################
################# Execution

ipanet_dir = '/home/dhimmels/Documents/serg/ipanet/'
network_id = '130116-1'
network_dir = os.path.join(ipanet_dir, 'networks', network_id)
graph_dir = os.path.join(network_dir, 'graphs')

pkl_path_prepared = os.path.join(graph_dir, 'prepared-graph.pkl')
if not os.path.exists(pkl_path_prepared):
    print 'preparing network for feature computation.'
    pkl_path = os.path.join(graph_dir, 'graph.pkl')
    g = networkx.read_gpickle(pkl_path)
    
    ### Variables to set
    max_path_length = 6
    edge_kind_tuple = ('drug', 'indication', 'disease')
    num_pos, num_neg = 2000, 2000

    nxext.prepare_feature_computation(g, max_path_length, edge_kind_tuple, num_pos, num_neg)
    pkl_path_learning = os.path.join(graph_dir, 'learning-graph.pkl')
    networkx.write_gpickle(g, pkl_path_learning)
    nxext.prepare_feature_optimizations(g)
    networkx.write_gpickle(g, pkl_path_prepared)
    print 'prepared graph saved as pickle'
else:
    print 'reading', pkl_path_prepared
    g = networkx.read_gpickle(pkl_path_prepared)

feature_dir = os.path.join(ipanet_dir, 'networks', network_id, 'features')
networks.features.write_partitioned_features(g, feature_dir)

"""
# Create a subset with learning 
subset_fxn = lambda feature_dict: feature_dict['status'] in ['0', '1']
subset_feature_file('features', 'features-learning', subset_fxn)

# Create a subset with Multiple Sclerosis 
subset_fxn = lambda feature_dict: feature_dict['target'] == 'Multiple Sclerosis'
subset_feature_file('features', 'features-ms', subset_fxn)
"""
