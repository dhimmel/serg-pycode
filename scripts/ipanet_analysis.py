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

pkl_path_prepared = os.path.join(ipanet_dir, 'networks', network_id, 'prepared-graph.pkl')
if not os.path.exists(pkl_path_prepared):
    print 'preparing network for feature computation.'
    pkl_path = os.path.join(ipanet_dir, 'networks', network_id, 'graph.pkl')
    g = networkx.read_gpickle(pkl_path)
    
    ###########################################################################
    ### Variables to set
    max_path_length = 6
    edge_kind_tuple = ('drug', 'indication', 'disease')
    num_pos, num_neg = 2000, 2000

    nxext.prepare_for_feature_computation(g, max_path_length, edge_kind_tuple, num_pos, num_neg)
    print 'writing prepared gpickle'
    networkx.write_gpickle(g, pkl_path_prepared)
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
