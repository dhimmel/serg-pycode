import argparse
import collections
import csv
import os
import random
import shutil
import string
import sys

import networkx

import heteronets.features
import heteronets.metapaths
import heteronets.schema
import heteronets.nxutils

################################################################################
################# Execution


parser = argparse.ArgumentParser()
parser.add_argument('--ipadir', type=os.path.expanduser, default=
    '~/Documents/serg/ipanet/')
parser.add_argument('--network-id', required=True)
parser.add_argument('--max-path-length', type=int, required=True)
parser.add_argument('--num-pos', type=int, default=2000)
parser.add_argument('--num-neg', type=int, default=2000)
parser.add_argument('--source-kind', required=True)
parser.add_argument('--target-kind', required=True)
parser.add_argument('--edge-key', required=True)
parser.add_argument('--exclude-edges', type=set, help='set of tuples indicating\
    edges to omit from feature computation.', default=set())
parser.add_argument('--exclude-all-source-target-edges', action='store_true')
parser.add_argument('--reprepare', action='store_true')
args = parser.parse_args()


network_dir = os.path.join(args.ipadir, 'networks', args.network_id)
graph_dir = os.path.join(network_dir, 'graphs')


pkl_path_prepared = os.path.join(graph_dir, 'prepared-graph.pkl')
if not os.path.exists(pkl_path_prepared) or args.reprepare:
    print 'preparing network for feature computation.'
    pkl_path = os.path.join(graph_dir, 'graph.pkl')
    g = networkx.read_gpickle(pkl_path)
    
    g.graph['source_kind'] = args.source_kind
    g.graph['target_kind'] = args.target_kind
    g.graph['edge_key'] = args.edge_key
    g.graph['max_path_length'] = args.max_path_length
    
    metapaths = heteronets.schema.extract_metapaths(
        g.graph['schema'], g.graph['source_kind'],
        g.graph['target_kind'], g.graph['max_path_length'],
        args.exclude_all_source_target_edges, args.exclude_edges)
    g.graph['metapaths'] = metapaths
    
    # Filter nodes
    print 'Before filtering'
    heteronets.nxutils.print_node_kind_counts(g)
    heteronets.nxutils.print_edge_kind_counts(g)
    heteronets.metapaths.total_path_counts(g)
    nodes_to_remove = heteronets.metapaths.nodes_outside_metapaths(g)
    g.remove_nodes_from(nodes_to_remove)
    print 'After filtering'
    heteronets.metapaths.total_path_counts(g)
    heteronets.nxutils.print_node_kind_counts(g)
    heteronets.nxutils.print_edge_kind_counts(g)


    positives, negatives = heteronets.metapaths.learning_edge_subset(
        g, args.num_pos, args.num_neg, seed=0)
    g.graph['positives'] = positives
    g.graph['negatives'] = negatives
    required_source_to_targets = dict()
    for source, target, edge_kind in positives + negatives:
        required_source_to_targets.setdefault(source, set()).add(target)
    g.graph['required_source_to_targets'] = required_source_to_targets
    
    pkl_path_learning = os.path.join(graph_dir, 'learning-graph.pkl')
    networkx.write_gpickle(g, pkl_path_learning)
    heteronets.metapaths.prepare_feature_optimizations(g)
    networkx.write_gpickle(g, pkl_path_prepared)
    print 'prepared graph saved as pickle'
else:
    print 'reading', pkl_path_prepared
    g = networkx.read_gpickle(pkl_path_prepared)

feature_dir = os.path.join(network_dir, 'features')
heteronets.features.write_partitioned_features(g, feature_dir)

"""
# Create a subset with learning 
subset_fxn = lambda feature_dict: feature_dict['status'] in ['0', '1']
subset_feature_file('features', 'features-learning', subset_fxn)

# Create a subset with Multiple Sclerosis 
subset_fxn = lambda feature_dict: feature_dict['target'] == 'Multiple Sclerosis'
subset_feature_file('features', 'features-ms', subset_fxn)
"""
