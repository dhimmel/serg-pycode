import os
import csv
import sys
import random
import collections

import networkx

import networks.networkx_extensions as nxext


def edge_status(g, source, target):
    """Returns the status for an edge.
    0 - Negative evaluation edge
    1 - Positive evaluation edge
    2 - Non existent network edge
    3 - Existent network edge
    """
    edge = source, target
    if edge in g.graph['negatives']:
        status = 0
    elif edge in g.graph['positives']:
        status = 1
    else: 
        status = int(g.has_edge(source, target, g.graph['edge_kind'])) + 2

def feature_generator(g):
    """
    Generates features (only NPC for now) based on the graph specifications.
    """
    schema = g.graph['schema']
    metapaths = g.graph['metapaths']
    metapath_abbrevs = [schema.path_as_abbrev_str(metapath) for metapath in metapaths]
    kind_to_nodes = nxext.get_kind_to_nodes(g)
    sources = kind_to_nodes[g.graph['source_kind']]
    for source in sources:
        print source
        target_to_metapath_to_npc = nxext.normalized_path_counter(g, metapaths, source)
        for target, metapath_to_npc in target_to_metapath_to_npc.iteritems():
            metapath_to_npc = {schema.path_as_abbrev_str(key): value for key, value in metapath_to_npc.items()}
            feature_dict = dict.fromkeys(metapath_abbrevs, None)
            feature_dict.update(metapath_to_npc)
            feature_dict['source'] = source
            feature_dict['target'] = target
            status = edge_status(g, source, target)
            yield feature_dict    


def write_features(g, path, status_subset = None):
    """
    If status_subset is None all edges are written despite status.
    """
    schema = g.graph['schema']
    metapaths = g.graph['metapaths']
    feature_file = open(feature_file_path, 'w')
    fieldnames = ['source', 'target', 'status'] + [schema.path_as_abbrev_str(metapath) for metapath in metapaths]
    dict_writer = csv.DictWriter(feature_file, fieldnames=fieldnames, delimiter='\t')
    dict_writer.writeheader()

    for feature_dict in feature_generator(g):
        if status_subset is None or feature_dict['status'] in status_subset:
            dict_writer.writerow(feature_dict)
    feature_file.close()
        

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

feature_file_path = os.path.join(ipanet_dir, 'networks', network_id, 'features.txt')
write_features(g, feature_file_path)


