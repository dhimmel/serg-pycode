import os
import csv
import sys
import random
import collections

import networkx

import networks.networkx_extensions as nxext




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


###########################################
## Compute features for all pairs

schema = g.graph['schema']
metapaths = g.graph['metapaths']
shortcuts = g.graph['shortcuts']


feature_file_path = os.path.join(ipanet_dir, 'networks', network_id, 'features.txt')
feature_file = open(feature_file_path, 'w')
fieldnames = ['source', 'target', 'status'] + [schema.path_as_abbrev_str(metapath) for metapath in metapaths]
dict_writer = csv.DictWriter(feature_file, fieldnames=fieldnames, delimiter='\t')
dict_writer.writeheader()

positives = set(g.graph['positives'])
negatives = set(g.graph['negatives'])

kind_to_nodes = nxext.get_kind_to_nodes(g)
drugs = kind_to_nodes['drug']

metapath_abbrevs = map(schema.path_as_abbrev_str, metapaths)
for source in drugs:
    print source
    target_to_metapath_to_npc = nxext.normalized_path_counter(g, metapaths, source, shortcuts)
    for target, metapath_to_npc in target_to_metapath_to_npc.iteritems():
        metapath_to_npc = {schema.path_as_abbrev_str(key): value for key, value in metapath_to_npc.items()}
        row_dict = dict.fromkeys(metapath_abbrevs, None)
        row_dict.update(metapath_to_npc)
        row_dict['source'] = source
        row_dict['target'] = target
        edge = source, target
        if edge in negatives:
            status = 2
        elif edge in positives:
            status = 3
        else: 
            status = int(g.has_edge(source, target, 'indication'))
        row_dict['status'] = status
        if status < 2:
            continue
        dict_writer.writerow(row_dict)
feature_file.close()
 
