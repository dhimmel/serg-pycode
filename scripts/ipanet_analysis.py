import os
import csv
import sys
import random
import collections

import networkx

import networks.networkx_extensions as nxext


################################################################################
################################## Functions ###################################

def choose_learning_indications(g, pos_num=6000, neg_num=6000):
    # Select positives and negatives
    kind_to_nodes = nxext.get_kind_to_nodes(g)
    kind_to_edges = nxext.get_kind_to_edges(g)
    drugs = kind_to_nodes['drug']
    diseases = kind_to_nodes['disease']
    ordered_indication = lambda edge: edge if edge[0] in drugs else (edge[1], edge[0])
    indications = list(kind_to_edges['indication'])
    indications = map(ordered_indication, indications)
    positives = random.sample(indications, pos_num)
    negatives = list()
    while len(negatives) < neg_num:
        drug = random.choice(indications)[0]
        disease = random.choice(indications)[1]
        if not g.has_edge(disease, drug):
            edge = drug, disease
            negatives.append(edge)
    # delete positives edges from the network
    g.remove_edges_from(positives)
    return positives, negatives

def prepare(g):
    positives, negatives = choose_learning_indications(g, 2000, 2000)
    g.graph['positives'] = positives
    g.graph['negatives'] = negatives
    # compute shorcuts
    print 'computing shortcuts'
    nxext.compute_shortcuts(g, shortcuts)
    print 'computing total path counts'
    nxext.total_path_counts(g, depth)
    
################################################################################
################################## Execution ###################################


depth = 6
ipanet_dir = '/home/dhimmels/Documents/serg/ipanet/'
pkl_with_shortcuts_path = os.path.join(ipanet_dir, 'ipanet-with-shortcuts.pkl')
if not os.path.exists(pkl_with_shortcuts_path):
    pkl_path = os.path.join(ipanet_dir, 'ipanet.pkl')
    g = networkx.read_gpickle(pkl_path)
    schema = g.graph['schema']
    source_kind, target_kind = 'drug', 'disease'
    metapaths = schema.metapaths(source_kind, target_kind, depth)
    shortcuts = nxext.shortcuts_for_metapaths(metapaths)
    prepare(g)
    networkx.write_gpickle(g, pkl_with_shortcuts_path)
    print 'Shortcut gpickle written'
else:
    g = networkx.read_gpickle(pkl_with_shortcuts_path)
    schema = g.graph['schema']
    source_kind, target_kind = 'drug', 'disease'
    #### MUST BE UNDONE
    metapaths = schema.metapaths(source_kind, target_kind, 5)
    shortcuts = nxext.shortcuts_for_metapaths(metapaths)
    metapaths = schema.metapaths(source_kind, target_kind, depth)

########################
## Compute features
"""
positives = g.graph['positives']
negatives = g.graph['negatives']
feature_file_path = os.path.join(ipanet_dir + 'features.txt')
feature_file = open(feature_file_path, 'w')
fieldnames = ['source', 'target', 'status'] + [schema.path_as_abbrev_str(metapath) for metapath in metapaths]
dict_writer = csv.DictWriter(feature_file, fieldnames=fieldnames, delimiter='\t')
dict_writer.writeheader()
for status in (0, 1):
    edges = positives if status else negatives
    for source, target in edges:
        print source + '\t' + target
        npc_dict = nxext.focussed_normalized_path_counter(g, metapaths, source, target, shortcuts=shortcuts)
        npc_dict = {schema.path_as_abbrev_str(key): value for key, value in npc_dict.items()}
        npc_dict['source'] = source
        npc_dict['target'] = target
        npc_dict['status'] = status
        dict_writer.writerow(npc_dict)
feature_file.close()
"""
#######################
## Compute features for all pairs


feature_file_path = os.path.join(ipanet_dir + 'features-len-6.txt')
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
 

"""
nxext.total_path_counts(g, depth)
print g.node['interferon beta-1a']['all_paths']
print g.node['multiple sclerosis']['all_paths']

#source = 'interferon beta-1a'
#target = 'multiple sclerosis'

#print nxext.path_counter(g, metapaths, source, target, shortcuts=shortcuts)
#print nxext.path_counter(g, metapaths, source, target=None, shortcuts=shortcuts)
"""
