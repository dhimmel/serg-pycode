import os
import csv
import sys
import random

import networkx

import networks.networkx_extensions as nxext


ipanet_dir = '/home/dhimmels/Documents/serg/ipanet/'
pkl_path = os.path.join(ipanet_dir, 'ipanet.pkl')
"""
pkl_with_shortcuts_path = os.path.join(ipanet_dir, 'ipanet-with-shortcuts.pkl')
if os.path.exists(pkl_with_shortcuts_path):
    pkl_path = os.path.join(ipanet_dir, 'ipanet-with-shortcuts.pkl')
"""
g = networkx.read_gpickle(pkl_path)
print 'loaded network from pickle'

# prepare schema and metapaths
schema = g.graph['schema']
source_kind = 'drug'
target_kind = 'disease'
metapaths = schema.metapaths(source_kind, target_kind, 4)
#for metapath in metapaths: print schema.path_as_str(metapath)
#nxext.print_edge_kind_counts(g)





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

positives, negatives = choose_learning_indications(g, 2000, 2000)

# compute shorcuts
shortcuts = nxext.shortcuts_for_metapaths(metapaths)
print 'computing shortcuts'
nxext.compute_shortcuts(g, shortcuts)
#print shortcuts
#networkx.write_gpickle(g, pkl_with_shortcuts_path)


feature_file_path = os.path.join(ipanet_dir + 'features.txt')
feature_file = open(feature_file_path, 'w')
fieldnames = ['status'] + [schema.path_as_str(metapath) for metapath in metapaths]
dict_writer = csv.DictWriter(feature_file, fieldnames=fieldnames, delimiter='\t')
dict_writer.writeheader()
for status in (0, 1):
    edges = positives if status else negatives
    for source, target in edges:
        print source + '\t' + target
        npc_dict = nxext.normalized_path_counter(g, metapaths, source, target, shortcuts=shortcuts)
        npc_dict['status'] = status
        dict_writer.writerow(npc_dict)
feature_file.close()



#source = 'interferon beta-1a'
#target = 'multiple sclerosis'

#print nxext.path_counter(g, metapaths, source, target, shortcuts=shortcuts)
#print nxext.path_counter(g, metapaths, source, target=None, shortcuts=shortcuts)

