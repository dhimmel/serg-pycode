import itertools

import nxutils
import schema
import metapaths
import features

"""
Module to test performance of heteronets modules with smallscale known examples.
"""

edge_metapaths = [('compound', 'gene', 'target'),
                  ('compound', 'disease', 'indication'),
                  ('disease', 'gene', 'risk')]
#kind_to_abbrev = {'compound': 'C', 'disease': 'D', 'gene': 'G',
#                  'risk': 'r', 'indication':'i', 'target': 't'}
g = nxutils.create_undirected_network(edge_metapaths)

nodes = [('clomipramine', 'compound'),
         ('duloxetine', 'compound'),
         ('vanoxerine', 'compound'),
         ('SLC6A2', 'gene'),
         ('SLC6A3', 'gene'),
         ('HLA-DRB1', 'gene'),
         ('depression', 'disease'),
         ('multiple sclerosis', 'disease'),
         ('cocaine dependence', 'disease')]

edges = [('clomipramine', 'depression', 'indication'),
         ('duloxetine', 'depression', 'indication'),
         ('duloxetine', 'multiple sclerosis', 'indication'),
         ('vanoxerine', 'cocaine dependence', 'indication'),
         ('clomipramine', 'SLC6A2', 'target'),
         ('duloxetine', 'SLC6A2', 'target'),
         ('duloxetine', 'SLC6A3', 'target'),
         ('vanoxerine', 'SLC6A3', 'target'),
         ('SLC6A2', 'depression', 'risk'),
         ('SLC6A2', 'multiple sclerosis', 'risk'),
         ('SLC6A3', 'depression', 'risk'),
         ('SLC6A3', 'multiple sclerosis', 'risk'),
         ('SLC6A3', 'cocaine dependence', 'risk'),
         ('HLA-DRB1', 'multiple sclerosis', 'risk')]

for name, kind in nodes:
    g.add_node(name, kind=kind)

for node, neighbor, key in edges:
    g.add_edge(node, neighbor, key=key)

#nxutils.print_node_kind_counts(g)
#nxutils.print_edge_kind_counts(g)


g.graph['source_kind'] = 'compound'
g.graph['target_kind'] = 'disease'
g.graph['edge_key'] = 'indication'
g.graph['max_path_length'] = 2
g.graph['negatives'] = []
g.graph['positives'] = []

g.graph['metapaths'] = schema.extract_metapaths(
    g.graph['schema'], g.graph['source_kind'],
    g.graph['target_kind'], g.graph['max_path_length'],
    exclude_all_source_target_edges=True)
mpaths = g.graph['metapaths']
print mpaths


#metapaths.total_path_counts(g)
#print metapaths.features_for_metapath(g, 'clomipramine', 'multiple sclerosis', 'indication', mpaths[0])


kind_to_nodes = nxutils.get_kind_to_nodes(g)
pairs = list(itertools.product(kind_to_nodes[g.graph['source_kind']],
                                   kind_to_nodes[g.graph['target_kind']]))

for source, target in pairs:
    print source, '---', target
    metapath_to_metric_dict = metapaths.features_for_metapaths(g, source, target, g.graph['edge_key'], mpaths)
    print metapaths.flatten_feature_dict(metapath_to_metric_dict)

# compute and print normalized path counts
#for compound in kind_to_nodes['compound']:
#    print compound, metapaths.normalized_path_counter(g, g.graph['metapaths'], compound)
#"""

# compute and print features
#for feature_dict in features.feature_generator(g):
#    print feature_dict
