import itertools

import nxutils
import schema
import metapaths

"""
Module to test performance of heteronets modules with smallscale known examples.
"""

edge_metapaths = [('compound', 'gene', 'target'),
                  ('compound', 'disease', 'indication'),
                  ('disease', 'gene', 'risk')]
kind_to_abbrev = {'compound': 'C', 'disease': 'D', 'gene': 'G',
                  'risk': 'r', 'indication':'i', 'target': 't'}
g = nxutils.create_undirected_network(edge_metapaths, kind_to_abbrev)

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

nxutils.print_node_kind_counts(g)
nxutils.print_edge_kind_counts(g)


g.graph['source_kind'] = 'compound'
g.graph['target_kind'] = 'disease'
g.graph['edge_key'] = 'indication'
g.graph['max_path_length'] = 2
    
g.graph['metapaths'] = schema.extract_metapaths(
    g.graph['schema'], g.graph['source_kind'],
    g.graph['target_kind'], g.graph['max_path_length'],
    exclude_all_source_target_edges=True)

metapaths.total_path_counts(g)

kind_to_nodes = nxutils.get_kind_to_nodes(g)

# test g.graph['required_source_to_targets']
if True:
    pairs = list(itertools.product(kind_to_nodes[g.graph['source_kind']],
                                   kind_to_nodes[g.graph['target_kind']]))
    required_source_to_targets = dict()
    for source, target in pairs:
        required_source_to_targets.setdefault(source, set()).add(target)
    g.graph['required_source_to_targets'] = required_source_to_targets


# compute and print normalized path counts
for compound in kind_to_nodes['compound']:
    print compound, metapaths.normalized_path_counter(g, g.graph['metapaths'], compound)
