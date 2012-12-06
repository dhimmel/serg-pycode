import os

import networkx

import create_ipanet
import networks.networkx_extensions as nxext

g = networkx.read_gpickle(create_ipanet.pkl_path)
schema = g.graph['schema']

source_kind = 'drug'
target_kind = 'disease'
metapaths = schema.metapaths(source_kind, target_kind, 4)
for metapath in metapaths: print schema.path_as_str(metapath)
#nxext.print_edge_kind_counts(g)

path_heads, path_tails = nxext.get_precompute_metapaths(metapaths)


source = 'interferon beta-1a'
target = 'multiple sclerosis'
#print g.edges(source, keys=True)


for metapath in metapaths:
    break
    print schema.path_as_str(metapath)
    paths = nxext.get_paths(g, source, target, metapath)
    print paths
    print len(paths)
    print '----------------------------------------------------------'