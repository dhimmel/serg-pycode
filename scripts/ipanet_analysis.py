import os
import sys

import networkx

import networks.networkx_extensions as nxext


ipanet_dir = '/home/dhimmels/Documents/serg/ipanet/'
pkl_path = os.path.join(ipanet_dir, 'ipanet.pkl')
pkl_with_shortcuts_path = os.path.join(ipanet_dir, 'ipanet-with-shortcuts.pkl')
if os.path.exists(pkl_with_shortcuts_path):
    pkl_path = os.path.join(ipanet_dir, 'ipanet-with-shortcuts.pkl')
g = networkx.read_gpickle(pkl_path)
print 'loaded network from pickle'

# prepare schema and metapaths
schema = g.graph['schema']
source_kind = 'drug'
target_kind = 'disease'
metapaths = schema.metapaths(source_kind, target_kind, 4)
#for metapath in metapaths: print schema.path_as_str(metapath)
#nxext.print_edge_kind_counts(g)


# compute shorcuts
shortcuts = nxext.shortcuts_for_metapaths(metapaths)
#nxext.compute_shortcuts(g, shortcuts)
print shortcuts
#networkx.write_gpickle(g, pkl_with_shortcuts_path)


source = 'interferon beta-1a'
target = 'multiple sclerosis'

print nxext.path_counter(g, metapaths, source, target, shortcuts=shortcuts)
print nxext.path_counter(g, metapaths, source, target=None, shortcuts=shortcuts)


#print g.edges(source, keys=True)


for metapath in metapaths:
    break
    print schema.path_as_str(metapath)
    paths = nxext.get_paths(g, metapath, source, target)
    print paths
    print len(paths)
    print '----------------------------------------------------------'