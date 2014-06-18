import argparse
import os

import hetnet.readwrite.graph
import hetnet.permutation


parser = argparse.ArgumentParser()
parser.add_argument('--origin-dir', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/networks/XXXX')
parser.add_argument('--permuted-dir', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/networks/XXXX')
parser.add_argument('--overwrite', action='store_true')
args = parser.parse_args()


# Read unpermuted graph
pkl_path = os.path.join(args.origin_dir, 'graph', 'graph.pkl.gz')
graph = hetnet.readwrite.graph.read_pickle(pkl_path)

# Permute graph
permuted_graph = hetnet.permutation.permute_graph(graph, verbose=True)

# Create directories for permuted network and save graph
permuted_dir = args.permuted_dir
pgraph_dir = os.path.join(permuted_dir, 'graph')
ppkl_path = os.path.join(pgraph_dir, 'graph.pkl.gz')
for directory in (permuted_dir, pgraph_dir):
    if not os.path.isdir(directory):
        os.mkdir(directory)

if not args.overwrite:
    assert not os.path.exists(ppkl_path)

hetnet.readwrite.graph.write_pickle(permuted_graph, ppkl_path)









