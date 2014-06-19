import argparse
import os
import gzip
import csv

import hetnet.readwrite.graph
import hetnet.permutation


parser = argparse.ArgumentParser()
parser.add_argument('--origin-dir', type=os.path.expanduser,
    default='~/Documents/serg/gene-disease-hetnet/networks/XXXX')
parser.add_argument('--permuted-dir', type=os.path.expanduser,
    default='~/Documents/serg/gene-disease-hetnet/networks/XXXX')
parser.add_argument('--partition-path', type=os.path.expanduser,
    default='~/Documents/serg/gene-disease-hetnet/partitions.txt.gz')
parser.add_argument('--overwrite', action='store_true')
args = parser.parse_args()



# Read unpermuted graph
pkl_path = os.path.join(args.origin_dir, 'graph', 'graph.pkl.gz')
graph = hetnet.readwrite.graph.read_pickle(pkl_path)

metaedge_GaD = graph.metagraph.edge_dict[('gene', 'disease', 'association', 'both')]
metaedge_DaG = metaedge_GaD.inverse

# Read partition information and add testing associations (both negative and positives)
# to a set of excluded permutations
metaedge_to_excluded = {metaedge_GaD: set(), metaedge_DaG: set()}
part_file = gzip.open(args.partition_path)
reader = csv.DictReader(part_file, delimiter='\t')
for row in reader:
    if row['part'] != 'test':
        continue
    disease_code = row['disease_code']
    gene_symbol = row['gene_symbol']
    metaedge_to_excluded[metaedge_GaD].add((gene_symbol, disease_code))
    metaedge_to_excluded[metaedge_DaG].add((disease_code, gene_symbol))
part_file.close()

# Permute graph
permuted_graph = hetnet.permutation.permute_graph(graph, seed=0,
    metaedge_to_excluded=metaedge_to_excluded, verbose=True)

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









