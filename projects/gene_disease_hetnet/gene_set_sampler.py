import os
import csv
import gzip
import random

import yard.curve

import hetnet
import hetnet.algorithms

network_dir = os.path.expanduser('~/Documents/serg/gene-disease-hetnet/networks/140615-all-assoc/')
partition_path = os.path.expanduser('~/Documents/serg/gene-disease-hetnet/partitions.txt.gz')
output_path = os.path.join(network_dir, 'model', 'gene-set-subset-aurocs.txt')

damping_exponent = 0.4
repetitions = 10
random.seed(0)

print 'loading graph'
pkl_path = os.path.join(network_dir, 'graph', 'graph.pkl.gz')
graph = hetnet.readwrite.graph.read_pickle(pkl_path)
metagraph = graph.metagraph
print 'graph loaded'


metanode_to_nodes = graph.get_metanode_to_nodes()
metaedge_to_edges = graph.get_metaedge_to_edges()


metaedge_GaD = metagraph.get_edge(('gene', 'disease', 'association', 'both'))
metaedge_DaG = metaedge_GaD.inverse
metaedge_GiG = metagraph.get_edge(('gene', 'gene', 'interaction', 'both'))


gene_set_metapaths = list()
for metapath in metagraph.extract_metapaths('gene', 'gene', max_length=2):
    if len(metapath) != 2:
        continue
    if metaedge_GiG in metapath.edges:
        continue
    if metaedge_GaD in metapath.edges:
        continue
    metapath = metagraph.get_metapath(metapath.edges + (metaedge_GaD, ))
    gene_set_metapaths.append(metapath) 


part_rows = dict()
with gzip.open(partition_path) as part_file:
    part_reader = csv.DictReader(part_file, delimiter='\t')
    for row in part_reader:
        #row['disease_node']
        part_rows.setdefault(row['status'], list()).append(row)

part_rows.keys()


def generate_dgs_tuples(negative_prob = 0.005):
    """ """
    
    positives = part_rows['HC_primary']
    negatives = part_rows['negative']
    
    k = int(round(len(negatives) * negative_prob))
    negatives = random.sample(population = negatives, k = k)
    
    for row in negatives + positives:
        disease_node = graph.node_dict[row['disease_code']]
        gene_node = graph.node_dict[row['gene_symbol']]
        yield disease_node, gene_node, int(row['status_int'])


def get_auroc(dgs_tuples, cache = dict()):
    """Cache is dg_tuple_to_paths"""

    statuses = list()
    dwpcs = list()

    for disease, gene, status in dgs_tuples:

        # exclude edges
        if status:
            edge = graph.edge_dict[(gene.id_, disease.id_, 'association', 'both')]
            exclude_edges = {edge, edge.inverse}
        else:
            exclude_edges = set()

        # caching
        dg_tuple = disease, gene
        if dg_tuple not in cache:
            paths = graph.paths_between_tree(gene, disease, metapath, masked=True, exclude_edges=exclude_edges)
            cache[dg_tuple] = paths 
        paths = cache[dg_tuple]

        # calculate DWPC
        paths = [path for path in paths if not path.is_masked()]
        dwpc = hetnet.algorithms.DWPC(paths, damping_exponent=damping_exponent, exclude_edges=exclude_edges, exclude_masked=True)

        statuses.append(status)
        dwpcs.append(dwpc)

    roc_curve = yard.curve.ROCCurve(data=zip(dwpcs, statuses))
    return roc_curve.auc()

#graph.unmask()

# initiate writer
fieldnames = ['metanode', 'metaedge', 'total_nodes', 'total_edges',
              'mask_type', 'repetition', 'nodes', 'edges', 'auroc']
write_file = open(output_path, 'w')
writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=fieldnames)
writer.writeheader()

for metapath in gene_set_metapaths:

    metanode = metapath.get_nodes()[1]
    metaedge = metapath.edges[1]

    all_nodes = metanode_to_nodes[metanode]
    all_edges = metaedge_to_edges[metaedge]

    n_nodes = len(all_nodes)
    n_edges = len(all_edges)

    node_numbers = [int(round(0.1 * x * n_nodes)) for x in range(1, 11)] + [186]
    edge_numbers = [int(round(0.1 * x * n_edges)) for x in range(1, 11)] + [4456]

    dgs_tuples = list(generate_dgs_tuples())
    cache = dict()

    for node_number in node_numbers:

        for repetition in range(repetitions):

            # randomly calculate nodes to mask
            mask_nodes = random.sample(all_nodes, n_nodes - node_number)

            # mask nodes
            for node in mask_nodes:
                node.mask()

            # calculate number of edges without masked source or target nodes
            edge_number = n_edges - sum(edge.source.is_masked() or edge.target.is_masked() for edge in all_edges)

            # calculate auroc
            auroc = get_auroc(dgs_tuples, cache=cache)

            # unmask nodes
            for node in mask_nodes:
                node.unmask()

            result = {'metanode': metanode, 'metaedge': metaedge,
                      'total_nodes': n_nodes, 'total_edges': n_edges,
                      'mask_type': 'node', 'repetition': repetition,
                      'nodes': node_number, 'edges': edge_number,
                      'auroc': auroc}

            writer.writerow(result)
            print result


    for edge_number in edge_numbers:

        for repetition in range(repetitions):

            # randomly calculate edges to mask
            mask_edges = random.sample(all_edges, n_edges - edge_number)

            # mask edges
            for edge in mask_edges:
                edge.mask()
                edge.inverse.mask()

            # calculate number of nodes with any unmasked edges
            node_number = sum(bool(node.get_edges(metaedge, exclude_masked=True)) for node in all_nodes)

            # calculate auroc
            auroc = get_auroc(dgs_tuples, cache=cache)

            # unmask edges
            for edge in mask_edges:
                edge.unmask()
                edge.inverse.unmask()

            result = {'metanode': metanode, 'metaedge': metaedge,
                      'total_nodes': n_nodes, 'total_edges': n_edges,
                      'mask_type': 'edge', 'repetition': repetition,
                      'nodes': node_number, 'edges': edge_number,
                      'auroc': auroc}

            writer.writerow(result)
            print result

write_file.close()
