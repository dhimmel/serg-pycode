import argparse
import os
import gzip
import csv
import random
import itertools
import operator

import bioparser.data
import hetnet.readwrite.graph
import hetnet.permutation
import hetnet.display

part_fieldnames = ['disease_code', 'disease_name', 'gene_code',
    'gene_symbol', 'status', 'status_int', 'percentile', 'part']


def get_part_rows(graph, percent_training=0.75, min_genes=10, seed=0):

    # Calculate gene-disease pairs
    symbol_to_gene = bioparser.data.Data().hgnc.get_symbol_to_gene()
    metaedge_DaG = graph.metagraph.edge_dict[('disease', 'gene', 'association', 'both')]
    metanode_to_nodes = graph.get_metanode_to_nodes()
    genes = sorted(metanode_to_nodes[metaedge_DaG.target])
    diseases = sorted(disease for disease in metanode_to_nodes[metaedge_DaG.source])
    diseases = [disease for disease in diseases if len(disease.edges[metaedge_DaG]) >= min_genes]
    status_to_rows = dict()
    all_rows = list()
    for disease, gene in itertools.product(diseases, genes):
        edge = graph.edge_dict.get((gene.id_, disease.id_, 'association', 'both'))
        #network_status = int(any(edge.target == disease for edge in gene.edges[metaedge_GaD]))
        status = 'HC_primary' if edge else 'negative'
        status_int = int(bool(edge))
        gene_symbol = gene.id_
        gene_code = symbol_to_gene[gene_symbol].hgnc_id
        disease_code = disease.id_
        disease_name = disease.data['name']
        row = {'disease_code': disease_code, 'disease_name': disease_name,
               'gene_code': gene_code, 'gene_symbol': gene_symbol,
               'status': status, 'status_int': status_int}
        status_to_rows.setdefault(status, list()).append(row)
        all_rows.append(row)

    # Assign associations as either testing or training
    random.seed(seed)
    for status, rows in status_to_rows.iteritems():
        n = len(rows)
        rindexes = range(n)
        random.shuffle(rindexes)
        for rindex, row in zip(rindexes, rows):
            percentile = round(float(rindex + 1) / n, 9)
            row['percentile'] = percentile
            row['part'] = 'train' if percentile <= percent_training else 'test'

    all_rows.sort(key=operator.itemgetter('disease_name', 'gene_symbol'))
    return all_rows, diseases


def write_partition_file(all_rows, directory):
    partition_path = os.path.join(directory, 'partitions.txt.gz')
    partition_file = gzip.open(partition_path, 'w')
    writer = csv.DictWriter(partition_file, delimiter='\t', fieldnames=part_fieldnames)
    writer.writeheader()
    writer.writerows(all_rows)
    partition_file.close()


def write_partition_files(all_rows, directory, diseases):
    partition_dir = os.path.join(directory, 'disease-partitions')
    if not os.path.isdir(partition_dir):
        os.mkdir(partition_dir)

    doid_to_writer = dict()
    write_files = list()
    for disease in diseases:
        disease_code = disease.id_
        path = os.path.join(partition_dir, '{}.txt.gz'.format(disease_code.replace(':', '_')))
        write_file = gzip.open(path, 'w')
        write_files.append(write_file)
        writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=part_fieldnames)
        writer.writeheader()
        doid_to_writer[disease_code] = writer

    for row in all_rows:
        doid_to_writer[row['disease_code']].writerow(row)

    for write_file in write_files:
        write_file.close()


def get_testing_edges(graph, all_rows):
    # create set of testing edges to exclude while saving network
    testing_edges = set()
    for row in all_rows:
        if row['status'] != 'HC_primary':
            continue
        if row['part'] != 'test':
            continue
        edge = graph.edge_dict[(row['gene_symbol'], row['disease_code'], 'association', 'both')]
        testing_edges.add(edge)
        testing_edges.add(edge.inverse)
    return testing_edges


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--origin-dir', type=os.path.expanduser,
        default='~/Documents/serg/gene-disease-hetnet/networks/XXXXXX-all-assoc')
    parser.add_argument('--networks-dir', type=os.path.expanduser,
        default='~/Documents/serg/gene-disease-hetnet/networks/')
    parser.add_argument('--percent-training', type=float, default=0.75)
    parser.add_argument('--min-genes', type=int, default=10)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--network-number', type=int, default=3)
    parser.add_argument('--multiplier', type=int, default=10)
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--prefix', default='140615')

    args = parser.parse_args()

    print 'Reading unpermuted graph'
    origin_pkl_path = os.path.join(args.origin_dir, 'graph', 'graph.pkl.gz')
    graph = hetnet.readwrite.graph.read_pickle(origin_pkl_path)

    current_graph = graph

    # Permute graph
    for i in range(args.network_number):
        print 'Permuting network {}'.format(i).center(60, '#')
        network_id = '{}-permuted{}'.format(args.prefix, i)
        network_ids = [network_id]

        training = i == 0
        if training:
            training_id = '{}-permuted{}-training'.format(args.prefix, i)
            network_ids.append(training_id)

        pkl_paths = list()
        network_dirs = list()
        for netid in network_ids:
            network_dir = os.path.join(args.networks_dir, netid)
            network_dirs.append(network_dir)
            graph_dir = os.path.join(network_dir, 'graph')
            for directory in (network_dir, graph_dir):
                if not os.path.isdir(directory):
                    os.mkdir(directory)
            pkl_path = os.path.join(network_dir, 'graph.pkl.gz')
            if not args.overwrite:
                assert not os.path.exists(pkl_path)
            pkl_paths.append(pkl_path)

        current_graph = hetnet.permutation.permute_graph(
            current_graph, multiplier=args.multiplier, seed=args.seed, verbose=True)
        print 'Permuted graph metrics'.center(60, '-')
        print hetnet.display.graph_metrics(current_graph)

        # compute and save partitions
        all_rows, diseases = get_part_rows(
            current_graph, percent_training=args.percent_training, min_genes=args.min_genes)
        write_partition_file(all_rows, network_dirs[0])
        write_partition_files(all_rows, network_dirs[0], diseases=diseases)

        # save permuted graph
        print 'saving permuted graph'
        hetnet.readwrite.graph.write_pickle(current_graph, pkl_paths[0])
        # Save training-testing graph
        if training:
            testing_edges = get_testing_edges(current_graph, all_rows)
            print '{} testing edges'.format(len(testing_edges))
            hetnet.readwrite.graph.write_pickle(
                current_graph, pkl_paths[0], exclude_edges=testing_edges)
            write_partition_file(all_rows, network_dirs[1])
            print 'saving permuted training graph'
            write_partition_files(all_rows, network_dirs[1], diseases=diseases)

