import ast
import argparse
import csv
import os
import collections
import itertools
import logging
import ConfigParser
import pprint
import gzip

import bioparser.gxa
import bioparser.data
import hetnet
import hetnet.readwrite

def create_graph(associations_path, doidprocess_path, pathophys_path, partition_path):
    data = bioparser.data.Data()
    doid_remove, doid_pop = bioparser.gwas_plus.GwasCatalog.read_ontprocess_info(doidprocess_path)
    exclude_doids = doid_remove | set(doid_pop)

    msigdb = bioparser.data.Data().msigdb
    msig_set_types = ['c1.all', 'c2.cgp', 'c2.cp.biocarta', 'c2.cp.kegg', 'c2.cp.reactome',
                      'c3.mir', 'c3.tft', 'c4.cgn', 'c4.cm', 'c5.bp', 'c5.cc', 'c5.mf',
                      'c6.all', 'c7.all']

    # Define and initialize graph
    metaedge_tuples = [('disease', 'gene', 'association', 'both'),
                       ('gene', 'gene', 'interaction', 'both'),
                       ('disease', 'pathophysiology', 'membership', 'both')]
    metaedge_tuples.extend([('gene', set_type, 'membership', 'both') for set_type in msig_set_types])
    metagraph = hetnet.MetaGraph.from_edge_tuples(metaedge_tuples)
    graph = hetnet.Graph(metagraph)

    # Add genes from HGNC
    logging.info('Adding HGNC gene nodes.')
    for gene in data.hgnc.get_genes():
        if not gene.coding:
            continue
        node_data = {'name': gene.name}
        graph.add_node(gene.symbol, 'gene', node_data)

    # Add diseases from DOID
    logging.info('Adding DOID disease nodes.')
    doid_onto = data.doid.get_ontology()
    for doid_id, nx_data in doid_onto.graph.nodes(data=True):
        if doid_id in exclude_doids:
            continue
        node_data = {'name': nx_data['name']}
        graph.add_node(doid_id, 'disease', node_data)

    # Add pathophysiology nodes
    exclude_pathophys = {'unspecific', 'idiopathic'}
    with open(pathophys_path) as read_file:
        reader = csv.DictReader(read_file, delimiter='\t')
        pathophys_rows = list(reader)
    pathophys_rows = [row for row in pathophys_rows
                      if row['pathophysiology'] not in exclude_pathophys]
    pathophys_rows = [row for row in pathophys_rows
                      if row['disease_code'] not in exclude_doids]
    pathophys_terms = {row['pathophysiology'] for row in pathophys_rows}
    for pathophys_term in pathophys_terms:
        graph.add_node(pathophys_term, 'pathophysiology')

    # Add (disease, pathophysiology, membership, both) edges
    for pathophys_row in pathophys_rows:
        doid_code = pathophys_row['disease_code']
        pathophys_term = pathophys_row['pathophysiology']
        graph.add_edge(doid_code, pathophys_term, 'membership', 'both')

    # Add (disease, gene, association, both) edges
    with gzip.open(partition_path) as part_file:
        reader = csv.DictReader(part_file, delimiter='\t')
        part_rows = [row for row in reader if row['status'] == 'HC_primary']
    assoc_to_part = {(row['disease_code'], row['gene_symbol']): row['part']
                     for row in part_rows}

    logging.info('Adding GWAS catalog disease-gene associations.')
    associations_file = open(associations_path)
    associations_reader = csv.DictReader(associations_file, delimiter='\t')
    doids_with_associations = set()
    for association in associations_reader:
        disease_code = association['disease_code']
        gene_symbol = association['gene_symbol']
        assoc_tuple = disease_code, gene_symbol
        if association['status'] != 'HC_primary':
            continue
        part = assoc_to_part.get(assoc_tuple, 'excluded')
        if part == 'test':
            continue
        graph.add_edge(disease_code, gene_symbol, 'association', 'both')
        doids_with_associations.add(disease_code)
    associations_file.close()


    # Add (gene, gene, interaction, both) edges
    logging.info('Adding ppiTrim gene-gene interaction.')
    interactions = data.ppitrim.collapsed_binary_interactions()
    edge_keys = ['pubmed', 'method_id', 'interaction_type_id']
    for interaction in interactions:
        edge_data = {k: v for k, v in interaction.items() if k in edge_keys}
        source = interaction['source'].symbol
        target = interaction['target'].symbol
        try:
            graph.add_edge(source, target, 'interaction', 'both', edge_data)
        except KeyError:
            pass

    # Add MSigDB gene sets
    for set_type in msig_set_types:
        for name, description, genes in msigdb.gene_set_generator(set_type):
            unique_name = 'MSigDB_{}:{}'.format(set_type, name)
            node_data = {'description': description}
            graph.add_node(unique_name, set_type, node_data)
            for gene in genes:
                try:
                    graph.add_edge(gene.symbol, unique_name, 'membership', 'both')
                except KeyError:
                    pass

    # Print metanode counter
    logging.info('MetaNode Counts')
    for metanode, nodes in graph.get_metanode_to_nodes().items():
        line = '{}: {}'.format(metanode, len(nodes))
        print line
        logging.info(line)

    # Print metaedge counter
    logging.info('MetaEdge Counts')
    for metaedge, edges in graph.get_metaedge_to_edges(exclude_inverts=True).items():
        line = '{}: {}'.format(metaedge, len(edges))
        print line
        logging.info(line)

    return graph




if __name__ == '__main__':
    # Parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--network-dir', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/networks/140614-metricsweep')
    parser.add_argument('--doidprocess-path', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/data-integration/doid-ontprocess-info.txt')
    parser.add_argument('--pathophys-path', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/data-integration/pathophysiology.txt')
    parser.add_argument('--partition-path', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/partitions.txt.gz')
    parser.add_argument('--associations-id', default='processed')
    parser.add_argument('--create', action='store_true')
    args = parser.parse_args()
    network_dir = args.network_dir
    graph_dir = os.path.join(network_dir, 'graph')

    for directory in (network_dir, graph_dir):
        if not os.path.isdir(directory):
            os.mkdir(directory)

    associations_path = os.path.join(
        bioparser.data.Data().gwas_plus.directory,
        args.associations_id, 'association-statuses.txt')

    if args.create:
        
        # Create the graph
        log_path = os.path.join(graph_dir, 'creation.log')
        logging.basicConfig(filename=log_path, level=logging.INFO,
                            filemode='w', format='%(levelname)s:%(message)s')
        graph = create_graph(associations_path=associations_path,
                             doidprocess_path=args.doidprocess_path,
                             pathophys_path=args.pathophys_path,
                             partition_path=args.partition_path)
        
        # Save the graph
        pkl_path = os.path.join(network_dir, 'graph', 'graph.pkl.gz')
        hetnet.readwrite.graph.write_pickle(graph, pkl_path)
