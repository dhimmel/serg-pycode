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
import hetnet.agents
from projects.gene_disease_hetnet.data_integration import copub_analysis



def create_graph(associations_path, doidprocess_path, partition_path):
    data = bioparser.data.Data()
    doid_remove, doid_pop = bioparser.gwas_plus.GwasCatalog.read_ontprocess_info(doidprocess_path)
    exclude_doids = doid_remove | set(doid_pop)


    # Define and initialize networkx graph
    metaedge_tuples = [('disease', 'gene', 'association', 'both'),
                       ('gene', 'tissue', 'expression', 'both'),
                       ('disease', 'tissue', 'cooccurrence', 'both')]
    metagraph = hetnet.MetaGraph.from_edge_tuples(metaedge_tuples)
    graph = hetnet.Graph(metagraph)

    # Add genes from HGNC
    logging.info('Adding HGNC gene nodes.')
    for gene in data.hgnc.get_genes():
        if not gene.coding:
            continue
        node_data = {'name': gene.name}
        graph.add_node(gene.symbol, 'gene', node_data)


    # Add tissues from BTO
    logging.info('Adding BTO tissue nodes.')
    bto_graph = data.bto.get_animal_tissue_subgraph()
    for bto_id, nx_data in bto_graph.nodes(data=True):
        node_data = {'name': nx_data['name']}
        graph.add_node(bto_id, 'tissue', node_data)


    # Add diseases from DOID
    logging.info('Adding DOID disease nodes.')
    doid_onto = data.doid.get_ontology()
    for doid_id, nx_data in doid_onto.graph.nodes(data=True):
        if doid_id in exclude_doids:
            continue
        node_data = {'name': nx_data['name']}
        graph.add_node(doid_id, 'disease', node_data)


    # Add (disease, gene, association, both) edges
    with gzip.open(partition_path) as part_file:
        reader = csv.DictReader(part_file, delimiter='\t')
        part_rows = [row for row in reader if row['status'] == 'assoc_high']
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
        if association['status'] != 'assoc_high':
            continue
        part = assoc_to_part.get(assoc_tuple, 'excluded')
        if part == 'test':
            continue
        graph.add_edge(disease_code, gene_symbol, 'association', 'both')
        doids_with_associations.add(disease_code)
    associations_file.close()


    # Add (gene, tissue, expression, both) edges
    logging.info('Adding GNF gene-tissue expression.')
    log10_expr_cutoff = 0.5
    logging.info('log10_expression_cutoff: {}'.format(log10_expr_cutoff))
    expressions = data.gnf.expression_generator()
    for expression in expressions:
        if expression['gene'].locus_group != 'protein-coding gene':
            continue
        if expression['log10_expr'] < log10_expr_cutoff:
            continue
        edge_data = {'log10_expr': expression['log10_expr']}
        graph.add_edge(expression['bto_id'], expression['gene'].symbol, 'expression', 'both', edge_data)


    # Add (disease, tissue, cooccurrence, both)
    logging.info('Adding CoPub disease-tissue cooccurrence.')
    r_scaled_cutoff = 10
    logging.info('r_scaled_cutoff: {}'.format(r_scaled_cutoff))
    coocc_gen = copub_analysis.doid_bto_cooccurrence_generator()
    for row in coocc_gen:
        doid_id = row['doid_id']
        if doid_id not in doids_with_associations:
            continue
        bto_id = row['bto_id']
        r_scaled = row['r_scaled']
        if r_scaled < r_scaled_cutoff:
            continue
        edge_data = {'r_scaled': r_scaled}
        graph.add_edge(doid_id, bto_id, 'cooccurrence', 'both', edge_data)

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
        '~/Documents/serg/gene-disease-hetnet/networks/140518-thresholdsweep')
    parser.add_argument('--doidprocess-path', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/data-integration/doid-ontprocess-info.txt')
    parser.add_argument('--partition-path', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/partitions.txt.gz')
    parser.add_argument('--associations-id', default='processed')
    parser.add_argument('--create', action='store_true')
    args = parser.parse_args()
    network_dir = args.network_dir
    graph_agent = hetnet.agents.GraphAgent(network_dir)
    graph_dir = graph_agent.graph_dir
    
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
                             partition_path=args.partition_path)
        
        # Save the graph
        graph_agent = hetnet.agents.GraphAgent(network_dir)
        graph_agent.set(graph)
        graph_agent.write_additional_formats()
