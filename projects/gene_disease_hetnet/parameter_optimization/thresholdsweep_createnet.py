import ast
import argparse
import csv
import os
import collections
import itertools
import logging
import ConfigParser
import pprint

import bioparser.gxa
import bioparser.data
import hetnet
import hetnet.agents
from projects.gene_disease_hetnet.data_integration import copub_analysis



def create_graph():
    data = bioparser.data.Data()

    msigdb = bioparser.data.Data().msigdb
    #msig_set_types = msigdb.abbrev_to_name.keys()
    # http://www.broadinstitute.org/gsea/msigdb/collections.jsp
    #msig_set_types = ['c1.all', 'c2.cgp', 'c2.cp.all', 'c3.mir', 'c3.tft',
    #                  'c4.cgn', 'c4.cm', 'c5.bp', 'c5.cc', 'c5.mf', 'c6.all', 'c7.all']
    msig_set_types = list()

    # Define and initialize networkx graph
    metaedge_tuples = [('disease', 'gene', 'association', 'both'),
                       ('gene', 'gene', 'interaction', 'both'),
                       ('gene', 'gene', 'function', 'both'),
                       ('gene', 'tissue', 'expression', 'both'),
                       ('disease', 'tissue', 'cooccurrence', 'both'),
                       ('disease', 'disease', 'similarity', 'both')]
    metaedge_tuples.extend([('gene', set_type, 'membership', 'both') for set_type in msig_set_types])
    metagraph = hetnet.MetaGraph.from_edge_tuples(metaedge_tuples)
    graph = hetnet.Graph(metagraph)

    # Add genes from HGNC
    logging.info('Adding HGNC gene nodes.')
    for gene in data.hgnc.get_genes():
        if gene.locus_group != 'protein-coding gene':
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
        node_data = {'name': nx_data['name']}
        graph.add_node(doid_id, 'disease', node_data)

    # Add (disease, gene, association, both) edges
    exclude_doids = {'DOID:0050589', 'DOID:2914'} # IBD and immune system disease
    logging.info('Adding GWAS catalog disease-gene associations.')
    associations_path = os.path.join(data.gwas_plus.directory, 'associations.txt')
    associations_file = open(associations_path)
    associations_reader = csv.DictReader(associations_file, delimiter='\t')
    doids_with_associations = set()
    for association in associations_reader:
        doid_code = association['doid_code']
        if doid_code in exclude_doids:
            continue
        graph.add_edge(doid_code, association['symbol'], 'association', 'both')
        doids_with_associations.add(doid_code)
    associations_file.close()

    # Add (disease, disease, similarity, both)
    logging.info('Adding Intrinsic Semantic Similarity disease-disease edges.')
    lin_cutoff = 0.0
    logging.info('lin_cutoff: {}'.format(lin_cutoff))
    similarity_generator = doid_onto.pairwise_similarities(doids_with_associations)
    similarities = list(similarity_generator)
    for similarity in similarities:
        lin_similarity = similarity['lin_similarity']
        if lin_similarity < lin_cutoff:
            continue
        source = similarity['source']
        target = similarity['target']
        edge_data = {'lin_similarity': lin_similarity}
        graph.add_edge(source, target, 'similarity', 'both', edge_data)


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

    # (gene, gene, function, both) information
    logging.info('Adding IMP gene-gene relationships.')
    processed_file = None #'positives_and_predictions_0.5'
    include_positives = True
    include_predictions = True
    probability_cutoff = 0.2

    if processed_file:
        logging.info('processed_file: {}'.format(processed_file))
        relationships = data.frimp.read_processed_relationships(processed_file)
    else:
        logging.info('include_positives: {}'.format(include_positives))
        logging.info('include_predictions: {}'.format(include_predictions))
        logging.info('probability_cutoff: {}'.format(probability_cutoff))
        relationships = data.frimp.get_relationships(
            predictions = include_predictions,
            positives = include_positives,
            prob_cutoff = probability_cutoff)
    for symbol_a, symbol_b, prob in relationships:
        edge_data = {'probability': prob}
        try:
            graph.add_edge(symbol_a, symbol_b, 'function', 'both', edge_data)
        except KeyError:
            pass

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
        '~/Documents/serg/gene-disease-hetnet/networks/140313-thresholdsweep')
    parser.add_argument('--config', action='store_true')
    parser.add_argument('--create', action='store_true')
    args = parser.parse_args()
    network_dir = args.network_dir
    graph_agent = hetnet.agents.GraphAgent(network_dir)
    graph_dir = graph_agent.graph_dir
    
        
    if args.create:
        
        # Create the graph
        log_path = os.path.join(graph_dir, 'creation.log')
        logging.basicConfig(filename=log_path, level=logging.INFO,
                            filemode='w', format='%(levelname)s:%(message)s')
        graph = create_graph()
        
        # Save the graph
        graph_agent = hetnet.agents.GraphAgent(network_dir)
        graph_agent.set(graph)
        graph_agent.write_additional_formats()
