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

import copub_analysis
import mappings



def create_graph():
    data = bioparser.data.Data()
    
    symbol_to_gene = data.hgnc.get_symbol_to_gene()
    hgnc_symbols = set(symbol_to_gene)
    
    
    # Define and initialize networkx graph
    metaedge_tuples = [('disease', 'gene', 'association', 'both'),
                       ('gene', 'gene', 'interaction', 'both'),
                       ('gene', 'gene', 'function', 'both'),
                       ('gene', 'tissue', 'expression', 'both'),
                       ('disease', 'tissue', 'cooccurrence', 'both'),
                       ('disease', 'disease', 'similarity', 'both')]
    metagraph = hetnet.MetaGraph.from_edge_tuples(metaedge_tuples)
    graph = hetnet.Graph(metagraph)
    
    # Add genes from HGNC
    logging.info('Adding HGNC gene nodes.')
    for gene in data.hgnc.get_genes():
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
    
    # Print metanode counter
    logging.info('MetaNode Counts')
    for metanode, nodes in graph.get_metanode_to_nodes().items():
        line = '{}: {}'.format(metanode, len(nodes))
        print line
        logging.info(line)

    # Add (disease, gene, association, both) edges
    logging.info('Adding GWAS catalog disease-gene associations.')
    fdr_cutoff = 0.05
    mapped_term_cutoff = 1
    logging.info('fdr_cutoff: {}'.format(fdr_cutoff))
    logging.info('mapped_term_cutoff: {}'.format(mapped_term_cutoff))
    doid_id_to_genes = data.gwas_catalog.get_doid_id_to_genes(
        fdr_cutoff=fdr_cutoff, mapped_term_cutoff=mapped_term_cutoff)
    for doid_id, genes in doid_id_to_genes.iteritems():
        for gene in genes:
            symbol = gene.symbol
            graph.add_edge(doid_id, symbol, 'association', 'both')
    doids_with_associations = [doid_id for doid_id, genes in
                               doid_id_to_genes.iteritems() if len(genes) > 0]

    # Add (disease, disease, similarity, both)
    logging.info('Adding Instrinsic Semantic Similarity disease-disease edges.')
    lin_cutoff = 0.45
    logging.info('lin_cutoff: {}'.format(lin_cutoff))
    similarity_generator = doid_onto.pairwise_similarities(doids_with_associations)
    for similarity in similarity_generator:
        lin_similarity = similarity['lin_similarity']
        if lin_similarity < lin_cutoff:
            continue
        source = similarity['source']
        target = similarity['target']
        edge_data = {'lin_similarity': lin_similarity}
        graph.add_edge(source, target, 'similarity', 'both', edge_data)
    
    # Add (gene, tissue, expression, both) edges
    logging.info('Adding BodyMap2 gene-tissue expression.')
    fpkm_cutoff = 75.0
    logging.info('fpkm_cutoff: {}'.format(fpkm_cutoff))
    edge_tuples = data.bodymap2.get_edges(fpkm_cutoff)
    for symbol, tissue, fpkm in edge_tuples:
        edge_data = {'fpkm': fpkm}
        graph.add_edge(tissue, symbol, 'expression', 'both', edge_data)
    
    # Add (gene, gene, interaction, both) edges
    logging.info('Adding ppiTrim gene-gene interaction.')
    interactions = data.ppitrim.collapsed_binary_interactions()
    edge_keys = ['pubmed', 'method_id', 'interaction_type_id']
    for interaction in interactions:
        edge_data = {k: v for k, v in interaction.items() if k in edge_keys}
        source = interaction['source'].symbol
        target = interaction['target'].symbol
        graph.add_edge(source, target, 'interaction', 'both', edge_data)
        
    # (gene, gene, function) information
    logging.info('Adding IMP gene-gene relationships.')
    processed_file = None #'positives_and_predictions_0.5'
    include_positives = True
    include_predictions = True
    probability_cutoff = 0.5

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
        graph.add_edge(symbol_a, symbol_b, 'function', 'both', edge_data)

    # Add (disease, tissue, cooccurrence, both)
    logging.info('Adding CoPub disease-tissue cooccurrence.')
    r_scaled_cutoff = 30.0
    logging.info('r_scaled_cutoff: {}'.format(r_scaled_cutoff))
    coocc_gen = copub_analysis.doid_bto_cooccurrence_generator()
    for row in coocc_gen:
        doid_id = row['doid_id']
        bto_id = row['bto_id']
        r_scaled = row['r_scaled']
        if r_scaled < r_scaled_cutoff:
            continue
        edge_data = {'r_scaled': r_scaled}
        graph.add_edge(doid_id, bto_id, 'pathology', 'both', edge_data)
    

    # Print metaedge counter
    logging.info('MetaEdge Counts')
    for metaedge, edges in graph.get_metaedge_to_edges().items():
        line = '{}: {}'.format(metaedge, len(edges))
        print line
        logging.info(line)
    

    return graph




if __name__ == '__main__':
    # Parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--network-dir', type=os.path.expanduser, default=
        '~/Documents/serg/ashg13/131010-1')
    parser.add_argument('--config', action='store_true')
    parser.add_argument('--create', action='store_true')
    args = parser.parse_args()
    network_dir = args.network_dir
    graph_agent = hetnet.agents.GraphAgent(network_dir)
    graph_dir = graph_agent.graph_dir
    
        
    if args.create:
        
        """
        # Read configuration file
        config = ConfigParser.SafeConfigParser()
        config_path = os.path.join(graph_dir, 'graph.cfg')
        config.read(config_path)
        """
        
        # Create the graph
        log_path = os.path.join(graph_dir, 'creation.log')
        logging.basicConfig(filename=log_path, level=logging.INFO,
                            filemode='w', format='%(levelname)s:%(message)s')
        graph = create_graph()
        
        # Save the graph
        graph_agent = hetnet.agents.GraphAgent(network_dir)
        graph_agent.set(graph)
        graph_agent.write_additional_formats()
