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
                       ('disease', 'tissue', 'pathology', 'both'),
                       ('disease', 'factor', 'involvement', 'both')]
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
    doid_graph = data.doid.get_graph()
    for doid_id, nx_data in doid_graph.nodes(data=True):
        node_data = {'name': nx_data['name']}
        graph.add_node(doid_id, 'disease', node_data)
    
    # Print metanode counter
    for metanode, nodes in graph.get_metanode_to_nodes().items():
        print metanode, len(nodes)

    # Add (disease, gene, association, both) edges
    logging.info('Adding GWAS catalog disease-gene associations.')
    doid_id_to_genes = data.gwas_catalog.get_doid_id_to_genes(
        fdr_cutoff=0.05, mapped_term_cutoff=1)
    for doid_id, genes in doid_id_to_genes.iteritems():
        for gene in genes:
            symbol = gene.symbol
            graph.add_edge(doid_id, symbol, 'association', 'both')

    # Add (gene, tissue, expression, both) edges
    logging.info('Adding BodyMap2 gene-tissue expression.')
    fpkm_cutoff = 50.0
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
        
    
    # Print metaedge counter
    for metanode, edges in graph.get_metaedge_to_edges().items():
        print metanode, len(edges)
        
    """

    # Add Factor Nodes
    for factor in data.etiome.get_factors():
        factor_id = 'factor: ' + factor
        graph.add_node(factor_id, kind='factor')

    # Add (disease, factor, involvement) edges
    disease_to_factors = data.etiome.get_disease_to_factors()
    etiome_to_efo_id = mappings.get_mappings('etiome_name', 'efo_id')
    for etiome_disease, efo_id in etiome_to_efo_id.items():
        factors = disease_to_factors[etiome_disease]
        for factor in factors:
            factor_id = 'factor: ' + factor
            graph.add_edge(efo_id, factor_id, 'involvement', 'both')

    # Add (disease, gene, association) edges
    section = "('disease', 'gene', 'association', 'both')"
    exclude_pmids = set(ast.literal_eval(config.get(section, 'exclude_pmids')))# {'21833088'}
    efo_id_to_genes = data.gwas_catalog.get_efo_id_to_genes(
        fdr_cutoff=config.getfloat(section, 'fdr_cutoff'),
        mapped_term_cutoff=config.getint(section, 'mapped_term_cutoff'),
        exclude_pmids=exclude_pmids)
    for efo_id, gcat_symbols in efo_id_to_genes.iteritems():
        if efo_id not in graph.node_dict:
            continue
        matched_symbols = gcat_symbols & hgnc_symbols
        for gcat_symbol in matched_symbols:
            gene_symbol = symbol_to_gene[gcat_symbol].symbol
            graph.add_edge(efo_id, gene_symbol, 'association', 'both')

    # Add (gene, tissue, specificity) edges
    gene_to_tissues = data.tiger.get_gene_to_tissues()
    for symbol, tissues in gene_to_tissues.iteritems():
        gene = symbol_to_gene.get(symbol)
        if gene is None:
            continue
        hgnc_symbol = gene.symbol
        for tissue in tissues:
            graph.add_edge(hgnc_symbol, tissue, 'specificity', 'both')

    # Add (disease, tissue, pathology) edges
    section = "('disease', 'tissue', 'pathology', 'both')"
    r_scaled_cutoff = config.getfloat(section, 'r_scaled_cutoff')
    coocc_gen = copub_analysis.tiger_efo_cooccurrence_generator()
    for row in coocc_gen:
        disease = row['efo_disease_id']
        tissue = row['tiger_tissue']
        r_scaled = row['r_scaled']
        if r_scaled < r_scaled_cutoff:
            continue
        edge_data = {'r_scaled': r_scaled}
        graph.add_edge(disease, tissue, 'pathology', 'both', edge_data)
    
    
    # Add (disease, gene, regulation) edges
    # Add (disease, gene, up-regulation) edges
    # Add (disease, gene, down-regulation) edges
    gxa_reader = bioparser.gxa.Reader()
    for efo_id in disease_terms:
        genes_tuple = gxa_reader.get_genes(efo_id, p_cutoff=0.05)
        if genes_tuple is None:
            continue
        down_symbols, up_symbols = genes_tuple
        down_symbols = {symbol_to_gene[gxa_symbol].symbol
                        for gxa_symbol in down_symbols & hgnc_symbols}
        up_symbols = {symbol_to_gene[gxa_symbol].symbol
                      for gxa_symbol in up_symbols & hgnc_symbols}
        for symbol in down_symbols | up_symbols:
            assert efo_id in g and g.node[efo_id]['kind'] == 'disease'
            assert symbol in g and g.node[symbol]['kind'] == 'gene'
            g.add_edge(efo_id, symbol, key='regulation')
        for symbol in down_symbols:
            g.add_edge(efo_id, symbol, key='down-regulation')
        for symbol in up_symbols:
            g.add_edge(efo_id, symbol, key='up-regulation')

    # (gene, gene, interaction) information:
    for interaction in data.iref.get_interactions(min_publications=2):
        symbol_a, symbol_b = interaction
        graph.add_edge(symbol_a, symbol_b, 'interaction', 'both')

    # (gene, gene, function) information
    section = "('gene', 'gene', 'function', 'both')"
    processed_file = config.get(section, 'processed_file')
    include_positives = config.getboolean(section, 'include_positives')
    include_predictions = config.getboolean(section, 'include_predictions')
    probability_cutoff = config.getfloat(section, 'probability_cutoff')

    if processed_file:
        relationships = data.frimp.read_processed_relationships(processed_file)
    else:
        relationships = data.frimp.get_relationships(
            predictions = include_predictions,
            positives = include_positives,
            prob_cutoff = probability_cutoff)
    for symbol_a, symbol_b, prob in relationships:
        edge_data = {'probability': prob}
        graph.add_edge(symbol_a, symbol_b, 'function', 'both', edge_data)
    """

    return graph




if __name__ == '__main__':
    # Parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--network-dir', type=os.path.expanduser, default=
        '~/Documents/serg/ashg13/131007-1')
    parser.add_argument('--config', action='store_true')
    parser.add_argument('--create', action='store_true')
    args = parser.parse_args()
    network_dir = args.network_dir
    graph_agent = hetnet.agents.GraphAgent(network_dir)
    graph_dir = graph_agent.graph_dir
    
    
    if args.config:
        write_config(graph_dir)
    
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
