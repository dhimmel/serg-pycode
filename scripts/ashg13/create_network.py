import ast
import argparse
import csv
import os
import collections
import itertools
import logging
import ConfigParser

import bioparser.gxa
import bioparser.data
import hetnet
import hetnet.agents

import copub_analysis
import mappings


def write_config(graph_dir):
    """Writes a default configuration file"""
    config = ConfigParser.SafeConfigParser()
    
    # Set defulat config file
    section = 'disease'
    config.add_section(section)
    config.set(section, 'include', 'True')
    config.set(section, 'disease_root', 'EFO_0000408')

    section = 'gene'
    config.add_section(section)
    config.set(section, 'include', 'True')

    section = 'tissue'
    config.add_section(section)
    config.set(section, 'include', 'True')
    
    section = 'factor'
    config.add_section(section)
    config.set(section, 'include', 'True')
    
    section = "('disease', 'gene', 'association', 'both')"
    config.add_section(section)
    config.set(section, 'include', 'True')
    config.set(section, 'fdr_cutoff', '0.05')
    config.set(section, 'mapped_term_cutoff', '1')
    config.set(section, 'exclude_pmids', '[]')

    section = "('gene', 'gene', 'interaction', 'both')"
    config.add_section(section)
    config.set(section, 'include', 'True')

    section = "('gene', 'gene', 'function', 'both')"
    config.add_section(section)
    config.set(section, 'include', 'True')

    section = "('gene', 'tissue', 'specificity', 'both')"
    config.add_section(section)
    config.set(section, 'include', 'True')

    section = "('disease', 'tissue', 'pathology', 'both')"
    config.add_section(section)
    config.set(section, 'include', 'True')
    config.set(section, 'r_scaled_cutoff', '30.0')

    section = "('disease', 'factor', 'involvement', 'both')"
    config.add_section(section)
    config.set(section, 'include', 'True')

    # Write defualt config file
    config_path = os.path.join(graph_dir, 'graph.cfg')
    assert not os.path.exists(config_path)
    with open(config_path, 'w') as config_file:
        config.write(config_file)
    

def create_graph(config):
    data = bioparser.data.Data()
    
    symbol_to_gene = data.hgnc.get_symbol_to_gene()
    hgnc_symbols = set(symbol_to_gene)
    
    
    # Define and initialize networkx graph
    metaedge_tuples = [('disease', 'gene', 'association', 'both'),
                       ('gene', 'gene', 'interaction', 'both'),
                       ('gene', 'gene', 'function', 'both'),
                       ('gene', 'tissue', 'specificity', 'both'),
                       ('disease', 'tissue', 'pathology', 'both'),
                       ('disease', 'factor', 'involvement', 'both')]
    metagraph = hetnet.MetaGraph.from_edge_tuples(metaedge_tuples)
    graph = hetnet.Graph(metagraph)
    
    # Add Gene Nodes
    logging.info('Adding HGNC genes.')
    for gene in data.hgnc.get_genes():
        node_data = {'name': gene.name}
        graph.add_node(gene.symbol, 'gene', node_data)
    
    # Add TIGER Tissue Nodes
    logging.info('Adding tiger tissue nodes.')
    for tissue in data.tiger.get_tissues():
        graph.add_node(tissue, 'tissue')
    
    # Add Disease Nodes
    section = 'disease'
    disease_root = config.get(section, 'disease_root')
    efo_id_to_name = data.efo.get_id_to_name()
    efo_graph = data.efo.get_graph()
    disease_terms = data.efo.get_non_neoplastic_diseases()
    for disease_term in disease_terms:
        name = efo_id_to_name[disease_term]
        node_data = {'name': name}
        graph.add_node(disease_term, 'disease', node_data)

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
    
    
    """
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
    """

    # (gene, gene, interaction) information:
    for interaction in data.iref.get_interactions(min_publications=2):
        symbol_a, symbol_b = interaction
        graph.add_edge(symbol_a, symbol_b, 'interaction', 'both')

    # (gene, gene, function) information
    relationships = data.frimp.read_processed_relationships('positives_and_predictions_0.5')
    for symbol_a, symbol_b, prob in relationships:
        edge_data = {'probability': prob}
        graph.add_edge(symbol_a, symbol_b, 'function', 'both', edge_data)
    
    #may want to remove unconnected nodes here
    return graph




if __name__ == '__main__':
    # Parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--network-dir', type=os.path.expanduser, default=
        '~/Documents/serg/ashg13/130814-1')
    parser.add_argument('--config', action='store_true')
    parser.add_argument('--create', action='store_true')
    args = parser.parse_args()
    network_dir = args.network_dir
    graph_dir = os.path.join(network_dir, 'graph')
    
    
    if args.config:
        write_config(graph_dir)
    
    if args.create:
        
        # Read configuration file
        config = ConfigParser.SafeConfigParser()
        config_path = os.path.join(graph_dir, 'graph.cfg')
        config.read(config_path)
    
        # Create the graph
        log_path = os.path.join(graph_dir, 'creation.log')
        logging.basicConfig(filename=log_path, level=logging.INFO,
                            filemode='w', format='%(levelname)s:%(message)s')
        graph = create_graph(config)
        
        # Save the graph
        metagraph_agent = hetnet.agents.MetaGraphAgent(network_dir)
        metagraph_agent.write(graph.metagraph)
        
        graph_agent = hetnet.agents.GraphAgent(network_dir)
        graph_agent.write(graph)

