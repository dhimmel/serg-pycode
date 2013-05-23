import argparse
import csv
import os
import collections
import itertools

import networkx

import bioparser.gxa
import bioparser.data
import heteronets.nxutils
import heteronets.schema
import heteronets.metapaths
import heteronets.features

import ashg13_copub



def create_graph(network_id):
    data = bioparser.data.Data()
    
    symbol_to_gene = data.hgnc.get_symbol_to_gene()
    hgnc_symbols = set(symbol_to_gene)
    
    
    # Define and initialize networkx graph
    edge_metapaths = [('disease', 'gene', 'association'),
                      ('gene', 'gene', 'interaction'),
                      ('gene', 'tissue', 'specificity'),
                      ('disease', 'tissue', 'pathology')]
    g = heteronets.nxutils.create_undirected_network(edge_metapaths)
    g.graph['name'] = 'ashg-net'
    g.graph['network_id'] = network_id
    g.graph['description'] = 'Network designed for predicted disease associated genes.'
    #heteronets.schema.print_schema(g.graph['schema'])
    
    # Add Gene Nodes
    for gene in data.hgnc.get_genes():
        g.add_node(gene.symbol, name=gene.name, kind='gene')

    # Add TIGER Tissue Nodes
    for tissue in data.tiger.get_tissues():
        g.add_node(tissue, kind='tissue')
    
    # Add Disease Nodes
    disease_root = 'EFO_0000408'
    efo_id_to_name = data.efo.get_id_to_name()
    efo_graph = data.efo.get_graph()
    disease_terms = list(networkx.dfs_postorder_nodes(efo_graph, source=disease_root))
    for disease_term in disease_terms:
        name = efo_id_to_name[disease_term]
        g.add_node(disease_term, name=name, kind='disease')
    
    # Add (disease, gene, association) edges
    efo_id_to_genes = data.gwas_catalog.get_efo_id_to_genes(fdr_cutoff=0.05, mapped_term_cutoff=1)
    for efo_id, gcat_symbols in efo_id_to_genes.iteritems():
        if efo_id not in g:
            continue
        matched_symbols = gcat_symbols & hgnc_symbols
        for gcat_symbol in matched_symbols:
            gene_symbol = symbol_to_gene[gcat_symbol].symbol
            assert efo_id in g and g.node[efo_id]['kind'] == 'disease'
            assert gene_symbol in g and g.node[gene_symbol]['kind'] == 'gene'
            g.add_edge(efo_id, gene_symbol, key='association')
    
    # Add (gene, tissue, specificity) edges
    gene_to_tissues = data.tiger.get_gene_to_tissues()
    for symbol, tissues in gene_to_tissues.iteritems():
        gene = symbol_to_gene.get(symbol)
        if gene is None:
            continue
        hgnc_symbol = gene.symbol
        for tissue in tissues:
            assert tissue in g
            g.add_edge(hgnc_symbol, tissue, key='specificity')

    # Add (disease, tissue, pathology) edges
    r_scaled_cutoff = 35.0
    coocc_gen = ashg13_copub.tiger_efo_cooccurrence_generator()
    for row in coocc_gen:
        disease = row['efo_disease_id']
        tissue = row['tiger_tissue']
        r_scaled = row['r_scaled']
        if r_scaled < r_scaled_cutoff:
            continue
        assert disease in g
        assert tissue in g
        g.add_edge(disease, tissue, key='pathology')
        
        
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
    for interaction in data.iref.get_interactions(min_publications=1):
        symbol_a, symbol_b = interaction
        assert symbol_a in g and g.node[symbol_a]['kind'] == 'gene'
        assert symbol_b in g and g.node[symbol_b]['kind'] == 'gene'
        g.add_edge(symbol_a, symbol_b, key='interaction')
    
    heteronets.nxutils.remove_unconnected_nodes(g)
    print 'After filtering unconnected nodes'
    heteronets.nxutils.print_node_kind_counts(g)
    heteronets.nxutils.print_edge_kind_counts(g)
    return g



def get_parser_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--networks-dir', type=os.path.expanduser, default=
        '~/Documents/serg/networks/')
    parser.add_argument('--project-dir', type=os.path.expanduser, default=
        '~/Documents/serg/ashg13/')
    parser.add_argument('--network-id', required=True)
    #parser.add_argument('--description', required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_parser_args()

    # Create network_dir if necessary
    network_dir = os.path.join(args.project_dir, args.network_id)
    if not os.path.isdir(network_dir):
        os.mkdir(network_dir)

    g = create_graph(args.network_id)
    pkl_path = os.path.join(network_dir, 'raw-graph.pkl')
    gml_path = os.path.join(network_dir, 'raw-graph.gml')
    heteronets.nxutils.export_as_gml(g, gml_path)
    heteronets.nxutils.write_gpickle(g, pkl_path)
