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



def create_graph(network_id):
    data = bioparser.data.Data()
    
    symbol_to_gene = data.hgnc.get_symbol_to_gene()
    hgnc_symbols = set(symbol_to_gene)
    
    
    # Define and initialize networkx graph
    edge_metapaths = [('disease', 'gene', 'up-regulation'),
                      ('disease', 'gene', 'down-regulation'),
                      ('disease', 'gene', 'regulation'),
                      ('disease', 'gene', 'association'),
                      ('disease', 'disease', 'comorbidity'),
                      ('gene', 'gene', 'interaction')]
    g = heteronets.nxutils.create_undirected_network(edge_metapaths)
    g.graph['name'] = 'ashg-net'
    g.graph['network_id'] = network_id
    g.graph['description'] = 'Network designed for predicted disease associated genes.'
    #heteronets.schema.print_schema(g.graph['schema'])
    
    # Add Gene Nodes
    for gene in data.hgnc.get_genes():
        g.add_node(gene.symbol, name=gene.name, kind='gene')
    
    # Add Disease Nodes
    disease_root = 'EFO_0000408'
    efo_id_to_name = data.efo.get_id_to_name()
    efo_graph = data.efo.get_graph()
    disease_terms = list(networkx.dfs_postorder_nodes(efo_graph, source=disease_root))
    for disease_term in disease_terms:
        name = efo_id_to_name[disease_term]
        g.add_node(disease_term, name=name, kind='disease')
    
    # Add (disease, gene, association) edges
    gcat = data.gwas_catalog
    efo_id_to_genes = gcat.get_efo_id_to_genes(fdr_cutoff=0.05, mapped_term_cutoff=1)
    for efo_id, gcat_symbols in efo_id_to_genes.iteritems():
        if efo_id not in g:
            continue
        matched_symbols = gcat_symbols & hgnc_symbols
        for gcat_symbol in matched_symbols:
            gene_symbol = symbol_to_gene[gcat_symbol].symbol
            assert efo_id in g and g.node[efo_id]['kind'] == 'disease'
            assert gene_symbol in g and g.node[gene_symbol]['kind'] == 'gene'
            g.add_edge(efo_id, gene_symbol, key='association')
    
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
        assert symbol_a in g and g.node[symbol_a]['kind'] == 'gene'
        assert symbol_b in g and g.node[symbol_b]['kind'] == 'gene'
        g.add_edge(symbol_a, symbol_b, key='interaction')
    
    heteronets.nxutils.remove_unconnected_nodes(g)
    """
    kind_to_node = heteronets.nxutils.get_kind_to_nodes(g)
    for node in kind_to_node['disease']:
        degree_counter = heteronets.nxutils.node_degree_counter(g, node)
        #print degree_counter
        if ('disease', 'gene', 'regulation') in degree_counter and ('disease', 'gene', 'association') in degree_counter:
            print '---------------------------'
            print g.node[node]['name']
            for edge_kind, degree in degree_counter.items():
                print edge_kind, degree
        else:
            g.remove_nodes_from([node])
    
    for node in kind_to_node['gene']:
        degree_counter = heteronets.nxutils.node_degree_counter(g, node)
    #    if ('gene', 'disease', 'regulation') not in degree_counter and ('gene', 'disease', 'association') not in degree_counter:
        if ('gene', 'disease', 'association') not in degree_counter:
           g.remove_nodes_from([node])
    """
    print 'After filtering unconnected nodes'
    heteronets.nxutils.print_node_kind_counts(g)
    heteronets.nxutils.print_edge_kind_counts(g)
    return g

def degree_printer(g):
    degree_tuples = sorted(g.degree_iter(), key=lambda x: x[1], reverse=True)
    for node, degree in degree_tuples[:10]:
        degree_counter = heteronets.nxutils.node_degree_counter(g, node)
        print '------------------'
        print g.node[node]['name']
        print degree_counter
    

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

def select_positives(g):
    kind_to_node = heteronets.nxutils.get_kind_to_nodes(g)

    # Find all diseases with expression and association edges
    diseases = kind_to_node['disease']
    learning_diseases = set()
    for disease in diseases:
        degree_counter = heteronets.nxutils.node_degree_counter(g, disease)
        if ('disease', 'gene', 'association') not in degree_counter:
            continue
        if ('disease', 'gene', 'regulation') not in degree_counter:
            continue
        learning_diseases.add(disease)
    
    # Take all associations edges from those nodes as positives
    positives = g.edges(learning_diseases, keys=True)
    positives = filter(lambda e: e[2] == 'association', positives)
    positives = sorted((e[1], e[0], e[2]) for e in positives)
    return positives


if __name__ == '__main__':
    args = get_parser_args()

    # Create network_dir if necessary
    network_dir = os.path.join(args.project_dir, args.network_id)
    if not os.path.isdir(network_dir):
        os.mkdir(network_dir)

    # Create raw graph if necesarry
    pkl_path = os.path.join(network_dir, 'raw-graph.pkl')
    if not os.path.exists(pkl_path):
        g = create_graph(args.network_id)
        degree_printer(g)
        gml_path = os.path.join(network_dir, 'raw-graph.gml')
        heteronets.nxutils.export_as_gml(g, gml_path)
        heteronets.nxutils.write_gpickle(g, pkl_path)
    else:
        g = heteronets.nxutils.read_gpickle(pkl_path)


    prepared_pkl_path = os.path.join(network_dir, 'prepared-graph.pkl')    
    if not os.path.exists(prepared_pkl_path):
        # Select positives
        g.graph['source_kind'] = 'gene'
        g.graph['target_kind'] = 'disease'
        g.graph['edge_key'] = 'association'
        positives = select_positives(g)
        print len(positives), 'positives'
        
        positive_to_exclusions, negative_to_exclusions = heteronets.metapaths.matched_negatives(positives)
        g.graph['positive_to_exclusions'] = positive_to_exclusions
        g.graph['negative_to_exclusions'] = negative_to_exclusions
        g.graph['negatives'] = set(negative_to_exclusions)
        g.graph['positives'] = set(positive_to_exclusions)
        
        
        schema = g.graph['schema']
        g.graph['metapaths'] = heteronets.schema.extract_metapaths(
            g.graph['schema'], g.graph['source_kind'], g.graph['target_kind'],
            max_length = 2, exclude_all_source_target_edges = False)
    
        print 'Writing prepared graph to pkl.'
        heteronets.nxutils.write_gpickle(g, prepared_pkl_path)
    else:
        g = heteronets.nxutils.read_gpickle(prepared_pkl_path)
    

    # Compute and learning features
    if False:
        positive_to_exclusions = g.graph['positive_to_exclusions']
        negative_to_exclusions = g.graph['negative_to_exclusions']
        edge_to_exclusions = dict(negative_to_exclusions.items() + positive_to_exclusions.items())
        path = os.path.join(network_dir, 'learning-features.txt')
        heteronets.features.write_features(g, edge_to_exclusions, path)


    # Compute and save multiple sclerosis features
    edge_to_exclusions = dict()
    genes = heteronets.nxutils.get_kind_to_nodes(g)['gene']
    pairs = list(itertools.product(genes, ['EFO_0003885']))
    for pair in pairs:
        edge = (pair[0], pair[1], 'association')
        exclusions = {edge, (edge[1], edge[0], edge[2])}
        edge_to_exclusions[edge] = exclusions    
    path = os.path.join(network_dir, 'features-multiple-sclerosis.txt')
    heteronets.features.write_features(g, edge_to_exclusions, path)

