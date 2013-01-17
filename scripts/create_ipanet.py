import collections
import os
import random

import networkx
#import sklearn.linear_model

import bioparser.data
import networks.schema


ipanet_dir = '/home/dhimmels/Documents/serg/ipanet/'
#gml_path = os.path.join(ipanet_dir, 'ipanet.gml')

def read_diseases_to_remove():
    path = os.path.join(ipanet_dir, 'manually-removed-disease-nodes.txt')
    with open(path) as f:
        diseases = {line.strip() for line in f}
    return diseases

def disease_function_subset(ipa, printing=False):
    """
    Returns a disease subset of functions. A function is considered a
    disease if its lowercase name is the same as its class and its name is
    not a function category. Build must be run first
    """
    disease_names = set()
    for function in ipa.functions:
        if function.name.lower() == function.function_class.lower():
            disease_names.add(function.name)
    diseases_to_remove = read_diseases_to_remove()
    disease_names -= diseases_to_remove
    
    
    disease_functions = {ipa.name_to_function[disease] for disease in disease_names}
    print len(disease_functions), 'diseases'
    
    ## print random sample of removed function names
    omitted_functions = {function.name for function in ipa.functions - disease_functions}
    if printing:
        for function_name in random.sample(omitted_functions, 20):
            print function_name
            
    return disease_functions


def build_networkx():
    """
    """
    ipa = bioparser.data.Data().ipa
    ipa.build()

    hgnc = bioparser.data.Data().hgnc
    entrez_to_hgnc = hgnc.get_entrez_to_gene()
    symbol_to_hgnc = hgnc.get_symbol_to_gene()
            
    g = networkx.MultiGraph(name='ipanet')
    ################################################################################
    ################################# Create Nodes #################################
    
    # Create drug nodes
    targets = set()
    for drug in ipa.drugs:
        g.add_node(drug.symbol, kind='drug')
        targets |= set(drug.targets)
    
    # Create disease nodes
    disease_functions = disease_function_subset(ipa, printing=False)
    for disease in disease_functions:
        g.add_node(disease.name, kind='disease')
    
    # Create gene nodes
    hugu_genes_added = set()
    for gene in ipa.genes:
        entrez_id = gene.entrez_id_human.split('|')[0]
        hgnc_gene = entrez_to_hgnc.get(entrez_id)
        if hgnc_gene:
            hugu_genes_added.add(hgnc_gene)
        hugu_symbol = hgnc_gene.symbol if hgnc_gene else None
        g.add_node(gene.symbol, hugu_symbol=hugu_symbol, kind='gene')
    for target in targets:
        if target not in g:
            hgnc_gene = symbol_to_hgnc.get(target)
            if hgnc_gene:
                hugu_genes_added.add(hgnc_gene)
            hugu_symbol = hgnc_gene.symbol if hgnc_gene else None
            g.add_node(target, hugu_symbol=hugu_symbol, kind='gene')
    del targets
    
    missing_hugu_genes = set(hgnc.get_genes()) - hugu_genes_added
    for hgnc_gene in missing_hugu_genes:
        if hgnc_gene.symbol in g:
            print hgnc_gene.symbol
            continue
            #raise Exception('pre-existing ipa symbol matching gene name')
        g.add_node(hgnc_gene.symbol, hugu_symbol=hgnc_gene.symbol, kind='gene')
    
    ################################################################################
    ################################# Create Edges #################################
    # Create drug-gene links from drug target annotations.
    for drug in ipa.drugs:
        for target in drug.targets:
            g.add_edge(drug.symbol, target, key='target')
    
    # Create disease-gene and disease-drug links from ipa function annotations.
    for disease in disease_functions:
        for effect, molecules in disease.molecules.items():
            for molecule in molecules:
                if disease.name in g and molecule in g:
                    if g.node[molecule]['kind'] == 'drug':
                        if effect == 'decreases':
                            kind = 'indication'
                        elif effect == 'increases':
                            kind = 'increase'
                        else:
                            kind = 'effect'
                    else:
                        kind = 'risk'
                    g.add_edge(disease.name, molecule, effect=effect, key=kind)
    
    
    # Define schema for network
    node_kinds = {'drug', 'disease', 'gene'}
    edge_tuples = [('drug', 'gene', 'target'),
                   ('gene', 'disease', 'risk'),
                   ('drug', 'disease', 'indication')]
    schema = networks.schema.UndirectedSchema(node_kinds, edge_tuples)
    g.graph['schema'] = schema
    g.graph['prepared'] = False
    return g


def purify(g):
    """Keep only edges of specified kinds and then remove unconnected nodes."""
    # Remove improper edge kinds
    valid_edge_kinds = {'target', 'risk', 'indication', 'increase', 'effect'}
    for node, neighbor, key in g.edges_iter(keys=True):
        if key not in valid_edge_kinds:
            g.remove_edge(node, neighbor, key)
    
    # Remove unconnected nodes
    unconnected_nodes = (node for node, degree in g.degree_iter() if not degree)
    g.remove_nodes_from(unconnected_nodes)        

def node_degree_list(g, printing=False):
    degree_list = networkx.degree(g).items()
    degree_list.sort(key=lambda x: x[1], reverse=True)
    if printing:
        for i in range(100):
            print degree_list[i]
    return degree_list

def save_as_pickle(g, network_id):
    pkl_path = os.path.join(ipanet_dir, 'networks', network_id, 'graph.pkl')
    networkx.write_gpickle(g, pkl_path)
    print 'IPA network saved to', pkl_path
    

if __name__ == '__main__':
    g = build_networkx()
    purify(g)
    print networkx.info(g)
    save_as_pickle(g, '130116-1')


