import collections
import os
import random

import networkx
#import sklearn.linear_model

import bioparser.data
import networks.schema


ipanet_dir = '/home/dhimmels/Documents/serg/ipanet/'
pkl_path = os.path.join(ipanet_dir, 'ipanet.pkl')
gml_path = os.path.join(ipanet_dir, 'ipanet.gml')

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
    for disease in ipa.functions:
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
            raise Exception('pre-existing ipa symbol matching gene name')
        g.add_node(hgnc_gene.symbol, hugu_symbol=hgnc_gene.symbol, kind='gene')
    
    ################################################################################
    ################################# Create Edges #################################
    # Create drug-gene links from drug target annotations.
    for drug in ipa.drugs:
        for target in drug.targets:
            g.add_edge(drug.symbol, target, key='target')
    
    # Create disease-gene and disease-drug links from ipa function annotations.
    for disease in ipa.functions:
        for effect, molecules in disease.molecules.items():
            for molecule in molecules:
                if disease.name in g and molecule in g:
                    if g.node[molecule]['kind'] == 'drug':
                        if effect == 'decreases':
                            kind = 'indication'
                        else:
                            kind = 'disease_modifying_drug'
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
    
    return g


def purify(g):
    """Keep only edges of specified kinds and then remove unconnected nodes."""
    # Remove improper edge kinds
    valid_edge_kinds = {'target', 'risk', 'indication'}
    for node, neighbor, key in g.edges_iter(keys=True):
        if key not in valid_edge_kinds:
            g.remove_edge(node, neighbor, key)
    
    # Remove unconnected nodes
    unconnected_nodes = (node for node, degree in g.degree_iter() if not degree)
    g.remove_nodes_from(unconnected_nodes)        

if __name__ == '__main__':
    g = build_networkx()
    print networkx.info(g)
    purify(g)
    print networkx.info(g)
    networkx.write_gml(g, gml_path)
    print 'IPA network written as GML'
    networkx.write_gpickle(g, pkl_path)
    print 'IPA network written as pickle'




"""
pkl_path = '/home/dhimmels/Documents/serg/ipanet/ipanet.pkl'
g = networkx.read_gpickle(pkl_path)


# Select positives and negatives
indications = list(kind_to_edges['indication'])
num_of_positives = len(indications) / 200
positives = random.sample(indications, num_of_positives)
drugs = kind_to_nodes['drug']
diseases = kind_to_nodes['disease']
negatives = list()
while len(negatives) < num_of_positives:
    disease = (set(random.choice(indications)) & diseases).pop()
    drug = (set(random.choice(indications)) & drugs).pop()
    if not g.has_edge(disease, drug):
        edge = drug, disease
        negatives.append(edge)

# delete positives edges from the network
g.remove_edges_from(positives)
total_path_counts(g)

# Create predictor and response arrays
training_edges = negatives + positives
y = numpy.repeat([0, 1], [len(negatives), len(positives)])

npcs_by_edge = [normalized_path_counts(*edge) for edge in training_edges]

metapaths = set()
for npc in npcs_by_edge:
    metapaths |= set(npc.keys())

metapaths = list(metapaths)
metapaths.sort(key=lambda x: len(x))


X = list()
for npc in npcs_by_edge:
    x_row = list()
    for metapath in metapaths:
        value = npc.get(metapath)
        if value is None:
            value = 0.0
        x_row.append(value)
    X.append(x_row)

X = numpy.array(X)

logreg = sklearn.linear_model.LogisticRegression()
logreg.fit(X, y)
y_predicted = logreg.predict_proba(X)[:,1]

fpr, tpr, thresholds = sklearn.metrics.roc_curve(y, y_predicted)
sklearn.metrics.auc(fpr, tpr)

feature_file = '/home/dhimmels/Documents/serg/ipanet/features.txt'
numpy.savetxt(feature_file, numpy.column_stack((X, y.T)))


#print g.node['multiple sclerosis']




        




source = 'multiple sclerosis'
target = 'interferon beta-1a'
normalized_path_counts(source, target)



path_counts = path_counts(source, target)

#print compute_metapath_counter(source, source)
#print compute_metapath_counter(target, target)




#print networkx.info(g)
#node_kind_counts = collections.Counter(data['kind'] for node, data in g.nodes_iter(data=True))
#print node_kind_counts


#print 'Number of connected components:', networkx.number_connected_components(g)
"""