import math

import hetnet
import hetnet.agents
import bioparser.data

import reader


def ppa_from_pvalue(p_value, prior):
    """
    Calculate Posterior probability of association from a p_value and
    prior probability. doi:10.1038/nrg2615
    """
    bayes_factor = -1.0 / (math.e * p_value * math.log(p_value))
    posterior_odds = bayes_factor * prior / (1.0 - prior)
    ppa = posterior_odds / (1.0 + posterior_odds)
    return ppa


msigdb = bioparser.data.Data().msigdb
#msig_set_types = msigdb.abbrev_to_name.keys()
msig_set_types = ['c2.all', 'c5.all', 'c7.all']

metaedge_tuples = [('drug', 'protein', 'target', 'both'),
                   ('protein', 'gene', 'product', 'both'),
                   ('gene', 'gene', 'interaction', 'both')]
                   #('gene', 'gene', 'function', 'both'),
                   #('gene', 'tissue', 'expression', 'both')]

metaedge_tuples.extend([('gene', set_type, 'membership', 'both') for set_type in msig_set_types])
metagraph = hetnet.MetaGraph.from_edge_tuples(metaedge_tuples)
graph = hetnet.Graph(metagraph)

# Add screened drugs
screen_reader = reader.ScreenReader()
screened_compounds = screen_reader.get_screened_compounds()
for compound in screened_compounds:
    name = compound['name']
    smiles = compound['canonical_smiles']
    graph.add_node(smiles, 'drug', compound)

# Add targets
targets = screen_reader.get_SEA_targets()
uniprot_to_protein = dict()
for target in targets:
    identifier, name, description = target
    node_data = {'name': name, 'description': description}
    graph.add_node(identifier, 'protein', node_data)
    if not identifier.startswith('sp_'):
        continue
    uniprot_id = identifier[3:]
    uniprot_to_protein[uniprot_id] = identifier

# Add genes from HGNC
gene_protein_edges = list()
hgnc = bioparser.data.Data().hgnc
for gene in hgnc.get_genes():
    uniprot_id = gene.uniprot_id
    if uniprot_id:
        gene_protein_edges.append((gene.symbol, uniprot_id))
    node_data = {'name': gene.name, 'locus_group': gene.locus_group}
    graph.add_node(gene.symbol, 'gene', node_data)

for symbol, uniprot_id in gene_protein_edges:
    protein_id = uniprot_to_protein.get(uniprot_id)
    if not protein_id:
        continue
    graph.add_edge(symbol, protein_id, 'product', 'both')


sea_p_cutoff = 1.0 / math.e
target_prior = 0.005
interactions = screen_reader.get_SEA_interactions()
for interaction in interactions:
    smiles = interaction['Query Smiles']
    target = interaction['Target ID']
    p_val = interaction['p-value']
    if p_val > sea_p_cutoff:
        continue

    probability = ppa_from_pvalue(p_val, target_prior)
    edge_data = {'p_value': p_val, 'max_tc': interaction['Max Tc'],
                 'affinity': interaction['Affinity Threshold (nM)'],
                 'probability': probability, 'log_p_value': math.log10(p_val)}
    graph.add_edge(smiles, target, 'target', 'both', edge_data)



# Add genes from HGNC
hgnc = bioparser.data.Data().hgnc
for gene in hgnc.get_genes():
    node_data = {'name': gene.name, 'locus_group': gene.locus_group}
    graph.add_node(gene.symbol, 'gene', node_data)

# Add (gene, gene, interaction, both) edges
ppitrim = bioparser.data.Data().ppitrim
interactions = ppitrim.collapsed_binary_interactions()
edge_keys = ['pubmed', 'method_id', 'interaction_type_id']
for interaction in interactions:
    edge_data = {k: v for k, v in interaction.items() if k in edge_keys}
    source = interaction['source'].symbol
    target = interaction['target'].symbol
    graph.add_edge(source, target, 'interaction', 'both', edge_data)

# Add MSigDB gene sets
for set_type in msig_set_types:
    for name, description, genes in msigdb.gene_set_generator(set_type):
        unique_name = 'MSigDB_{}:{}'.format(set_type, name)
        node_data = {'description': description}
        graph.add_node(unique_name, set_type, node_data)
        for gene in genes:
            graph.add_edge(gene.symbol, unique_name, 'membership', 'both')



# Print metanode counter
for metanode, nodes in graph.get_metanode_to_nodes().items():
    line = '{}: {}'.format(metanode, len(nodes))
    print line

# Print metaedge counter
for metaedge, edges in graph.get_metaedge_to_edges(exclude_inverts=True).items():
    line = '{}: {}'.format(metaedge, len(edges))
    print line

network_dir = '/home/dhimmels/Documents/serg/remyelination/networks/140117'
graph_agent = hetnet.agents.GraphAgent(network_dir)
graph_agent.set(graph)
graph_agent.write_additional_formats()