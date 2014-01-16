

import hetnet
import hetnet.agents
import bioparser.data

import reader

#msigdb = bioparser.data.Data().msigdb
#msig_set_types = msigdb.abbrev_to_name.keys()

metaedge_tuples = [('drug', 'protein', 'target', 'both')]#,
                   #('gene', 'gene', 'interaction', 'both'),
                   #('gene', 'gene', 'function', 'both'),
                   #('gene', 'tissue', 'expression', 'both')]

#metaedge_tuples.extend([('gene', set_type, 'membership', 'both') for set_type in msig_set_types])
metagraph = hetnet.MetaGraph.from_edge_tuples(metaedge_tuples)
graph = hetnet.Graph(metagraph)

# Add screened drugs
screen_reader = reader.ScreenReader()
screened_compounds = screen_reader.get_screened_compounds()
for compound in screened_compounds:
    name = compound['name']
    smiles = compound['canonical_smiles']
    graph.add_node(smiles, 'drug', compound)


targets = screen_reader.get_SEA_targets()
for target in targets:
    identifier, name, description = target
    node_data = {'name': name, 'description': description}
    graph.add_node(identifier, 'protein', node_data)

sea_p_cutoff = 1.0
interactions = screen_reader.get_SEA_interactions()
for interaction in interactions:
    smiles = interaction['Query Smiles']
    target = interaction['Target ID']
    p_val = interaction['p-value']
    if p_val > sea_p_cutoff:
        continue
    edge_data = {'p_value': p_val, 'max_tc': interaction['Max Tc'],
                 'affinity': interaction['Affinity Threshold (nM)']}
    graph.add_edge(smiles, target, 'target', 'both', edge_data)


"""
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

"""


# Print metanode counter
for metanode, nodes in graph.get_metanode_to_nodes().items():
    line = '{}: {}'.format(metanode, len(nodes))
    print line

# Print metaedge counter
for metaedge, edges in graph.get_metaedge_to_edges(exclude_inverts=True).items():
    line = '{}: {}'.format(metaedge, len(edges))
    print line

network_dir = '/home/dhimmels/Documents/serg/remyelination/networks/140106'
graph_agent = hetnet.agents.GraphAgent(network_dir)
graph_agent.set(graph)
graph_agent.write_additional_formats()