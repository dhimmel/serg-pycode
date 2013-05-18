import csv
import os

import networkx

import bioparser.gxa
import bioparser.data
import heteronets.nxutils
import heteronets.schema


project_dir = '/home/dhimmels/Documents/serg/ashg13'

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
g.graph['network_id'] = '130517-1'
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
for interaction in data.iref.all_interactions():
    symbol_a, symbol_b = interaction
    assert symbol_a in g and g.node[symbol_a]['kind'] == 'gene'
    assert symbol_b in g and g.node[symbol_b]['kind'] == 'gene'
    g.add_edge(symbol_a, symbol_b, key='interaction')

heteronets.nxutils.remove_unconnected_nodes(g)
print 'After filtering unconnected nodes'
heteronets.nxutils.print_node_kind_counts(g)
heteronets.nxutils.print_edge_kind_counts(g)
kind_to_node = heteronets.nxutils.get_kind_to_nodes(g)
for node in kind_to_node['disease']:
    degree_counter = heteronets.nxutils.node_degree_counter(g, node)
    #print degree_counter
    if ('disease', 'gene', 'regulation') in degree_counter and ('disease', 'gene', 'association') in degree_counter:
        print '---------------------------'
        print g.node[node]['name']
        for edge_kind, degree in degree_counter.items():
            print edge_kind, degree

## SAVE networkx.
"""
union_terms = gxa_terms & set(efo_id_to_genes)

print len(union_terms)

symbol_to_gene = data.hgnc.get_symbol_to_gene()
hgnc_symbols = set(symbol_to_gene)

diseases_path = os.path.join(project_dir, 'diseases.txt')
with open(diseases_path, 'w') as diseases_file:
    writer = csv.writer(diseases_file, delimiter='\t')
    writer.writerow(['efo_id', 'name', 'Associations', 'Regulations'])
    for efo_id in union_terms:
        name = efo_id_to_name[efo_id]
        associations = efo_id_to_genes[efo_id] & hgnc_symbols
        regulations = {row['gene'] for row in gxa_reader.read_processed_file(efo_id)} & hgnc_symbols
        writer.writerow([efo_id, name, len(associations), len(regulations)])

def read_efo_doid_mappings():
    efo_to_doid = dict()
    path = os.path.join(project_dir, 'efo-to-doid-mappings-manual.txt')
    mapping_file = open(path)
    reader = csv.DictReader(mapping_file, delimiter='\t')
    for row in reader:
        doid = row['doid']
        efo_to_doid[row['efo']] = doid
    mapping_file.close()
    return efo_to_doid


efo_to_doid = read_efo_doid_mappings()
doids = efo_to_doid.values()
pdn_graph = data.pdn.get_graph().subgraph(doids)
print len(pdn_graph.nodes())
print len(pdn_graph.edges())

"""




