import obo

import networkx


efo_obo_path = '/home/dhimmels/Documents/serg/data-sources/efo/120508-revision-249.obo.txt'
ontology = obo.OBOOntology(efo_obo_path)
g = ontology.to_directed_networkx()

replace_colon = lambda s: s.replace(':', '_')

# Chemicals
chemical_compounds = list(networkx.dfs_postorder_nodes(g, source='CHEBI:37577')) # CHEBI_37577 chemical compound
query_chemicals = filter(lambda x: g.out_degree(x) == 0, chemical_compounds)
query_chemicals = map(replace_colon, query_chemicals)

# Diseases
diseases = list(networkx.dfs_postorder_nodes(g, source='EFO:0000408')) # EFO_0000408 disease
leaf_diseases = filter(lambda x: g.out_degree(x) == 0, diseases)
query_diseases = set(leaf_diseases)
for leaf_disease in leaf_diseases:
    query_diseases |= set(g.predecessors(leaf_disease))

for node in list(query_diseases):
    num_descendents = len(list(networkx.dfs_postorder_nodes(g, source=node))) - 1
    if num_descendents > 5:
        predecessors.remove(node)

for node in list(query_diseases):
    descendents = set(networkx.dfs_postorder_nodes(g, source=node))
    descendents.remove(node)
    query_diseases -= descendents

query_diseases = map(replace_colon, query_diseases)