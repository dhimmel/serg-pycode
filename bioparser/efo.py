import obo

import networkx


efo_obo_path = '/home/dhimmels/Documents/serg/data-sources/gxa/efo/trunk/src/efoinobo/efo.obo'
ontology = obo.OBOOntology(efo_obo_path)
g = ontology.to_networkx()
#assert networkx.is_directed_acyclic_graph(g)