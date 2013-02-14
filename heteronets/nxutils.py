import networkx

"""
Graph Attributes

schema
kind_to_nodes
kind_to_edges
abbrev_to_path
"""

def get_kind_to_nodes(g, refresh=False):
    """Create a dictionary of node kind to nodes. The dictionary is returned
    and stored as a graph attribute named 'kind_to_nodes'. refresh specifies
    whether the dictionary should be recomputed if it already exists."""
    kind_to_nodes = g.graph.get('kind_to_nodes')
    if not refresh and kind_to_nodes:
        return kind_to_nodes
    kind_to_nodes = dict()
    for node, data in g.nodes_iter(data=True):
        kind = data['kind']
        kind_to_nodes.setdefault(kind, set()).add(node)
    g.graph['kind_to_nodes'] = kind_to_nodes
    return kind_to_nodes

def get_kind_to_edges(g, refresh=False):
    """Create a dictionary of edge kind to edges. The dictionary is returned
    and stored as a graph attribute named 'kind_to_edges'. refresh specifies
    whether the dictionary should be recomputed if it already exists."""
    kind_to_edges = g.graph.get('kind_to_edges')
    if not refresh and kind_to_edges:
        return kind_to_edges
    kind_to_edges = dict()
    for node, neighbor, key in g.edges_iter(keys=True):
        edge = node, neighbor
        kind_to_edges.setdefault(key, set()).add(edge)
    g.graph['kind_to_edges'] = kind_to_edges 
    return kind_to_edges

def print_node_kind_counts(g):
    kind_to_nodes = get_kind_to_nodes(g)
    for key, value in kind_to_nodes.items():
        print key, len(value)

def print_edge_kind_counts(g):
    kind_to_edges = get_kind_to_edges(g)
    for key, value in kind_to_edges.items():
        print key, len(value)

def create_undirected_network(edge_metapaths, kind_to_abbrev, **kwargs):
    """ """
    g = networkx.MultiGraph()
    g.graph.update(kwargs)
    g.graph['schema'] = create_undirected_schema(edge_metapaths, kind_to_abbrev)
