import networkx

import schema

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
        edge = node, neighbor, key
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
    """Create an undirected heterogeneous network encoded as a networkx
    MultiGraph. The graph schema is defined using edge_metapaths and
    kind_to_abbrev.
    """
    g = networkx.MultiGraph()
    g.graph.update(kwargs)
    g.graph['schema'] = schema.create_undirected_schema(edge_metapaths, kind_to_abbrev)
    return g


def filter_nodes_by_edge_kind(g, edge_kinds):
    """Remove nodes without a single edge of a specified kind.
    In the future should make method taking a node and returning a counter of
    edge kind to number
    """
    edge_kinds = set(edge_kinds)
    edge_kinds |= {(edge_kind[1], edge_kind[0], edge_kind[2])
                           for edge_kind in edge_kinds}
    nodes_to_remove = set()
    for node in g.nodes_iter():
        remove = True
        for u, v, key in g.edges(node, keys=True):
            edge_kind = g.node[u]['kind'], g.node[v]['kind'], key
            if edge in edge_kinds:
                remove = False
        if remove:
            nodes_to_remove.add(node)
    
    g.remove_nodes_from(nodes_to_remove)

def filter_edges_by_edge_kind(g, edges_to_keep):
    """Keep only edges of specified kinds and then remove unconnected nodes.
    UNTESTED
    """
    edges_to_keep = set(edges_to_keep)
    for node, neighbor, key in g.edges_iter(keys=True):
        edge_set = {(node, neighbor, key), (neighbor, node, key)}
        if not edge_set & edges_to_keep:
            g.remove_edge(node, neighbor, key)
    remove_unconnected_nodes(g)
    
def remove_unconnected_nodes(g):
    """Remove unconnected nodes"""
    unconnected_nodes = (node for node, degree in g.degree_iter() if not degree)
    g.remove_nodes_from(unconnected_nodes)        
