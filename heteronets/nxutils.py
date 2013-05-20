import collections

import networkx

import schema

"""
Making a networkx graph that conforms to the standards used in heteronets:
The function create_undirected_network creates a networkx MultiGraph instance.
To create a graph, first call this function and then add nodes and edges to
the returned graph. Each node must be given a node kind which is stored under
the the 'kind' key of that nodes data dictionary. Node kind is accessible for a
single node using g.node[node_id]['kind']. Each edge must have a kind which is
stored under the 'key' attribute of that edge. Edges are represented as tuples
formatted as (source_node_id, target_node_id, edge_kind). An edge describing
an edge kind is formatted as (source_kind, target_kind, edge_kind). All longer
metapaths are formatted exclusively as metapath objects.


Graph Attributes

schema
kind_to_nodes
kind_to_edges
abbrev_to_path
"""

def get_kind_to_nodes(g, refresh=True):
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

def get_kind_to_edges(g, refresh=True):
    """Create a dictionary of edge kind to edges. The dictionary is returned
    and stored as a graph attribute named 'kind_to_edges'. refresh specifies
    whether the dictionary should be recomputed if it already exists."""
    kind_to_edges = g.graph.get('kind_to_edges')
    if not refresh and kind_to_edges:
        return kind_to_edges
    kind_to_edges = dict()
    
    for node, neighbor, key in g.edges_iter(keys=True):
        node_kind = g.node[node]['kind']
        neighbor_kind = g.node[neighbor]['kind']
        
        edge_kind = node_kind, neighbor_kind, key
        edge = node, neighbor, key
        kind_to_edges.setdefault(edge_kind, set()).add(edge)
        
        edge_kind = neighbor_kind, node_kind, key
        edge = neighbor, node, key
        kind_to_edges.setdefault(edge_kind, set()).add(edge)
        
    g.graph['kind_to_edges'] = kind_to_edges 
    return kind_to_edges

def print_node_kind_counts(g):
    kind_to_nodes = get_kind_to_nodes(g)
    for key, value in kind_to_nodes.items():
        print key, len(value)

def print_edge_kind_counts(g):
    kind_to_edges = get_kind_to_edges(g)
    printed_edge_kinds = set()
    for edge_kind, edges in kind_to_edges.items():
        edge_set_kind = frozenset({edge_kind[0], edge_kind[1]}), edge_kind[2]
        if not edge_set_kind in printed_edge_kinds:
            printed_edge_kinds.add(edge_set_kind)
            print edge_kind, len(edges)

def create_undirected_network(edge_metapaths, kind_to_abbrev=None, **kwargs):
    """Create an undirected heterogeneous network encoded as a networkx
    MultiGraph. The graph schema is defined using edge_metapaths and
    kind_to_abbrev.
    """
    g = networkx.MultiGraph()
    g.graph.update(kwargs)
    g.graph['schema'] = schema.create_undirected_schema(edge_metapaths, kind_to_abbrev)
    return g


def node_degree_counter(g, node):
    """Returns a Counter object with edge_kind tuples as keys and the number
    of edges with the specified edge_kind incident to the node as counts.
    """
    degree_counter = collections.Counter()
    for node, neighbor, key in g.edges(node, keys=True):
        node_kind = g.node[node]['kind']
        neighbor_kind = g.node[neighbor]['kind']
        edge_kind = node_kind, neighbor_kind, key
        degree_counter[edge_kind] += 1
    return degree_counter


def filter_nodes_by_edge_kind(g, edge_kinds):
    """Remove nodes without a single edge of a specified kind."""
    edge_kinds = set(edge_kinds)
    edge_kinds |= {(edge_kind[1], edge_kind[0], edge_kind[2])
                           for edge_kind in edge_kinds}
    nodes_to_remove = set()
    for node in g.nodes_iter():
        if not set(node_degree_counter(g, node)) & edge_kinds:
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

def export_as_gml(g, path):
    """save g as a .gml file after modifications."""
    g = g.copy()
    for key in g.graph.keys():
        del g.graph[key]
    for node, data in g.nodes_iter(data=True):
        for key in data.keys():
            if key not in {'kind', 'name'}:
                del data[key]
        
    for node, neighbor, key in g.edges(keys=True):
        g.edge[node][neighbor][key]['kind'] = key
    networkx.write_gml(g, path)
    
def write_gpickle(g, path):
    kind_to_abbrev = schema.MetaPath.kind_to_abbrev
    g.graph['kind_to_abbrev'] = kind_to_abbrev
    networkx.write_gpickle(g, path)
    
def read_gpickle(path):
    g = networkx.read_gpickle(path)
    kind_to_abbrev = g.graph['kind_to_abbrev']
    schema.MetaPath.kind_to_abbrev = kind_to_abbrev
    return g

