import collections
import random

import networkx

def longest_matching_shortcut(metapath, shortcuts):
    """Returns the longest shortcut in shortcuts which left aligns to match
    metapath. If no metapath does not begin with any shortcut, None is returned.
    """
    metapath_len = len(metapath)
    match = None
    match_len = 0
    for shortcut in shortcuts:
        shortcut_len = len(shortcut)
        if (shortcut_len <= metapath_len
            and all(shortcut[i] == metapath[i] for i in xrange(shortcut_len))
            and match_len < shortcut_len):
            match = shortcut
            match_len = shortcut_len
    return match

def source_to_target_node_count(g, metapath, source, shortcuts=None):
    """Returns a counter with nodes as keys and number of paths of type metatype
    reaching node as the count. Use path_to_nodes dictionary to enable speedup."""
    if g.node[source]['kind'] != metapath[0]:
        return None

    metapath_position = 0
    counter = collections.Counter()
    counter[source] = 1
    
    while metapath_position < len(metapath) - 1:
        remaining_path = metapath[metapath_position: ]
        node_counter_temp = collections.Counter()        
        
        if len(remaining_path) > 3 and shortcuts:
            shortcut = longest_matching_shortcut(remaining_path, shortcuts)
            for node, count in counter.items():
                nodes_counter = g.node[node]['path_to_nodes_counter'][shortcut]
                for destination, destination_count in nodes_counter.iteritems():
                    node_counter_temp[destination] += destination_count * count
            metapath_position += len(shortcut) - 1
        
        else:
            node_kind = remaining_path[2]            
            edge_key = remaining_path[1]
            for node, count in counter.iteritems():
                for node, neighbor, key in g.edges(node, keys=True):
                    neighbor_kind = g.node[neighbor]['kind']
                    if key == edge_key and neighbor_kind == node_kind:
                        node_counter_temp[neighbor] += count   
            metapath_position += 2
            
        counter = node_counter_temp
        
    return counter



def path_counter(g, metapaths, source, shortcuts=None):
    """Count the number of paths between source and target for desired
    metapaths.
    """
    target_to_metapath_to_count = dict()
    for metapath in metapaths:
        target_to_count = source_to_target_node_count(g, metapath, source, shortcuts)
        for target, count in target_to_count.iteritems():
            counter = target_to_metapath_to_count.setdefault(target, collections.Counter())
            counter[metapath] = count
    return target_to_metapath_to_count


def normalized_path_counter(g, metapaths, source, shortcuts=None):
    """Compute normalized path count between source and all targets connected
    by atleast one metapath."""
    target_to_metapath_to_count = path_counter(g, metapaths, source, shortcuts)
    target_to_metapath_to_npc = dict()
    all_paths_source = g.node[source]['all_paths']
    source_denomenator = collections.Counter()
    for metapath in metapaths:
        source_denomenator[metapath] += all_paths_source[metapath]


    for target, metapath_to_count in target_to_metapath_to_count.iteritems():
        all_paths_target = g.node[target]['all_paths']
        denomenator = collections.Counter()
        for metapath in metapaths:
            reversed_metapath = tuple(reversed(metapath))
            denomenator[metapath] = source_denomenator[metapath] + all_paths_target[reversed_metapath]    

        metapath_to_npc = dict()
        
        for metapath, denom in denomenator.iteritems():
            numer = metapath_to_count[metapath]
            metapath_to_npc[metapath] = float(numer) / denom if denom else None
        target_to_metapath_to_npc[target] = metapath_to_npc
            
    return target_to_metapath_to_npc


def learning_edge_subset(g, num_pos, num_neg):
    """
    Returns a tuple of (num_pos positives edges, num_neg negative edges). 
    Selected edges are of kind g.graph['edge_kind'], start at the node kind
    specified by g.graph['source_kind'], and end at the node kind specified by
    g.graph['target_kind']. Negatives are randomly selected to reflect the 
    degree distribution of positives. Positives are randomly selected from all
    the edges of the specified kind. Edges selected as positives are removed
    from the graph.
    """
    assert g.graph['source_kind']
    assert g.graph['target_kind']
    assert g.graph['edge_kind']
    
    kind_to_nodes = get_kind_to_nodes(g)
    kind_to_edges = get_kind_to_edges(g)
    source_kind = g.graph['source_kind']
    target_kind = g.graph['target_kind']
    edge_kind = g.graph['edge_kind']
    sources = kind_to_nodes[source_kind]
    targets = kind_to_nodes[target_kind]
    edge_order_key = lambda edge: edge if edge[0] in sources else (edge[1], edge[0])
    edges = list(kind_to_edges[edge_kind])
    edges = map(edge_order_key, edges)
    # Randomly select negatices to mirror positive node degree
    negatives = list()
    while len(negatives) < num_neg:
        source = random.choice(edges)[0]
        target = random.choice(edges)[1]
        if not g.has_edge(source, target):
            edge = source, target
            negatives.append(edge)
    # Randomly select positives
    positives = random.sample(edges, num_pos)
    g.remove_edges_from(positives)    
    return positives, negatives

def prepare_for_feature_computation(g, max_path_length, edge_kind_tuple,
                                    num_pos, num_neg):
    """
    Compute metapaths, shortcuts, positives edges, and negative edges. Modify
    network by removing positive edges. Annotate nodes with computed shortcuts
    and path count counters.
    """
    source_kind, edge_kind, target_kind = edge_kind_tuple
    g.graph['source_kind'] = source_kind
    g.graph['target_kind'] = target_kind
    g.graph['edge_kind'] = edge_kind
    
    g.graph['max_path_length'] = max_path_length
    
    g.graph['metapaths'] = g.graph['schema'].metapaths(source_kind, target_kind, max_path_length)
    g.graph['shortcuts'] = shortcuts_for_metapaths(g.graph['metapaths'], 2)

    g.graph['positives'], g.graph['negatives'] = learning_edge_subset(g, num_pos, num_neg)
    print 'computing shortcuts'
    compute_shortcuts(g)
    print 'computing total path counts'
    total_path_counts(g)
    g.graph['prepared'] = True

    
def shortcuts_for_metapaths(metapaths, shortcut_length):
    """Compute desired shortcuts for faster computation of metapaths."""
    shorcuts = set()
    for metapath in metapaths:
        num_nodes = len(metapath) / 2
        depth = 0
        while depth + shortcut_length <= num_nodes:
            start = depth * 2
            end = start + shortcut_length * 2 + 1
            shortcut = metapath[start:end]
            shortcut = tuple(shortcut)
            shorcuts.add(shortcut)
            depth += shortcut_length
    return shorcuts

def compute_shortcuts(g):
    """Annotate each node in the graph with a dictionary named
    path_to_nodes_counter. Dictionary keys are relevant shortcut metapaths.
    The value corresponding to a shortcut metapath key is a counter with target
    nodes as keys and number of paths following the metapath as counts."""
    shortcuts = g.graph['shortcuts']
    for source_node in g.nodes_iter():
        path_to_nodes_counter = dict()
        for path in shortcuts:
            counter = source_to_target_node_count(g, path, source_node)
            if counter is not None:
                path_to_nodes_counter[path] = counter
        g.node[source_node]['path_to_nodes_counter'] = path_to_nodes_counter


def get_kind_to_nodes(g):
    """Create a dictionary of node kind to edges"""
    kind_to_nodes = dict()
    for node, data in g.nodes_iter(data=True):
        kind = data['kind']
        kind_to_nodes.setdefault(kind, set()).add(node)
    return kind_to_nodes

def get_kind_to_edges(g):
    """Create a dictionary of edge kind to edges."""
    kind_to_edges = dict()
    for node, neighbor, key in g.edges_iter(keys=True):
        edge = node, neighbor
        kind_to_edges.setdefault(key, set()).add(edge)    
    return kind_to_edges

def print_edge_kind_counts(g):
    kind_to_edges = get_kind_to_edges(g)
    for key, value in kind_to_edges.items():
        print key, len(value)

def total_path_counts(g):
    """Computes the total path counts ending and starting with each node. Saves
    the results in the data dictionary for the node under the keys:
    'ending_paths' and 'starting_paths'. Computation is dynamic to improve
    efficiency.
    """
    depth = g.graph['max_path_length']
    for node, data in g.nodes_iter(data=True):
        kind = data['kind']
        paths = collections.Counter([(kind, )])
        #data['ending_paths'] = paths
        data['all_paths'] = paths
        #data['temp_ending_paths'] = collections.Counter()
        data['temp_starting_paths'] = collections.Counter()
    
    for i in range(depth):
        
        for node, data in g.nodes_iter(data=True):
            kind = data['kind']
            neighbors = g.neighbors_iter(node)
            for node, neighbor, key in g.edges(node, keys=True):
                neighbor_data = g.node[neighbor]
                
                """
                # Ending Paths
                neighbor_ending_paths = neighbor_data['ending_paths']
                for neighbor_path, count in neighbor_ending_paths.iteritems():
                    path = list(neighbor_path)
                    path.append(key)
                    path.append(kind)
                    path = tuple(path)
                    data['temp_ending_paths'][path] += count
                """
                
                # Starting Paths            
                neighbor_starting_paths = neighbor_data['all_paths']
                for neighbor_path, count in neighbor_starting_paths.iteritems():
                    path = [kind, key] + list(neighbor_path)
                    path = tuple(path)
                    data['temp_starting_paths'][path] += count
        
        for node, data in g.nodes_iter(data=True):
            #data['ending_paths'] += data['temp_ending_paths']
            #data['temp_ending_paths'] = collections.Counter()
            data['all_paths'] += data['temp_starting_paths']
            data['temp_starting_paths'] = collections.Counter()
    
    # Delete temporary attibutes
    for node, data in g.nodes_iter(data=True):
        #del data['temp_ending_paths']
        del data['temp_starting_paths']
