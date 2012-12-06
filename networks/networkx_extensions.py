import networkx

def get_paths(g, source, target, metapath):
    """ """
    assert g.node[source]['kind'] == metapath[0]
    assert g.node[target]['kind'] == metapath[-1]
    metapath = list(metapath)
    metapath.pop(0)
    paths = [[source]]
    while metapath:
        edge_key = metapath.pop(0)
        node_kind = metapath.pop(0)
        current_depth_paths = list()
        while paths:
            preceeding_path = paths.pop()
            node = preceeding_path[-1]
            
            for node, neighbor, key in g.edges(node, keys=True):
                neighbor_kind = g.node[neighbor]['kind']
                if key == edge_key and neighbor_kind == node_kind:
                    path = preceeding_path + [edge_key, neighbor]
                    current_depth_paths.append(path)
        paths = current_depth_paths
    paths_from_source = len(paths)
    paths = filter(lambda path: path[-1] == target, paths)
    return paths


def get_precompute_metapaths(metapaths):
    """ """
    path_heads = set()
    path_tails = set()
    for metapath in metapaths:
        metapath_length = len(metapath) / 2
        if metapath_length > 2:
            path_heads.add(tuple(metapath[:5]))
        if metapath_length > 3:
            path_tails.add(tuple(metapath[-5:]))
    return path_heads, path_tails

 







def get_kind_to_nodes(g):
    # Create a dictionary of node kind to edges
    kind_to_nodes = dict()
    for node, data in g.nodes_iter(data=True):
        kind = data['kind']
        kind_to_nodes.setdefault(kind, set()).add(node)
    return kind_to_nodes

def get_kind_to_edges(g):
    # Create a dictionary of edge kind to edges
    kind_to_edges = dict()
    for node, neighbor, key in g.edges_iter(keys=True):
        edge = node, neighbor
        kind_to_edges.setdefault(key, set()).add(edge)    
    return kind_to_edges

def print_edge_kind_counts(g):
    kind_to_edges = get_kind_to_edges(g)
    for key, value in kind_to_edges.items():
        print key, len(value)




################################################################################
################################## Deprecate ###################################

def total_path_counts(g):
    """Computes the total path counts ending and starting with each node. Saves
    the results in the data dictionary for the node under the keys:
    'ending_paths' and 'starting_paths'. Computation is dynamic to improve
    efficiency.
    """
    for node, data in g.nodes_iter(data=True):
        kind = data['kind']
        paths = collections.Counter([(kind, )])
        data['ending_paths'] = paths
        data['starting_paths'] = paths
        data['temp_ending_paths'] = collections.Counter()
        data['temp_starting_paths'] = collections.Counter()
    
    for i in range(3):
        
        for node, data in g.nodes_iter(data=True):
            kind = data['kind']
            neighbors = g.neighbors_iter(node)
            for neighbor in neighbors:
                neighbor_data = g.node[neighbor]
                
                # Ending Paths
                neighbor_ending_paths = neighbor_data['ending_paths']
                for neighbor_path, count in neighbor_ending_paths.iteritems():
                    path = list(neighbor_path)
                    path.append(kind)
                    path = tuple(path)
                    data['temp_ending_paths'][path] += count
                
                # Starting Paths            
                neighbor_starting_paths = neighbor_data['starting_paths']
                for neighbor_path, count in neighbor_starting_paths.iteritems():
                    path = [kind] + list(neighbor_path)
                    path = tuple(path)
                    data['temp_starting_paths'][path] += count
        
        for node, data in g.nodes_iter(data=True):
            data['ending_paths'] += data['temp_ending_paths']
            data['temp_ending_paths'] = collections.Counter()
            data['starting_paths'] += data['temp_starting_paths']
            data['temp_starting_paths'] = collections.Counter()
    
    # Delete temporary attibutes
    for node, data in g.nodes_iter(data=True):
        del data['temp_ending_paths']
        del data['temp_starting_paths']

def path_counts(source, target, cutoff=3):
    metapath_counter = collections.Counter()
    paths = networkx.all_simple_paths(g, source, target, cutoff)
    for path in paths:
        metapath = tuple(g.node[node]['kind'] for node in path)
        metapath_counter[metapath] += 1
    return metapath_counter

def normalized_path_counts(source, target, cutoff=3):
    numerator = path_counts(source, target, cutoff) + path_counts(target, source, cutoff)
    denomenator = g.node[source]['starting_paths'] + g.node[target]['ending_paths']
    npc = dict()
    for metapath, counts in numerator.items():
        denom = denomenator[metapath]
        if denom != 0:
            npc[metapath] = float(counts) / denom
    return npc
