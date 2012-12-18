import networkx
import collections

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

def path_counter(g, metapaths, source, target=None, shortcuts=None):
    """Count the number of paths between source and target for desired
    metapaths.
    BUGGGGG: target to source doesn't work. reversed metapaths.
    """
    counter = collections.Counter()
    for metapath in metapaths:
        target_to_counts = source_to_target_node_count(g, metapath, source, shortcuts)
        #if target_to_counts is None:
        #    continue
        #print source, target
        #print target_to_counts
        if target:
            count = target_to_counts[target]
        else:
            count = sum(target_to_counts.itervalues())
        counter[metapath] = count
    return counter


def normalized_path_counter(g, metapaths, source, target, shortcuts=None):
    numerator = path_counter(g, metapaths, source, target, shortcuts)
    denomenator = collections.Counter()
    all_paths_source = g.node[source]['all_paths']
    all_paths_target = g.node[source]['all_paths']
    for metapath in metapaths:
        denomenator[metapath] += all_paths_source[metapath]
        reversed_metapath = tuple(reversed(metapath))
        denomenator[metapath] += all_paths_target[reversed_metapath]    
        
    npc = dict()
    for metapath, counts in numerator.items():
        denom = denomenator[metapath]
        npc[metapath] = float(counts) / denom if denom else 0.0
    return npc



def shortcuts_for_metapaths(metapaths, shortcut_len=2):
    """ """
    shorcuts = set()
    for metapath in metapaths:
        num_nodes = len(metapath) / 2
        depth = 0
        while depth + shortcut_len <= num_nodes:
            start = depth * 2
            end = start + shortcut_len * 2 + 1
            shortcut = metapath[start:end]
            shortcut = tuple(shortcut)
            shorcuts.add(shortcut)
            depth += shortcut_len
    return shorcuts

def compute_shortcuts(g, shortcuts):
    """ """
    for source_node in g.nodes_iter():
        path_to_nodes_counter = dict()
        for path in shortcuts:
            counter = source_to_target_node_count(g, path, source_node)
            if counter is not None:
                path_to_nodes_counter[path] = counter
        g.node[source_node]['path_to_nodes_counter'] = path_to_nodes_counter


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

def total_path_counts(g, depth):
    """Computes the total path counts ending and starting with each node. Saves
    the results in the data dictionary for the node under the keys:
    'ending_paths' and 'starting_paths'. Computation is dynamic to improve
    efficiency.
    """
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

###############################################################################

################################################################################
################################## Deprecate ###################################

def get_paths(g, metapath, source, target=None):
    """
    Get paths of kind metapath from the source to the target. If target is
    not specified, all paths from source of kind metpath are returned. If
    a source or target is specified which does not match the source or target
    node kind in the metapath, None is returned.
    """
    if g.node[source]['kind'] != metapath[0]:
        return None
    if target and g.node[target]['kind'] != metapath[-1]:
        return None
    metapath = collections.deque(metapath)
    metapath.popleft()
    paths = [[source]]
    while metapath:
        edge_key = metapath.popleft()
        node_kind = metapath.popleft()
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
    if target:
        paths = filter(lambda path: path[-1] == target, paths)
    return paths


def get_metapath_shortcuts_old(metapaths):
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



def compute_path_to_nodes_old(g, metapaths):
    for source_node in g.nodes_iter():
        path_to_nodes_counter = dict() #dict.fromkeys(paths, list())
        for path in paths:
            all_paths = get_paths(g, path, source_node)
            if all_paths is None:
                continue
            source_nodes = [one_path[-1] for one_path in all_paths]
            path_to_nodes_counter[path] = source_nodes
        g.node[source_node]['path_to_nodes'] = path_to_nodes
        print source_node




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
