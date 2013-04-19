import collections
import random
import os

import networkx

import nxutils
import schema

def longest_matching_shortcut(metapath, shortcuts):
    """Returns the longest shortcut in shortcuts which left aligns to match
    metapath. If no metapath does not begin with any shortcut, None is returned.
    """
    metapath_len = len(metapath)
    match = None
    match_len = 0
    for shortcut in shortcuts:
        shortcut_len = len(shortcut)
        if (match_len < shortcut_len and shortcut_len <= metapath_len and
            metapath.startswith(shortcut)):
            match = shortcut
            match_len = shortcut_len
    return match

def source_to_target_node_count(g, metapath, source, use_shortcuts=True):
    """Returns a counter with nodes as keys and number of paths of type metapath
    reaching node as the count. Use path_to_nodes dictionary to enable speedup."""
    if g.node[source]['kind'] != metapath.start():
        return None

    shortcuts = g.graph.get('shortcuts')
    
    metapath_position = 0
    counter = collections.Counter()
    counter[source] = 1
    
    while metapath_position < len(metapath.tuple_) - 1:
        remaining_path = metapath[metapath_position: ]
        remaining_path = schema.MetaPath(remaining_path)
        node_counter_temp = collections.Counter()        
        
        if use_shortcuts and shortcuts and len(remaining_path) > 1:
            shortcut = longest_matching_shortcut(remaining_path, shortcuts)
            for node, count in counter.items(): 
                if shortcut not in g.node[node]['path_to_nodes_counter']: # ADDED MUST DIAGNOSE
                    print 'Error', node, remaining_path, shortcut, g.node[node]['path_to_nodes_counter'].keys()
                    raise Exception
                    continue
                nodes_counter = g.node[node]['path_to_nodes_counter'][shortcut]                    
                for destination, destination_count in nodes_counter.iteritems():
                    node_counter_temp[destination] += destination_count * count
            metapath_position += len(shortcut.tuple_) - 1

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



def path_counter(g, metapaths, source, requirements=True):
    """Count the number of paths for a target starting at source which follow
    each metapath.
    """
    requirement_dict = g.graph.get('required_source_to_targets')
    if not requirement_dict:
        requirements = False
    
    target_to_metapath_to_count = dict()
    for metapath in metapaths:
        target_to_count = source_to_target_node_count(g, metapath, source)
        
        if requirements:
            for target in requirement_dict.get(source, set()):
                target_to_count[target] += 0
        
        for target, count in target_to_count.iteritems():            
            counter = target_to_metapath_to_count.setdefault(target, collections.Counter())
            counter[metapath] = count
    return target_to_metapath_to_count


def normalized_path_counter(g, metapaths, source):
    """Compute normalized path count between source and all targets connected
    by atleast one metapath."""
    target_to_metapath_to_count = path_counter(g, metapaths, source)
    target_to_metapath_to_npc = dict()
    all_paths_source = g.node[source]['all_paths']
    source_denomenator = collections.Counter()
    for metapath in metapaths:
        source_denomenator[metapath] += all_paths_source[metapath]


    for target, metapath_to_count in target_to_metapath_to_count.iteritems():
        all_paths_target = g.node[target]['all_paths']
        denomenator = collections.Counter()
        for metapath in metapaths:
            reversed_metapath = metapath.reverse()
            denomenator[metapath] = source_denomenator[metapath] + all_paths_target[reversed_metapath]    

        metapath_to_npc = dict()
        
        for metapath, denom in denomenator.iteritems():
            numer = 2 * metapath_to_count[metapath]
            metapath_to_npc[metapath] = float(numer) / denom if denom else None
        target_to_metapath_to_npc[target] = metapath_to_npc
            
    return target_to_metapath_to_npc


def learning_edge_subset(g, num_pos, num_neg, remove_positives=False, seed=None):
    """
    Returns a tuple of (num_pos positives edges, num_neg negative edges). 
    Selected edges are of kind g.graph['edge_key'], start at the node kind
    specified by g.graph['source_kind'], and end at the node kind specified by
    g.graph['target_kind']. Negatives are randomly selected to reflect the 
    degree distribution of positives. Positives are randomly selected from all
    the edges of the specified kind. Edges selected as positives are removed
    from the graph.
    """
    assert g.graph['source_kind']
    assert g.graph['target_kind']
    assert g.graph['edge_key']
    
    random.seed(seed)
    
    kind_to_nodes = nxutils.get_kind_to_nodes(g)
    kind_to_edges = nxutils.get_kind_to_edges(g)
    source_kind = g.graph['source_kind']
    target_kind = g.graph['target_kind']
    edge_key = g.graph['edge_key']
    
    edge_kind = source_kind, target_kind, edge_key
    
    sources = kind_to_nodes[source_kind]
    targets = kind_to_nodes[target_kind]
    edge_order_key = lambda edge: edge if edge[0] in sources else (edge[1], edge[0], edge[2])
    edges = list(kind_to_edges[edge_kind])
    edges = map(edge_order_key, edges)
    
    # Check if more positives are specified than exist
    if num_pos > len(edges):
        print 'Limiting number of positive learning edges to', len(edges)
        num_pos = len(edges)
        num_neg = len(edges)
    
    # Randomly select negatices to mirror positive node degree
    negatives = set()
    while len(negatives) < num_neg:
        source = random.choice(edges)[0]
        target = random.choice(edges)[1]
        edge = source, target, edge_key
        if not g.has_edge(*edge):
            negatives.add(edge)
    negatives = list(negatives)
    # Randomly select positives
    positives = random.sample(edges, num_pos)
    if remove_positives:
        g.remove_edges_from(positives)    
    return positives, negatives


def nodes_outside_metapaths(g):
    """Return the set of nodes that do not participate in any paths of a kind
    in g.graph['metapaths'].
    """
    metapaths = g.graph['metapaths']
    submetapath_tuples = set()
    for metapath in metapaths:
        for i in range(len(metapath)):
            for metapath in metapath, metapath.reverse():
                submetapath_tuples.add(metapath.split_by_index(i, reverse_head=True))
    
    outside_nodes = set()
    for node in g.nodes_iter():
        metapaths_started = set(g.node[node]['all_paths'])
        if not any(set(t) <= metapaths_started for t in submetapath_tuples):
            outside_nodes.add(node)
    return outside_nodes

def prepare_feature_optimizations(g):
    """Adds storage intensive attributes to the graph."""
    print 'computing total path counts'
    total_path_counts(g)
    
    print 'computing shortcuts'
    shortcuts = schema.MetaPath.shortcuts_for_metapaths(g.graph['metapaths'], 2)
    compute_shortcuts(g, shortcuts)
    g.graph['prepared'] = True

def compute_shortcuts(g, shortcuts):
    """Annotate each node in the graph with a dictionary named
    path_to_nodes_counter. Dictionary keys are relevant shortcut metapaths.
    The value corresponding to a shortcut metapath key is a counter with target
    nodes as keys and number of paths following the metapath as counts."""
    print 'shortcuts passed to metapaths.compute_shortcuts', shortcuts
    for source_node in g.nodes_iter():
        path_to_nodes_counter = dict()
        for path in shortcuts:
            counter = source_to_target_node_count(g, path, source_node, use_shortcuts=False)
            # counter is None if the source_node's kind does not match the path
            if counter is not None:
                path_to_nodes_counter[path] = counter
        g.node[source_node]['path_to_nodes_counter'] = path_to_nodes_counter
    g.graph['shortcuts'] = shortcuts


def total_path_counts(g):
    """Computes the total path counts starting with each node. Stores the
    resulting collection for each node as a node attribute named 'all_paths'. 
    Dynamic computation boosts efficiency.
    """
    
    depth = g.graph['max_path_length']
    for node, data in g.nodes_iter(data=True):
        node_kind = data['kind']
        node_as_path = schema.MetaPath((node_kind, ))
        paths = collections.Counter([node_as_path, ])
        data['all_paths'] = paths
        data['temp_starting_paths'] = collections.Counter()
    
    for i in range(depth):
        
        for node, data in g.nodes_iter(data=True):
            node_kind = data['kind']
            for node, neighbor, key in g.edges(node, keys=True):
                neighbor_data = g.node[neighbor]
                                
                # Starting Paths            
                neighbor_starting_paths = neighbor_data['all_paths']
                for neighbor_path, count in neighbor_starting_paths.iteritems():
                    path = neighbor_path.append(key, node_kind, prepend = True)
                    data['temp_starting_paths'][path] += count
        
        for node, data in g.nodes_iter(data=True):
            data['all_paths'] += data['temp_starting_paths']
            data['temp_starting_paths'] = collections.Counter()
    
    # Delete temporary attibutes
    for node, data in g.nodes_iter(data=True):
        del data['temp_starting_paths']


def get_paths(g, metapath, source, target=None):
    """
    Get paths of kind metapath from the source to the target. If target is
    not specified, all paths from source of kind metpath are returned. If
    a source or target is specified which does not match the source or target
    node kind in the metapath, None is returned.
    """
    if g.node[source]['kind'] != metapath.start():
        return None
    if target and g.node[target]['kind'] != metapath.end():
        return None
    metapath_deque = collections.deque(metapath.tuple_)
    metapath_deque.popleft()
    paths = [[source]]
    while metapath_deque:
        edge_key = metapath_deque.popleft()
        node_kind = metapath_deque.popleft()
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

def get_metapath_to_paths(g, source, target=None, metapaths=[], cutoff=None):
    """Return a dictionary of metapath to paths connecting source and target
    nodes. Cutoff removes metapaths exceeding the specified length.
    """
    if cutoff is not None:
        metapaths = filter(lambda x: len(x) <= cutoff, metapaths)
    metapath_to_paths = dict()
    for metapath in metapaths:
        paths = get_paths(g, metapath, source, target)
        metapath_to_paths[metapath] = paths
    return metapath_to_paths
    


def all_simple_paths(g, source, target, metapaths, cutoff=None):
    """Return dictionary of metapath to a list of paths following that metapath
    between the source and target. Only metapaths in metapaths are returned.
    If metapaths are provided cutoff is the max metapath length. Otherwise
    cutoff must be specified and all metapaths found with length less than or
    equal to cutoff are returned.
    """
    assert metapaths is not None or cutoff is not None
    if cutoff is None:
        cutoff = max(len(metapath) for metapath in metapaths)
    all_paths = networkx.all_simple_paths(g, source, target, cutoff)
    #print all_paths
    return all_paths

if __name__ =='__main__':
    ipanet_dir = '/home/dhimmels/Documents/serg/ipanet/'
    network_id = '130116-1'
    pkl_path_prepared = os.path.join(ipanet_dir, 'networks', network_id, 'prepared-graph.pkl')
    g = networkx.read_gpickle(pkl_path_prepared)



