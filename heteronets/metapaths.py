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

def compute_pathcount_upper_bounds(g):
    """Computes the total path counts starting with each node. Stores the
    resulting collection for each node as a node attribute named 'all_paths'. 
    Dynamic computation boosts efficiency.
    """
    
    depth = g.graph['max_path_length']
    for node, data in g.nodes_iter(data=True):
        node_kind = data['kind']
        node_as_path = schema.MetaPath((node_kind, ))
        paths = collections.Counter([node_as_path, ])
        data['PC_upper_bounds'] = paths
        data['PC_upper_bounds_temp'] = collections.Counter()
    
    for i in range(depth):
        
        for node, data in g.nodes_iter(data=True):
            node_kind = data['kind']
            for node, neighbor, key in g.edges(node, keys=True):
                neighbor_data = g.node[neighbor]
                                
                # Starting Paths            
                neighbor_starting_paths = neighbor_data['PC_upper_bounds']
                for neighbor_path, count in neighbor_starting_paths.iteritems():
                    path = neighbor_path.append(key, node_kind, prepend = True)
                    data['PC_upper_bounds_temp'][path] += count
        
        for node, data in g.nodes_iter(data=True):
            data['PC_upper_bounds'] += data['PC_upper_bounds_temp']
            data['PC_upper_bounds_temp'] = collections.Counter()
    
    # Delete temporary attibutes
    for node, data in g.nodes_iter(data=True):
        del data['PC_upper_bounds_temp']

def paths_from_source(g, source, metapath, excluded_edges=set()):
    """Returns a list of paths starting at source of kind metapath. Duplicate
    nodes and edges in excluded_edges are excluded. A path is a series of
    nodes separated by edge keys."""
    if g.node[source]['kind'] != metapath[0]:
        return None
    paths = [[source]]
    node_index = 2
    max_metapath_index = len(metapath) * 2
    while node_index <= max_metapath_index:
        metapath_node_kind = metapath[node_index]
        metapath_edge_key = metapath[node_index - 1]
        valid_paths = list()
        for path in paths:
            node = path[-1]
            for node, neighbor, key in g.edges(node, keys=True):
                if key != metapath_edge_key:
                    continue
                neighbor_kind = g.node[neighbor]['kind']
                if neighbor_kind != metapath_node_kind:
                    continue
                if neighbor in path[::2]:
                    continue
                if (node, neighbor, key) in excluded_edges or (neighbor, node, key) in excluded_edges:
                    continue
                valid_paths.append(path + [key, neighbor])
        paths = valid_paths
        node_index += 2
    return paths

def features_for_metapath(g, source, target, edge_key, metapath):
    """ """
    feature_dict = collections.OrderedDict()
    excluded_edges = {(source, target, edge_key), (target, source, edge_key)}
    source_paths = paths_from_source(g, source, metapath, excluded_edges)
    target_paths = paths_from_source(g, target, metapath.reverse(), excluded_edges)
    PCs = len(source_paths)
    PCt = len(target_paths)

    paths_st = filter(lambda path: path[-1] == target, source_paths)
    paths_ts = filter(lambda path: path[-1] == source, target_paths)
    PCst = len(paths_st)
    PCts = len(paths_ts)
    assert PCst == PCts
    PC = PCst
    NPC_denominator = PCs + PCt
    NPC = 2.0 * PC / NPC_denominator if NPC_denominator else None
    feature_dict['PC'] = PC
    feature_dict['PCs'] = PCs
    feature_dict['PCt'] = PCt
    feature_dict['NPC'] = NPC
    return feature_dict

def features_for_metapaths(g, source, target, edge_key, metapaths):
    """ """
    metapath_to_metric_dict = collections.OrderedDict()
    for metapath in metapaths:
        metapath_to_metric_dict[metapath] = features_for_metapath(
            g, source, target, edge_key, metapath)
    return metapath_to_metric_dict

def flatten_feature_dict(metapath_to_metric_dict):
    feature_dict = collections.OrderedDict()
    for metapath, metric_dict in metapath_to_metric_dict.iteritems():
        for metric, value in metric_dict.iteritems():
            key = metric + '_' + str(metapath)
            feature_dict[key] = value
    return feature_dict

if __name__ =='__main__':
    ipanet_dir = '/home/dhimmels/Documents/serg/ipanet/'
    network_id = '130116-1'
    pkl_path_prepared = os.path.join(ipanet_dir, 'networks', network_id, 'prepared-graph.pkl')
    g = networkx.read_gpickle(pkl_path_prepared)



