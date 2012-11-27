import collections

import networkx

import bioparser.data

pkl_path = '/home/dhimmels/Documents/serg/ipanet/ipanet.pkl'
g = networkx.read_gpickle(pkl_path)


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
    
    # Dele
    for node, data in g.nodes_iter(data=True):
        del data['temp_ending_paths']
        del data['temp_starting_paths']

total_path_counts(g)
#print g.node['multiple sclerosis']


kind_to_nodes = dict()
for node, data in g.nodes_iter(data=True):
    kind = data['kind']
    kind_to_nodes.setdefault(kind, set()).add(node)


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
        




source = 'multiple sclerosis'
target = 'interferon beta-1a'
normalized_path_counts(source, target)



path_counts = path_counts(source, target)

#print compute_metapath_counter(source, source)
#print compute_metapath_counter(target, target)






#print networkx.info(g)
#node_kind_counts = collections.Counter(data['kind'] for node, data in g.nodes_iter(data=True))
#print node_kind_counts


#print 'Number of connected components:', networkx.number_connected_components(g)


"""
for node, data in g.nodes_iter(data=True):
    print node
    print data
    print g.neighbors(node)
"""    







################################################################################
################################# Network Stats ################################
#execfile('create_ipanet.py')

if __name__ == '__main__':
    pass