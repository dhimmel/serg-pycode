import collections
import random

import networkx
import numpy
import sklearn.linear_model

import bioparser.data



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


pkl_path = '/home/dhimmels/Documents/serg/ipanet/ipanet.pkl'
g = networkx.read_gpickle(pkl_path)

# Create a dictionary of node kind to edges
kind_to_nodes = dict()
for node, data in g.nodes_iter(data=True):
    kind = data['kind']
    kind_to_nodes.setdefault(kind, set()).add(node)

# Delete improper edge kinds
valid_edge_kinds = {'target', 'disease_gene', 'indication'}
for node, neighbor, data in g.edges_iter(data=True):
    kind = data['kind']
    if kind not in valid_edge_kinds:
        g.remove_edge(node, neighbor)

# Create a dictionary of edge kind to edges
kind_to_edges = dict()
for node, neighbor, data in g.edges_iter(data=True):
    kind = data['kind']
    edge = node, neighbor
    kind_to_edges.setdefault(kind, set()).add(edge)

for key, value in kind_to_edges.items():
    print key, len(value)

# Select positives and negatives
indications = list(kind_to_edges['indication'])
num_of_positives = len(indications) / 200
positives = random.sample(indications, num_of_positives)
drugs = kind_to_nodes['drug']
diseases = kind_to_nodes['disease']
negatives = list()
while len(negatives) < num_of_positives:
    disease = (set(random.choice(indications)) & diseases).pop()
    drug = (set(random.choice(indications)) & drugs).pop()
    if not g.has_edge(disease, drug):
        edge = drug, disease
        negatives.append(edge)

# delete positives edges from the network
g.remove_edges_from(positives)
total_path_counts(g)

# Create predictor and response arrays
training_edges = negatives + positives
y = numpy.repeat([0, 1], [len(negatives), len(positives)])

npcs_by_edge = [normalized_path_counts(*edge) for edge in training_edges]

metapaths = set()
for npc in npcs_by_edge:
    metapaths |= set(npc.keys())

metapaths = list(metapaths)
metapaths.sort(key=lambda x: len(x))


X = list()
for npc in npcs_by_edge:
    x_row = list()
    for metapath in metapaths:
        value = npc.get(metapath)
        if value is None:
            value = 0.0
        x_row.append(value)
    X.append(x_row)

X = numpy.array(X)

logreg = sklearn.linear_model.LogisticRegression()
logreg.fit(X, y)
y_predicted = logreg.predict_proba(X)[:,1]

fpr, tpr, thresholds = sklearn.metrics.roc_curve(y, y_predicted)
sklearn.metrics.auc(fpr, tpr)

feature_file = '/home/dhimmels/Documents/serg/ipanet/features.txt'
numpy.savetxt(feature_file, numpy.column_stack((X, y.T)))


#print g.node['multiple sclerosis']




        




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