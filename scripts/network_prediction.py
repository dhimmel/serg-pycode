import argparse
import csv
import os

import numpy
import sklearn.linear_model

import heteronets.features

parser = argparse.ArgumentParser()
parser.add_argument('--networks-dir', type=os.path.expanduser, default=
    '~/Documents/serg/networks/')
parser.add_argument('--network-id', default='130403-1')
args = parser.parse_args()


network_dir = os.path.join(args.networks_dir, args.network_id)
feature_dir = os.path.join(network_dir, 'features')
leaning_features_path = os.path.join(network_dir, 'features', 'learning-features.txt')

feature_tuple = heteronets.features.numpy_read_features(leaning_features_path)
source, target, status, features, feature_names = feature_tuple


# logistic regression
logreg = sklearn.linear_model.LogisticRegression(C=0.0001)
logreg.fit(features, status)
logreg.score(features, status)

y_predicted = logreg.predict_proba(features)[:, 1]

svm = sklearn.svm.SVC(probability=True)
svm.fit(features, status)
svm.score(features, status)


y_predicted = svm.predict_proba(features)[:, 1]



svm = sklearn.svm.LinearSVC()
svm.fit(features, status)
svm.score(features, status)

y_predicted = svm.predict_proba(features)[:, 1]

"""
source = 'interferon beta-1a'
target = 'Multiple Sclerosis'

metapaths = g.graph['metapaths']
#metapaths.append(('drug', 'indication', 'disease'))
metapath_to_paths = get_metapath_to_paths(g, source, target=None, metapaths, cutoff=2)
print metapath_to_paths

reversed_metapaths = tuple(reversed(metapath) for metapath in metapaths)
"""
"""
pkl_path = '/home/dhimmels/Documents/serg/ipanet/ipanet.pkl'
g = networkx.read_gpickle(pkl_path)


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