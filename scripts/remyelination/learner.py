import collections
import pprint

import numpy
import sklearn
import sklearn.svm
import sklearn.metrics
import sklearn.linear_model
import sklearn.cross_validation
import sklearn.preprocessing
import sklearn.naive_bayes
import sklearn.ensemble
import sklearn.neighbors

import feature_reader

path = '/home/dhimmels/Documents/serg/remyelination/networks/140106/Anon-DWPC_0.5-SEA-pval-cutoff-1e-05.txt'
#path = '/home/dhimmels/Documents/serg/remyelination/networks/140106/Anon-PC-SEA-pval-cutoff-1e-05.txt'


#XY = numpy.loadtxt(path, delimiter='\t', skiprows=1)
#y = XY[:, 0].astype(int) - 1
#X = XY[:, 1:]
#X_full = X

column_prefixes = ['DWPC_0.5|D-t-P|', 'DWPC_0.5|D-t-P-p-G-m-C2|']
compounds, column_names, X_full, y = feature_reader.produce_data(column_prefixes=column_prefixes)
X = X_full
print 'X shape', X.shape

names_known0 = ['DWPC_0.5|D-t-P|sp_P10827', 'DWPC_0.5|D-t-P|sp_P10828', 'DWPC_0.5|D-t-P|sp_Q92731']
names_known1 = ['DWPC_0.5|D-t-P|sp_P11229', 'DWPC_0.5|D-t-P|sp_P10827', 'DWPC_0.5|D-t-P|sp_P10828', 'DWPC_0.5|D-t-P|sp_Q92731']
indices_known0 = numpy.where(numpy.in1d(column_names, names_known0))[0]
indices_known1 = numpy.where(numpy.in1d(column_names, names_known1))[0]

"""
keep_names = ['DWPC_0.5|D-t-P|sp_P11229', 'DWPC_0.5|D-t-P|sp_P10827', 'DWPC_0.5|D-t-P|sp_P10828', 'DWPC_0.5|D-t-P|sp_Q92731']
keep_names = ['DWPC_0.5|D-t-P|sp_P11229', 'DWPC_0.5|D-t-P|sp_P10827',
              'DWPC_0.5|D-t-P|sp_P10828', 'DWPC_0.5|D-t-P|sp_Q92731',
              #'DWPC_0.5|D-t-P|sp_P09960', # Leukotriene A4 hydrolase
              #'DWPC_0.5|D-t-P-p-G-m-C5|MSigDB_c5.all:AMINE_RECEPTOR_ACTIVITY',
              'DWPC_0.5|D-t-P-p-G-i-G|ARF1',
              'DWPC_0.5|D-t-P-p-G-i-G|ARF6',
              'DWPC_0.5|D-t-P-p-G-m-C2|MSigDB_c2.all:BILANGES_SERUM_AND_RAPAMYCIN_SENSITIVE_GENES',
              'DWPC_0.5|D-t-P-p-G-m-C5|MSigDB_c5.all:NEUROTRANSMITTER_RECEPTOR_ACTIVITY',
              'DWPC_0.5|D-t-P-p-G-m-C5|MSigDB_c5.all:ST_ERK1_ERK2_MAPK_PATHWAY',
              'DWPC_0.5|D-t-P-p-G-m-C5|MSigDB_c5.all:ACETYLCHOLINE_BINDING',
              'DWPC_0.5|D-t-P-p-G-m-C2|MSigDB_c2.all:FRASOR_TAMOXIFEN_RESPONSE_UP'
              ]


indices = numpy.where(numpy.in1d(column_names, keep_names))[0]
X = X_full[:, indices]
#indices = numpy.array([617, 619, 876, 1543]) - 1
#X = X[:, indices]
"""

sk_folds = sklearn.cross_validation.StratifiedKFold(y=y, n_folds=7)

feature_select = True
if feature_select:
    feature_counter = collections.Counter()
    feature_counter_RLR = collections.Counter()
    cv_features = list()
    feature_selection_C = 0.55
    print 'feature_selection_C: {}'.format(feature_selection_C)

scale = True
aucs = list()
for train_index, test_index in sk_folds:
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    if scale:
        scaler = sklearn.preprocessing.StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)
    if feature_select:
        randomized_logistic = sklearn.linear_model.RandomizedLogisticRegression(C=feature_selection_C)
        randomized_logistic.fit(X_train, y_train)
        rlr_indices = numpy.where(randomized_logistic.scores_ != 0)[0]
        feature_counter_RLR.update(column_names[rlr_indices])
        selected_indices = numpy.array(sorted(set(rlr_indices) | set(indices_known1)))
        features = column_names[selected_indices]
        feature_counter.update(features)
        cv_features.append(features)
        X_train = X_train[:, selected_indices]
        print '{} Selected Features'.format(len(features))
        X_test = X_test[:, selected_indices]
    #clf = sklearn.svm.SVC(C=0.0001, kernel='linear', class_weight='auto')
    clf = sklearn.svm.LinearSVC(C=0.5, class_weight='auto')
    #clf = sklearn.linear_model.SGDClassifier(class_weight='auto')
    #clf = sklearn.naive_bayes.BernoulliNB()
    #clf = sklearn.naive_bayes.GaussianNB()
    #clf = sklearn.ensemble.RandomForestClassifier(max_depth=4)
    #clf = sklearn.neighbors.KNeighborsClassifier(3)
    clf.fit(X_train, y_train)
    if hasattr(clf, 'decision_function'):
        scores = clf.decision_function(X_test)
    else:
        scores = clf.predict_proba(X_test)[:, 1]
    #decisions = svc.decision_function(X_test_scaled)
    #probability_array = svc.predict_proba(X_test_scaled)
    #probabilities = probability_array[:, 1]
    #print numpy.column_stack((y_test, probabilities))
    fpr, tpr, thresholds = sklearn.metrics.roc_curve(y_test, scores)
    roc_auc = sklearn.metrics.auc(fpr, tpr)
    aucs.append(roc_auc)
    print("Area under the ROC curve : %f" % roc_auc)

print 'AUC', numpy.mean(aucs)
print 'feature counter'
pprint.pprint(feature_counter)
print 'feature_counter_RLR'
pprint.pprint(feature_counter_RLR)

#numpy.median(aucs)
