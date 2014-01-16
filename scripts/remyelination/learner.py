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

path = '/home/dhimmels/Documents/serg/remyelination/networks/140106/Anon-DWPC_0.5-SEA-pval-cutoff-1e-05.txt'
#path = '/home/dhimmels/Documents/serg/remyelination/networks/140106/Anon-PC-SEA-pval-cutoff-1e-05.txt'


XY = numpy.loadtxt(path, delimiter='\t', skiprows=1)
y = XY[:, 0].astype(int) - 1
X = XY[:, 1:]




indices = numpy.array([617, 619, 876, 1543]) - 1
X = X[:, indices]

sk_folds = sklearn.cross_validation.StratifiedKFold(y=y, n_folds=10)

feature_select = True
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
        randomized_logistic = sklearn.linear_model.RandomizedLogisticRegression()
        randomized_logistic.fit(X_train, y_train)
        selection = randomized_logistic.scores_ != 0
        X_train = X_train[:, selection]
        print X_train.shape
        X_test = X_test[:, selection]
    #clf = sklearn.svm.SVC(C=0.0001, kernel='linear', class_weight='auto')
    clf = sklearn.svm.LinearSVC(C=0.1, class_weight='auto')
    #clf = sklearn.linear_model.SGDClassifier(class_weight='auto')
    #clf = sklearn.naive_bayes.BernoulliNB()
    #clf = sklearn.naive_bayes.GaussianNB()
    #clf = sklearn.ensemble.RandomForestClassifier(max_depth=4)
    #clf = sklearn.neighbors.KNeighborsClassifier(3)
    clf.fit(X_train, y_train)
    if hasattr(clf, "decision_function"):
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
    #print("Area under the ROC curve : %f" % roc_auc)

numpy.mean(aucs)
#numpy.median(aucs)



for train_index, test_index in sk_folds:
    print train_index, test_index