import os
import csv

import pandas
import pandas.io.parsers
import numpy
import numpy.core.defchararray

import reader


def csv_dict_rows(path, delimiter='\t', **kwargs):
    read_file = open(path)
    reader = csv.DictReader(read_file, delimiter=delimiter, **kwargs)
    rows = list(reader)
    read_file.close()
    return rows


def read_triformat(directory, x_filename):
    """
    X - 2d numpy array with rows representing compounds and columns represnting features
    compounds - pandas DataFrame with compound information
    features - pandas DataFrame with features information
    """
    features_path = os.path.join(directory, 'features.txt')
    features = pandas.io.parsers.read_table(features_path)

    compounds_path = os.path.join(directory, 'compounds.txt')
    compounds = pandas.io.parsers.read_table(compounds_path)

    x_path = os.path.join(directory, x_filename)
    X = numpy.loadtxt(x_path, delimiter='\t')

    triformat = {'X': X, 'compounds': compounds, 'features': features}
    return triformat

def compound_subset(triformat, bool_select):
    bool_select = numpy.array(bool_select)
    subset = dict()
    subset['X'] = triformat['X'][bool_select, :]
    subset['compounds'] = triformat['compounds'].loc[bool_select, :]
    subset['features'] = triformat['features']
    return subset

def feature_subset(triformat, bool_select):
    bool_select = numpy.array(bool_select)
    subset = dict()
    subset['X'] = triformat['X'][:, bool_select]
    subset['features'] = triformat['features'].loc[bool_select, :]
    subset['compounds'] = triformat['compounds']
    return subset



def read_features(path):
    read_file = open(path)
    csv_reader = csv.reader(read_file, delimiter='\t')
    column_names = csv_reader.next()[1:]
    feature_rows = list()
    row_names = list()
    for row in csv_reader:
        row_names.append(row[0])
        feature_rows.append(row[1:])
    read_file.close()
    X = numpy.array(feature_rows, dtype=float)
    row_names = numpy.array(row_names)
    column_names = numpy.array(column_names)
    return row_names, column_names, X


def produce_data(
        feature_path = '/home/dhimmels/Documents/serg/remyelination/networks/140117/features/SEA-pval-cutoff-1e-05.txt',
        row_filter = lambda y: y > 0,
        column_prefixes=None):

    # Read features
    row_names, column_names, X = read_features(feature_path)
    # Filter columns
    if column_prefixes:
        valid_columns = [numpy.core.defchararray.startswith(column_names, prefix) for prefix in column_prefixes]
        valid_columns = reduce(numpy.logical_or, valid_columns)
        column_names = column_names[valid_columns]
        X = X[:, valid_columns]
    #X = X.astype(int)

    # Read screen compounds
    screen_reader = reader.ScreenReader()
    positives, negatives, omitted = screen_reader.binarize_screen()
    smiles_to_compound = screen_reader.get_smiles_to_compound()

    #
    y = list()
    compounds = list()
    for smiles in row_names:
        compound = smiles_to_compound[smiles]
        compounds.append(compound)
        y.append(compound['status'])

    y = numpy.array(y)
    compounds = numpy.array(compounds)

    if row_filter:
        valid_rows = row_filter(y)
        row_names = row_names[valid_rows]
        y = y[valid_rows]
        X = X[valid_rows, :]
        compounds = compounds[valid_rows]

    return compounds, column_names, X, y

def write_dataset(directory, features_input_path):
    if not os.path.exists(directory):
        os.mkdir(directory)

    compounds, column_names, X, y = produce_data(feature_path=features_input_path, row_filter=lambda y: y != -1)

    # write compounds
    fieldnames = ['name', 'canonical_smiles', 'mean_MBP', 'SEM_MBP', 'mean_PDGFR', 'SEM_PDGFR', 'status']
    compounds_file = open(os.path.join(directory, 'compounds.txt'), 'w')
    writer = csv.DictWriter(compounds_file, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(compounds)
    compounds_file.close()

    # write features
    features_path = os.path.join(directory, 'features.txt')
    features_file = open(features_path, 'w')
    fieldnames = ['name', 'metric', 'metapath', 'target']
    writer = csv.writer(features_file, delimiter='\t')
    writer.writerow(fieldnames)
    for column_name in column_names:
        row = [column_name] + column_name.split('|')
        writer.writerow(row)
    features_file.close()

    # write X
    X_path = os.path.join(directory, 'X.txt')
    numpy.savetxt(X_path, X, delimiter='\t', fmt='%s')


def write_anonymized(path, X, y):
    XY = numpy.column_stack((y, X))
    #path = '/home/dhimmels/Documents/serg/remyelination/networks/140117/Hetero-DWPC_0.5-SEA-pval-cutoff-1e-05.txt'
    numpy.savetxt(path, XY, delimiter='\t', fmt='%s')


if __name__ == '__main__':
    network_dir = '/home/dhimmels/Documents/serg/remyelination/networks/140117/'
    directory = os.path.join(network_dir, 'SEA-1e-05_DWPC')
    input_path = os.path.join(network_dir, 'features', 'SEA-pval-cutoff-1e-05.txt')
    write_dataset(directory, input_path)