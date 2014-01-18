import os
import csv

import numpy
import numpy.core.defchararray

import reader



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


def produce_data(column_prefixes=None):

    # Read features
    path = '/home/dhimmels/Documents/serg/remyelination/networks/140117/features/SEA-pval-cutoff-1e-05.txt'
    row_names, column_names, X = read_features(path)
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

    valid_rows = y != -1
    row_names = row_names[valid_rows]
    y = y[valid_rows]
    X = X[valid_rows, :]
    compounds = compounds[valid_rows]
    return compounds, column_names, X, y


def write_anonymized(path, X, y):
    XY = numpy.column_stack((y, X))
    #path = '/home/dhimmels/Documents/serg/remyelination/networks/140117/Hetero-DWPC_0.5-SEA-pval-cutoff-1e-05.txt'
    numpy.savetxt(path, XY, delimiter='\t', fmt='%s')