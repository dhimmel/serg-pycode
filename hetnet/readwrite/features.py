import csv
import collections
import ast


def write_text(features, path):
    """ """
    fieldnames = ['name', 'algorithm', 'arguments']
    with open(path, 'w') as write_file:
        writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(features)

def read_text(path):
    """ """
    features = list()
    with open(path) as readfile:
        reader = csv.DictReader(readfile, delimiter='\t')
        for row in reader:
            feature = collections.OrderedDict()
            feature['name'] = row['name']
            feature['algorithm'] = row['algorithm']
            feature['arguments'] = ast.literal_eval(row['arguments'])
            features.append(feature)
    return features


