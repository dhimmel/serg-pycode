import collections
import csv
import os
import random
import shutil
import string
import sys

import networkx

import schema
import metapaths
import nxutils


valid_file_chars = set("-_.()[]{}&#@!~+= %s%s" % (string.ascii_letters, string.digits))
def string_to_valid_filename(s):
    """Return the string with invalid filename characters excluded."""
    return filter(lambda x: x in valid_file_chars, s)
    

def edge_status(g, source, target):
    """Returns the status for an edge.
    0 - Negative evaluation edge
    1 - Positive evaluation edge
    2 - Non existent network edge
    3 - Existent network edge
    """
    edge = source, target, g.graph['edge_key']
    if edge in g.graph['negatives']:
        status = 0
    elif edge in g.graph['positives']:
        status = 1
    else:
        status = int(g.has_edge(*edge)) + 2
    return status

def feature_generator(g, edge_to_exclusions):
    """
    Generates features (only NPC for now) based on the graph specifications.
    """
    feature_paths = g.graph['metapaths']
    prediction_metapath = schema.MetaPath((g.graph['source_kind'], g.graph['edge_key'], g.graph['target_kind']))
    if prediction_metapath in feature_paths:
        feature_paths.remove(prediction_metapath)
    
    kind_to_nodes = nxutils.get_kind_to_nodes(g)

    edge_to_exclusions = collections.OrderedDict(sorted(edge_to_exclusions.items()))

    for edge, exclusions in edge_to_exclusions.iteritems():
        source, target, edge_key = edge
        metapath_to_metric_dict = metapaths.features_for_metapaths(
            g, source, target, edge_key, exclusions, feature_paths)
        feature_dict = metapaths.flatten_feature_dict(metapath_to_metric_dict)
        combined_dict = collections.OrderedDict()
        combined_dict['source'] = source
        #combined_dict['source_name'] = g.node[source]['name']
        combined_dict['target'] = target
        combined_dict['target_name'] = g.node[target]['name']
        combined_dict['status'] = edge_status(g, source, target)
        combined_dict.update(feature_dict)
        yield combined_dict

def write_features(g, edge_to_exclusions, path):
    """
    """
    schema = g.graph['schema']
    feature_file = open(path, 'w')
    
    initialize_writer = True
    for feature_dict in feature_generator(g, edge_to_exclusions):
        if initialize_writer:
            fieldnames = feature_dict.keys()
            dict_writer = csv.DictWriter(feature_file, fieldnames=fieldnames, delimiter='\t')
            dict_writer.writeheader()
            initialize_writer = False
        dict_writer.writerow(feature_dict)
        print feature_dict['source'], '\t', feature_dict['target']
    feature_file.close()

def write_partitioned_features(g, feature_dir):
    """
    Not UPDATED
    """
    source_kind = g.graph['source_kind']
    target_kind = g.graph['target_kind']
    source_dir = os.path.join(feature_dir, source_kind)
    target_dir = os.path.join(feature_dir, target_kind)
    
    # Remove source and target directories since files within get appended to.
    for dir_path in source_dir, target_dir:
        if os.path.isdir(dir_path):
            shutil.rmtree(dir_path)
    
    # Create directories
    for dir_path in feature_dir, source_dir, target_dir:
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)
    
    schema = g.graph['schema']
    feature_paths = g.graph['metapaths']
    fieldnames = ['source', 'target', 'status'] + map(str, feature_paths)

    learning_path = os.path.join(feature_dir, 'learning-features.txt')
    learning_file = open(learning_path, 'w')
    learning_writer = csv.DictWriter(learning_file, fieldnames=fieldnames, delimiter='\t')
    learning_writer.writeheader()
    
    features_path = os.path.join(feature_dir, 'features.txt')
    features_file = open(features_path, 'w')
    features_writer = csv.DictWriter(features_file, fieldnames=fieldnames, delimiter='\t')
    features_writer.writeheader()
    
    initiated_sources = set()
    initiated_targets = set()
    
    for feature_dict in feature_generator(g):

        # Write edge to features.txt
        features_writer.writerow(feature_dict)
        
        # Write learning edge
        status = feature_dict['status']
        if status == 0 or status == 1:
            learning_writer.writerow(feature_dict)
        
        # Write to source and target files
        for source_or_target, directory, initiated in (
            ['source', source_dir, initiated_sources],
            ['target', target_dir, initiated_targets]):
            
            name = string_to_valid_filename(feature_dict[source_or_target])
            path = os.path.join(directory, name + '.txt')
            f = open(path, 'a')
            dict_writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            if name not in initiated:
                dict_writer.writeheader()
                initiated.add(name)
            dict_writer.writerow(feature_dict)
            f.close()
    
    learning_file.close()
    features_file.close()


def read_features_dict(path):
    """Return a generator of dictionaries representing rows of a feature file."""
    # type_dict specifies the type conversion to be applied. Each key denotes
    # a column name and the value is the conversion. Columns not included are
    # converted to floats.
    type_dict = {'source': str, 'target': str, 'status': int}
    with open(path) as feature_file:
        reader = csv.DictReader(feature_file, delimiter='\t')
        for row in reader:
            yield {key: type_dict.get(key, float)(value) for key, value in row.items()}
                

def read_features_list(path):
    """ """
    # type_dict specifies the type conversion to be applied. Each key denotes
    # a column name and the value is the conversion. Columns not included are
    # converted to floats.
    type_dict = {'source': str, 'target': str, 'status': int}
    with open(path) as feature_file:
        reader = csv.reader(feature_file, delimiter='\t')
        fieldnames = reader.next()
        conversions = tuple(type_dict.get(field, float) for field in fieldnames)
        
        for row in reader:
            yield map(lambda conv, elem: conv(elem), conversions, row)
            
def numpy_read_features(path):
    """Return a tuple of (source, target, status, features)."""
    import numpy
    # read table as a structured array (each row is a tuple)
    feature_array = numpy.genfromtxt(path, delimiter='\t', names=True, dtype=None)
    source = feature_array['source']
    target = feature_array['target']
    status = feature_array['status']
    feature_names = numpy.array(feature_array.dtype.names[3: ])
    features = feature_array[feature_names]
    # convert from structured array to normal ndarray
    features = features.view((numpy.float, len(features.dtype.names)))
    return source, target, status, features, feature_names


def subset_feature_file(name, subset_name, filter_function):
    """ """
    path = os.path.join(ipanet_dir, 'networks', network_id, 'features', name + '.txt')
    feature_file = open(path)
    dict_reader = csv.DictReader(feature_file, delimiter='\t')

    subset_path = os.path.join(ipanet_dir, 'networks', network_id, subset_name + '.txt')
    subset_file = open(subset_path, 'w')
    dict_writer = csv.DictWriter(subset_file, delimiter='\t',
                                 fieldnames=dict_reader.fieldnames)
    dict_writer.writeheader()
    for row in dict_reader:
        if filter_function(row):
            dict_writer.writerow(row)
    feature_file.close()
    subset_file.close()
    print 'writing', subset_name, 'complete'
