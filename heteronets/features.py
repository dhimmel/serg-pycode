import collections
import csv
import os
import random
import shutil
import string
import sys

import networkx

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

def feature_generator(g):
    """
    Generates features (only NPC for now) based on the graph specifications.
    """
    schema = g.graph['schema']
    feature_paths = g.graph['metapaths']
    kind_to_nodes = nxutils.get_kind_to_nodes(g)
    sources = kind_to_nodes[g.graph['source_kind']]
    for source in sources:
        print source
        target_to_metapath_to_npc = metapaths.normalized_path_counter(g, feature_paths, source)
        for target, metapath_to_npc in target_to_metapath_to_npc.iteritems():
            feature_dict = dict.fromkeys(feature_paths, None)
            feature_dict.update(metapath_to_npc)
            feature_dict = {str(key): value for key, value in feature_dict.items()}
            feature_dict['source'] = source
            feature_dict['target'] = target
            status = edge_status(g, source, target)
            feature_dict['status'] = status
            yield feature_dict    


def write_features(g, path, status_subset = None):
    """
    If status_subset is None all edges are written despite status.
    """
    schema = g.graph['schema']
    feature_paths = g.graph['metapaths']
    feature_file = open(feature_file_path, 'w')
    fieldnames = ['source', 'target', 'status'] + map(str, feature_paths)
    dict_writer = csv.DictWriter(feature_file, fieldnames=fieldnames, delimiter='\t')
    dict_writer.writeheader()

    for feature_dict in feature_generator(g):
        if status_subset is None or feature_dict['status'] in status_subset:
            dict_writer.writerow(feature_dict)
    feature_file.close()

def write_partitioned_features(g, feature_dir):
    """
    If status_subset is None all edges are written despite status.
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
    feature_names = list(feature_array.dtype.names[3: ])
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
