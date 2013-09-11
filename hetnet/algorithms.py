import collections
import operator
import random

import hetnet
import hetnet.agents

        
def path_degree_product(path, damping_exponent, exclude_masked=True):
    """ """
    degrees = list()
    for edge in path:
        source_degree = len(edge.source.get_edges(edge.metaedge, exclude_masked))
        target_degree = len(edge.target.get_edges(edge.metaedge.inverse, exclude_masked))
        degrees.append(source_degree)
        degrees.append(target_degree)

    damped_degrees = [degree ** damping_exponent for degree in degrees]
    degree_product = reduce(operator.mul, damped_degrees)
    return degree_product

def PCs(paths):
    return len(paths['paths_s'])

def PCt(paths):
    return len(paths['paths_t'])

def PC(paths):
    return len(paths['paths_st'])

def NPC(paths):
    denom = len(paths['paths_s']) + len(paths['paths_t'])
    if denom:
        return 2.0 * len(paths['paths_st']) / denom
    else:
        return None

def DWPC(paths, damping_exponent):
    paths = paths['paths_st']
    degree_products = (path_degree_product(path, damping_exponent) for path in paths)
    path_weights = (1.0 / degree_product for degree_product in degree_products)
    dwpc = sum(path_weights)
    return dwpc

algorithm_to_path_types = {
    'PCs': {'paths_s'},
    'PCt': {'paths_t'},
    'PC': {'paths_st'},
    'NPC': {'paths_st', 'paths_s', 'paths_t'},
    'DWPC': {'paths_st'}}

algorithm_to_function = {
    'PCs': PCs,
    'PCt': PCt,
    'PC': PC,
    'NPC': NPC,
    'DWPC': DWPC}
    

def features_betweens(graph, source, target, metapaths, features, exclude_edges=set()):
    results = collections.OrderedDict()
    for metapath in metapaths:
        feature_dict = features_between(graph, source, target, metapath, features, exclude_edges)
        for name, value in feature_dict.iteritems():
            key = '{}:{}'.format(name, metapath)
            results[key] = value
    return results

def features_between(graph, source, target, metapath, features, exclude_edges=set()):
    path_types = set()
    
    for feature in features:
        algorithm = feature['algorithm']
        path_types |= algorithm_to_path_types[algorithm]
        
    paths = paths_between(graph, source, target, metapath, path_types, exclude_edges)
    
    results = collections.OrderedDict()
    for feature in features:
        function = algorithm_to_function[feature['algorithm']]
        results[feature['name']] = function(paths, **feature['arguments'])
    return results

def paths_between(graph, source, target, metapath, path_types, exclude_edges):
    paths = dict()
    
    if 'paths_s' in path_types:
        paths['paths_s'] = graph.paths_from(source, metapath, exclude_edges=exclude_edges)
        
    if 'paths_t' in path_types:
        paths['paths_t'] = graph.paths_from(target, metapath.inverse, exclude_edges=exclude_edges)
    
    if 'paths_st' in path_types:
        
        if 'paths_s' in path_types:
            paths['paths_st'] = [path for path in paths['paths_s'] if path.target() == target]
            
        elif 'paths_t' in path_types:
            paths_ts = [path for path in paths['paths_t'] if path.target() == source]
            paths['paths_st'] = [hetnet.Path(path.inverse_edges) for path in paths_ts]
            
        else:
            paths['paths_st'] = graph.paths_between(source, target, metapath, exclude_edges=exclude_edges)
        
    return paths

def get_features():
    """ """
    features = list()
    
    feature = collections.OrderedDict()
    feature['name'] = 'PC'
    feature['algorithm'] = 'PC'
    feature['arguments'] = dict()
    features.append(feature)

    feature = collections.OrderedDict()
    feature['name'] = 'PCs'
    feature['algorithm'] = 'PCs'
    feature['arguments'] = dict()
    features.append(feature)

    feature = collections.OrderedDict()
    feature['name'] = 'PCt'
    feature['algorithm'] = 'PCt'
    feature['arguments'] = dict()
    features.append(feature)

    feature = collections.OrderedDict()
    feature['name'] = 'NPC'
    feature['algorithm'] = 'NPC'
    feature['arguments'] = dict()
    features.append(feature)
    
    dwpc_exponents = [x / 10.0 for x in range(0, 11)]
    for damping_exponent in dwpc_exponents:
        feature = collections.OrderedDict()
        feature['name'] = 'DWPC_{}'.format(damping_exponent)
        feature['algorithm'] = 'DWPC'
        feature['arguments'] = {'damping_exponent': damping_exponent}
        features.append(feature)
    
    return features
    

def create_example_graph():
    
    metaedges = [('gene', 'disease', 'association', 'both'),
             ('gene', 'gene', 'function', 'both'),
             ('gene', 'tissue', 'expression', 'both'),
             ('disease', 'tissue', 'pathology', 'both'),
             ('gene', 'gene', 'transcription', 'forward')]
    metagraph = hetnet.MetaGraph.from_edge_tuples(metaedges)

    graph = hetnet.Graph(metagraph)
    graph.add_node('IL17', 'gene')
    graph.add_node('BRCA1', 'gene')
    graph.add_node('ANAL1', 'gene')
    graph.add_node('MS', 'disease')
    als = graph.add_node('ALS', 'disease')
    #als.masked = True
    graph.add_node('brain', 'tissue')
    
    graph.add_edge('MS', 'IL17', 'association', 'both')
    graph.add_edge('ALS', 'IL17', 'association', 'both')
    graph.add_edge('MS', 'brain', 'pathology', 'both')
    graph.add_edge('IL17', 'brain', 'expression', 'both')
    graph.add_edge('IL17', 'BRCA1', 'transcription', 'backward')
    return graph
    

def matched_negatives(positives, source_negatives=1, target_negatives=1, seed=0):
    """ positives is a list of edges"""
    random.seed(seed)
    
    positive_ids = set(positive.get_id() for positive in positives)
    learning_edges = collections.OrderedDict()
    
    for i, positive in enumerate(positives):
        positive_id = positive.get_id()
        source_id, target_id, kind, direction = positive_id
        excluded_edges = {positive, positive.inverse}
        
        learning_edge = collections.OrderedDict((
            ('status', 1), ('group', i),
            ('edge_id', positive_id), ('exclusions', excluded_edges)))
        learning_edges[positive_id] = learning_edge
        
        # Generate negatives with the same source as the positive
        sources_matched = 0
        while sources_matched < source_negatives:
            new_target = random.choice(positives).target
            negative_id = (source_id, new_target.id_, kind, direction)
            if negative_id in learning_edges or negative_id in positive_ids:
                continue
            edges = list(new_target.edges[positive.metaedge.inverse])
            exclude_edge = random.choice(edges)
            negative_edge_exlusions = {exclude_edge, exclude_edge.inverse}
            
            learning_edge = collections.OrderedDict((
                ('status', 0), ('group', i),
                ('edge_id', negative_id),
                ('exclusions', excluded_edges | negative_edge_exlusions)))
            learning_edges[negative_id] = learning_edge
            sources_matched += 1

        # Generate negatives with the same target as the positive
        targets_matched = 0
        while targets_matched < target_negatives:
            new_source = random.choice(positives).source
            negative_id = (new_source.id_, target_id, kind, direction)
            if negative_id in learning_edges or negative_id in positive_ids:
                continue
            edges = list(new_source.edges[positive.metaedge])
            exclude_edge = random.choice(edges)
            negative_edge_exlusions = {exclude_edge, exclude_edge.inverse}

            learning_edge = collections.OrderedDict((
                ('status', 0), ('group', i),
                ('edge_id', negative_id),
                ('exclusions', excluded_edges | negative_edge_exlusions)))
            learning_edges[negative_id] = learning_edge
            targets_matched += 1
    
    learning_edges = learning_edges.values()
    return learning_edges


    

if __name__ == '__main__':

    graph = create_example_graph()
    metapaths = graph.metagraph.extract_metapaths('gene', 'disease', 2)

    for metapath in metapaths:
        print metapath
        print features_between(graph, 'IL17', 'MS', metapath, masked=False)

