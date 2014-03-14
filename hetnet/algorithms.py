import collections
import itertools
import operator
import random

import hetnet
import hetnet.agents

        
def path_degree_product(path, damping_exponent, exclude_edges=set(), exclude_masked=True):
    """ """
    degrees = list()
    for edge in path:
        source_edges = edge.source.get_edges(edge.metaedge, exclude_masked)
        target_edges = edge.target.get_edges(edge.metaedge.inverse, exclude_masked)
        if exclude_edges:
            source_edges -= exclude_edges
            target_edges -= exclude_edges
        source_degree = len(source_edges)
        target_degree = len(target_edges)
        degrees.append(source_degree)
        degrees.append(target_degree)

    damped_degrees = [degree ** damping_exponent for degree in degrees]
    degree_product = reduce(operator.mul, damped_degrees)
    return degree_product

def PCs(paths_s):
    return len(paths_s)

def PCt(paths_t):
    return len(paths_t)

def PC(paths):
    return len(paths)

def NPC(paths_s, paths_t):
    if len(paths_t) == 0:
        paths = list()
    else:
        target = paths_t[0].source()
        paths = [path for path in paths_s if path.target() == target]
    denom = len(paths_s) + len(paths_t)
    if denom:
        return 2.0 * len(paths) / denom
    else:
        return None

def DWPC(paths, damping_exponent, exclude_edges=set(), exclude_masked=True):
    degree_products = (path_degree_product(path, damping_exponent, exclude_edges=exclude_edges, exclude_masked=exclude_masked) for path in paths)
    path_weights = (1.0 / degree_product for degree_product in degree_products)
    dwpc = sum(path_weights)
    return dwpc

def get_metrics():
    """ """
    metrics = list()

    metric = collections.OrderedDict()
    metric['name'] = 'PC'
    metric['fxn'] = PC
    metric['algorithm'] = 'PC'
    metric['arguments'] = {'paths': None}
    metrics.append(metric)

    metric = collections.OrderedDict()
    metric['name'] = 'PCs'
    metric['fxn'] = PCs
    metric['algorithm'] = 'PCs'
    metric['arguments'] = {'paths_s': None}
    metrics.append(metric)

    metric = collections.OrderedDict()
    metric['name'] = 'PCt'
    metric['algorithm'] = 'PCt'
    metric['fxn'] = PCt
    metric['arguments'] = {'paths_t': None}
    metrics.append(metric)

    metric = collections.OrderedDict()
    metric['name'] = 'NPC'
    metric['algorithm'] = 'NPC'
    metric['fxn'] = NPC
    metric['arguments'] = {'paths_s': None, 'paths_t': None}
    metrics.append(metric)

    dwpc_exponents = [x / 10.0 for x in range(0, 11)]
    for damping_exponent in dwpc_exponents:
        metric = collections.OrderedDict()
        metric['name'] = 'DWPC_{}'.format(damping_exponent)
        metric['algorithm'] = 'DWPC'
        metric['fxn'] = DWPC
        metric['arguments'] = {'paths': None, 'damping_exponent': damping_exponent, 'exclude_edges': None}
        metrics.append(metric)

    return metrics




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

def product_learning_edges(graph, metaedge, sources, targets):
    """ """
    sources = sorted(sources, key = lambda x: x.id_)
    targets = sorted(targets, key = lambda x: x.id_)

    learning_edges = list()
    group = 0
    for source, target in itertools.product(sources, targets):
        edge_id = source.id_, target.id_, metaedge.kind, metaedge.direction
        status = int(edge_id in graph.edge_dict)
        exclusions = set()
        learning_edge = collections.OrderedDict((
        ('status', status), ('group', group),
        ('edge_id', edge_id), ('exclusions', exclusions)))
        learning_edges.append(learning_edge)
        group += 1
    return learning_edges

if __name__ == '__main__':

    graph = create_example_graph()
    metapaths = graph.metagraph.extract_metapaths('gene', 'disease', 2)

    for metapath in metapaths:
        print metapath
        print features_between(graph, 'IL17', 'MS', metapath, masked=False)

