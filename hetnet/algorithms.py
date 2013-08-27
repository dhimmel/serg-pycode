import collections
import operator
import random

import hetnet


def degree_weighted_path_count(paths, damping_exponent=0.5):
    """ """
    degree_products = (path_degree_product(path, damping_exponent) for path in paths)
    path_weights = (1.0 / degree_product for degree_product in degree_products)
    dwpc = sum(path_weights)
    return dwpc
        
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

def features_between(graph, source, target, metapath, **kwargs):
    """ """
    if not isinstance(source, hetnet.Node):
        source = graph.node_dict[source]
    if not isinstance(target, hetnet.Node):
        target = graph.node_dict[target]

    feature_dict = collections.OrderedDict()
    
    # compute normalized path count
    paths_source = graph.paths_from(source, metapath, **kwargs)
    print paths_source
    paths_target = graph.paths_from(target, metapath.inverse, **kwargs)
    PCs = len(paths_source)
    PCt = len(paths_target)
    paths = [path for path in paths_source if path.target() == target]
    PC = len(paths)
    assert PC == len([path for path in paths_target if path.target() == source])
    NPC = 2.0 * PC / (PCs + PCt) if PCs + PCt else None
    feature_dict['PC'] = PC
    feature_dict['PCs'] = PCs
    feature_dict['PCt'] = PCt
    feature_dict['NPC'] = NPC
    
    # compute degree weighted path count
    dwpc_exponents = [x / 10.0 for x in range(0, 11)]
    for dwpc_exponent in dwpc_exponents:
        feature_name = 'DWPC_' + str(dwpc_exponent)
        feature_dict[feature_name] = degree_weighted_path_count(paths, dwpc_exponent)
    
    return feature_dict

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
    
    negative_to_exclusions = dict()
    positive_to_exclusions = dict()
    
    for positive in positives:
        positive_id = positive.get_id()
        source_id, target_id, kind, direction = positive_id
        excluded_edges = {positive, positive.inverse}
        positive_to_exclusions[positive_id] = excluded_edges
        
        # Generate negatives with the same source as the positive
        sources_matched = 0
        while sources_matched < source_negatives:
            new_target = random.choice(positives).target
            negative_id = (source_id, new_target.id_, kind, direction)
            if (negative_id in positive_to_exclusions or 
                negative_id in negative_to_exclusions):
                continue
            edges = list(new_target.edges[positive.metaedge.inverse])
            exclude_edge = random.choice(edges)
            negative_edge_exlusions = {exclude_edge, exclude_edge.inverse}
            negative_to_exclusions[negative_id] = excluded_edges | negative_edge_exlusions
            sources_matched += 1

        # Generate negatives with the same target as the positive
        targets_matched = 0
        while targets_matched < target_negatives:
            new_source = random.choice(positives).source
            negative_id = (new_source.id_, target_id, kind, direction)
            if (negative_id in positive_to_exclusions or 
                negative_id in negative_to_exclusions):
                continue
            edges = list(new_source.edges[positive.metaedge])
            exclude_edge = random.choice(edges)
            negative_edge_exlusions = {exclude_edge, exclude_edge.inverse}
            negative_to_exclusions[negative_id] = excluded_edges | negative_edge_exlusions
            targets_matched += 1
    
    return positive_to_exclusions, negative_to_exclusions


if __name__ == '__main__':

    graph = create_example_graph()
    metapaths = graph.metagraph.extract_metapaths('gene', 'disease', 2)

    for metapath in metapaths:
        print metapath
        print features_between(graph, 'IL17', 'MS', metapath, masked=False)

