import argparse
import csv
import os
import collections
import itertools

import networkx

import heteronets.nxutils
import heteronets.schema
import heteronets.metapaths
import heteronets.features

def get_parser_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--project-dir', type=os.path.expanduser, default=
        '~/Documents/serg/ashg13/')
    parser.add_argument('--network-id', required=True)
    parser.add_argument('--learning-edges', action='store_true')
    parser.add_argument('--MS-edges', action='store_true')
    parser.add_argument('--diseases', nargs='*')
    #parser.add_argument('--description', required=True)
    args = parser.parse_args()
    return args


def select_positives(g):
    kind_to_node = heteronets.nxutils.get_kind_to_nodes(g)

    # Find all diseases with expression and association edges
    diseases = kind_to_node['disease']
    learning_diseases = set()
    for disease in diseases:
        degree_counter = heteronets.nxutils.node_degree_counter(g, disease)
        if ('disease', 'gene', 'association') not in degree_counter:
            continue
        if ('disease', 'tissue', 'pathology') not in degree_counter:
            continue
        learning_diseases.add(disease)
    
    # Take all associations edges from those nodes as positives
    positives = g.edges(learning_diseases, keys=True)
    positives = filter(lambda e: e[2] == 'association', positives)
    positives = sorted((e[1], e[0], e[2]) for e in positives)
    return positives



if __name__ == '__main__':
    args = get_parser_args()
    network_dir = os.path.join(args.project_dir, args.network_id)

    prepared_pkl_path = os.path.join(network_dir, 'prepared-graph.pkl')    
    if not os.path.exists(prepared_pkl_path):
        raw_pkl_path = os.path.join(network_dir, 'raw-graph.pkl')
        g = heteronets.nxutils.read_gpickle(raw_pkl_path)

        # Select positives
        g.graph['source_kind'] = 'gene'
        g.graph['target_kind'] = 'disease'
        g.graph['edge_key'] = 'association'
        positives = select_positives(g)
        print len(positives), 'positives'
        
        positive_to_exclusions, negative_to_exclusions = heteronets.metapaths.matched_negatives(positives)
        g.graph['positive_to_exclusions'] = positive_to_exclusions
        g.graph['negative_to_exclusions'] = negative_to_exclusions
        g.graph['negatives'] = set(negative_to_exclusions)
        g.graph['positives'] = set(positive_to_exclusions)
        
        
        schema = g.graph['schema']
        g.graph['metapaths'] = heteronets.schema.extract_metapaths(
            g.graph['schema'], g.graph['source_kind'], g.graph['target_kind'],
            max_length = 3, exclude_all_source_target_edges = False)
    
        print 'Writing prepared graph to pkl.'
        heteronets.nxutils.write_gpickle(g, prepared_pkl_path)
    else:
        g = heteronets.nxutils.read_gpickle(prepared_pkl_path)

    # Compute and learning features
    if args.learning_edges:
        positive_to_exclusions = g.graph['positive_to_exclusions']
        negative_to_exclusions = g.graph['negative_to_exclusions']
        edge_to_exclusions = dict(negative_to_exclusions.items() + positive_to_exclusions.items())
        path = os.path.join(network_dir, 'features', 'learning-features.txt')
        heteronets.features.write_features(g, edge_to_exclusions, path)

    
    # Compute and save disease features
    if args.diseases:
        genes = heteronets.nxutils.get_kind_to_nodes(g)['gene']
        for disease in args.diseases:
            name = g.node[disease]['name']
            print disease, name
            pairs = list(itertools.product(genes, [disease]))
            edge_to_exclusions = dict()
            for pair in pairs:
                edge = (pair[0], pair[1], 'association')
                exclusions = {edge, (edge[1], edge[0], edge[2])}
                edge_to_exclusions[edge] = exclusions    
            path = os.path.join(network_dir, 'features', 'disease', name + '.txt')
            heteronets.features.write_features(g, edge_to_exclusions, path)
