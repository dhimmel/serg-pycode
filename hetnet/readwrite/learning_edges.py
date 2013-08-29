import ast
import collections
import csv




def write_text(learning_edges, path):
    fieldnames = ['status', 'group', 'edge_id', 'exclusions']
    with open(path, 'w') as write_file:
        writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for learning_edge in learning_edges:
            learning_edge = learning_edge.copy()
            exclusions = learning_edge['exclusions']
            learning_edge['exclusions'] = [edge.get_id() for edge in exclusions]
            writer.writerow(learning_edge)
            
    
    
def read_text(path, graph):
    """ """
    with open(path) as read_file:
        reader = csv.DictReader(read_file, delimiter='\t')
        for row in reader:
            learning_edge = collections.OrderedDict()
            for key in ('status', 'group'):
                learning_edge[key] = int(row[key])
                
            learning_edge['edge_id'] = ast.literal_eval(row['edge_id'])
            exclusions = ast.literal_eval(row['exclusions'])
            exclusions = set(graph.edge_dict[edge_id] for edge_id in exclusions)
            learning_edge['exclusions'] = exclusions
            yield learning_edge

