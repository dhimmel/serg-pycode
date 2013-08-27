import ast
import collections
import csv

import metaedge

def as_dictionaries(metapaths):
    
    for metapath in metapaths:
        dictionary = collections.OrderedDict()
        dictionary['abbreviation'] = str(metapath)
        dictionary['edges'] = str([edge.get_id() for edge in metapath.edges])
        yield dictionary


def write_text(metapaths, path):
    fieldnames = ['abbreviation', 'edges']
    row_generator = as_dictionaries(metapaths)
    with open(path, 'w') as write_file:
        writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(row_generator)

def read_text(path, metagraph):
    print path
    with open(path) as read_file:
        reader = csv.DictReader(read_file, delimiter='\t')
        metapaths = list()
        for row in reader:
            #print row['edges']
            metaedge_id_list = ast.literal_eval(row['edges'])
            metaedges = tuple(metagraph.edge_dict[metaedge_id]
                              for metaedge_id in metaedge_id_list)
            metapath = metagraph.get_metapath(metaedges)
            metapaths.append(metapath)
    return metapaths

    


