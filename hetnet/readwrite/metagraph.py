import collections
import gzip
import json



def write_json(metagraph, path):
    """ """
    
    metanode_kinds = metagraph.node_dict.keys()
    
    metaedge_tuples = [edge.get_id() for edge in
                       metagraph.get_edges(exclude_inverts=True)]
    
    combined = collections.OrderedDict()
    combined['metanode_kinds'] = metanode_kinds
    combined['metaedge_tuples'] = metaedge_tuples
    
    if path.endswith('.gz'):
        write_file = gzip.open(path, 'w')
    else:
        write_file = open(path, 'w')
    json.dump(combined, write_file, indent=2)
    write_file.close()

def read_json(path):
    raise Exception('Incomplete')
