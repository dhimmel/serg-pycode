import ast
import collections

def get_representation_dict(metaedge):
    dictionary = collections.OrderedDict()
    #dictionary['abbreviation'] = str(metaedge)
    dictionary['path_str'] = str(metaedge)
    dictionary['edge_id'] = str(metaedge.get_id())
    return dictionary


def write_text(metaedge, path):
    """ """
    representation_dict = get_representation_dict(metaedge)
    with open(path, 'w') as write_file:
        write_file.write(representation_dict['edge_id'])

def read_text(metagraph, path):
    """ """
    with open(path) as read_file:
        line = read_file.readline().rstrip()
    metaedge_id = ast.literal_eval(line)
    metaedge = metagraph.get_edge(metaedge_id)
    return metaedge

