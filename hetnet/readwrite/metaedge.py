import ast

def get_representation_dict(metaedge):
    dictionary = collections.OrderedDict()
    #dictionary['abbreviation'] = str(metaedge)
    dictionary['path'] = str(metaedge)
    dictionary['edge_id'] = str(metaedge.get_id())
    return dictionary


def write_text(metaedge, path):
    """ """
    representation_dict = get_representation_dict(metaedge)
    with open(path, 'w') as write_file:
        file.write(representation_dict['edge_id'])
    