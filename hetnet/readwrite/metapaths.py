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
        writer = csv.DictWriter(write_file, delimiter='\t')
        writer.writeheader()
        writer.writerows(row_generator)

def read_text(path):
    raise Exception('Incomplete')

    


