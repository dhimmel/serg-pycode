import os
import csv

ashg_dir = os.path.expanduser('~/Documents/serg/ashg13/')

def get_mappings(key, value):
    path = os.path.join(ashg_dir, 'disease-mappings.txt')
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        map_dict = {row[key]: row[value] for row in reader}
    return map_dict
    