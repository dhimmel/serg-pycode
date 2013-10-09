
import csv

def row_generator(path):
    with open(path) as read_file:
        reader = csv.DictReader(read_file, delimiter='\t')
        for row in reader:
            yield row

def get_mapping_tuples(path, field_0, field_1):
    for row in row_generator(path):
        mapping_tuple = row[field_0], row[field_1]
        yield mapping_tuple

def get_mapping_dict(path, key_field, value_field, plural=True):
    """values are sets"""
    mapping_dict = dict()
    for row in row_generator(path):
        key = row[key_field]
        value = row[value_field]
        if plural:
            mapping_dict.setdefault(key, set()).add(value)
        else:
            mapping_dict[key] = value
    return mapping_dict