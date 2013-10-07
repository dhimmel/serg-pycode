
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

