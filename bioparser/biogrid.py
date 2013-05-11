import csv
import collections
import os


def psimitab_reader(path):
    """ """
    f = open(path)
    f.next() # skip header line
    fieldnames = ['id_a', 'id_b', 'alt_id_a', 'alt_id_b', 'alias_a', 'alias_b',
                  'method', 'author', 'pub_id', 'taxon_a', 'taxon_a', 'type',
                  'source', 'interaction_id', 'confidence']
    reader = csv.DictReader(f, fieldnames=fieldnames, delimiter='\t')
    for row in reader:
        for fieldname, field in row.items():
            row[fieldname] = field.split('|')
        yield row
    f.close()

g = psimitab_reader('/home/dhimmels/Documents/serg/data-sources/biogrid/BIOGRID-ORGANISM-Homo_sapiens-3.2.100.mitab.txt')
interactions = list(g)

method_counter = collections.Counter(interaction['method'][0] for interaction in interactions)
for k, v in method_counter.items():
    print k.split('(')[1][:-1], '\t', v