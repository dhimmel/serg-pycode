import os
import csv

import ashg13

import bioparser.data
import bioparser.copub


data = bioparser.data.Data()

ashg_dir = os.path.expanduser('~/Documents/serg/ashg13/')

def map_efo_to_copub():
    """ """
    efo_id_to_name = data.efo.get_id_to_name()
    path = os.path.join(ashg_dir, 'efo-to-copub-mappings-auto.txt')
    mapping_file = open(path, 'w')
    writer = csv.writer(mapping_file, delimiter='\t')
    efo_ids = ashg13.calculate_disease_subset()
    fieldnames = ['efo_id', 'efo_name', 'bi_id', 'copub_name']
    writer.writerow(fieldnames)
    for efo_id in efo_ids:
        efo_name = efo_id_to_name[efo_id]
        bi_id, copub_name = bioparser.copub.get_top_keyword(efo_name, 'disease')
        row = [efo_id, efo_name, bi_id, copub_name]
        writer.writerow(row)
    mapping_file.close()

def read_efo_to_copub():
    """ """
    path = os.path.join(ashg_dir, 'efo-to-copub-mappings-manual.txt')
    mapping_file = open(path)
    reader = csv.DictReader(mapping_file, delimiter='\t')
    efo_to_copub = dict()
    for row in reader:
        efo_id = row['efo_id']
        bi_id = int(row['bi_id'])
        efo_to_copub[efo_id] = bi_id
    return efo_to_copub



#print read_efo_to_copub()

