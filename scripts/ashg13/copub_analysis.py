import os
import csv

import bioparser.data
import bioparser.copub


data = bioparser.data.Data()

ashg_dir = os.path.expanduser('~/Documents/serg/ashg13/')

def map_efo_to_copub():
    """ ashg13 module has been removed. function will not work"""
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

def map_tiger_to_copub():
    """ """
    path = os.path.join(ashg_dir, 'tiger-tissue-to-copub-mappings-auto.txt')
    mapping_file = open(path, 'w')
    writer = csv.writer(mapping_file, delimiter='\t')
    tiger_tissues = data.tiger.get_tissues()
    fieldnames = ['tiger_tissue', 'bi_id', 'copub_name']
    writer.writerow(fieldnames)
    for tiger_tissue in tiger_tissues:
        modified_tiger_tissue = tiger_tissue.replace('_', ' ')
        bi_id, copub_name = bioparser.copub.get_top_keyword(modified_tiger_tissue, 'tissue')
        row = [tiger_tissue, bi_id, copub_name]
        writer.writerow(row)
    mapping_file.close()

def read_tiger_to_copub():
    """ """
    path = os.path.join(ashg_dir, 'tiger-tissue-to-copub-mappings-manual.txt')
    mapping_file = open(path)
    reader = csv.DictReader(mapping_file, delimiter='\t')
    tiger_to_copub = dict()
    for row in reader:
        tiger_tissue = row['tiger_tissue']
        bi_id = int(row['bi_id'])
        tiger_to_copub[tiger_tissue] = bi_id
    return tiger_to_copub

def compute_tiger_efo_cooccurrence():
    """ """
    tiger_to_copub = read_tiger_to_copub()
    efo_to_copub = read_efo_to_copub()
    copub_to_tiger = {v: k for k, v in tiger_to_copub.items()}
    copub_to_efo = {v: k for k, v in efo_to_copub.items()}
    copub_tissues = tiger_to_copub.values()
    copub_diseases = efo_to_copub.values()    
    
    cooccurrence_gen = bioparser.copub.term_set_cooccurrences(
        copub_diseases, copub_tissues)
    rows = list()
    efo_id_to_name = data.efo.get_id_to_name()
    for copub_disease, copub_tissue, r_scaled in cooccurrence_gen:
        efo_disease_id = copub_to_efo[copub_disease]
        efo_disease_name = efo_id_to_name[efo_disease_id]
        tiger_tissue = copub_to_tiger[copub_tissue]
        row = efo_disease_id, efo_disease_name, tiger_tissue, r_scaled
        rows.append(row)
        print row
    rows.sort(key=lambda x: (x[1], x[3]))
    fieldnames = ['efo_disease_id', 'efo_disease_name', 'tiger_tissue', 'r_scaled']
    path = os.path.join(ashg_dir, 'tiger-tissue-to-disease-cooccurrences.txt')
    write_file = open(path, 'w')
    writer = csv.writer(write_file, delimiter='\t')
    writer.writerow(fieldnames)
    writer.writerows(rows)
    write_file.close()

def tiger_efo_cooccurrence_generator():
    """ """
    path = os.path.join(ashg_dir, 'tiger-tissue-to-disease-cooccurrences.txt')
    read_file = open(path)
    reader = csv.DictReader(read_file, delimiter='\t')
    for row in reader:
        row['r_scaled'] = float(row['r_scaled'])
        yield row
    read_file.close()
    


if __name__ == '__main__':
    #map_tiger_to_copub()
    #print read_tiger_to_copub()
    compute_tiger_efo_cooccurrence()

