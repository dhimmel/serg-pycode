import os
import csv

import bioparser.data
import bioparser.copub
import mapping.manual_reader

data = bioparser.data.Data()

project_dir = os.path.expanduser('~/Documents/serg/gene-disease-hetnet/')
mapping_dir = os.path.join(project_dir, 'data-integration', 'mappings')
copub_dir = os.path.join(project_dir, 'data-integration', 'copub')


def automap_bto_to_copub_gnf():
    """ """

    gnf_tissues = set(x['bto_id'] for x in data.gnf.expression_generator())
    animal_tissues = set(data.bto.get_animal_tissues())
    gnf_animal_tissues = gnf_tissues & animal_tissues


    gnf_animal_tissues = sorted(gnf_animal_tissues)

    bto_graph = data.bto.get_graph()

    bioparser.copub.load_client()

    path = os.path.join(mapping_dir, 'bto-to-copub-mappings-gnf-auto.txt')
    mapping_file = open(path, 'w')
    writer = csv.writer(mapping_file, delimiter='\t')

    fieldnames = ['bto_id', 'bto_name', 'bi_id', 'copub_name']
    writer.writerow(fieldnames)
    for bto_id in gnf_animal_tissues:
        bto_name = bto_graph.node[bto_id]['name']
        bi_id, copub_name = bioparser.copub.get_top_keyword(bto_name, 'tissue')
        row = [bto_id, bto_name, bi_id, copub_name]
        writer.writerow(row)
    mapping_file.close()


def map_doid_to_copub():
    """ """
    bioparser.copub.load_client()
    
    path = '/home/dhimmels/Documents/serg/data-mapping/manual/efo_doid/gwas-pairs-editted.tsv'
    gcat_mapped_doid_rows = list(mapping.manual_reader.row_generator(path))
    
    path = os.path.join(mapping_dir, 'doid-to-copub-mappings-auto.txt')
    mapping_file = open(path, 'w')
    writer = csv.writer(mapping_file, delimiter='\t')

    fieldnames = ['doid_id', 'doid_name', 'bi_id', 'copub_name']
    writer.writerow(fieldnames)
    for doid_row in gcat_mapped_doid_rows:
        doid_id = doid_row['doid_id']
        doid_name = doid_row['doid_name']
        bi_id, copub_name = bioparser.copub.get_top_keyword(doid_name, 'disease')
        row = [doid_id, doid_name, bi_id, copub_name]
        writer.writerow(row)
    mapping_file.close()

def read_doid_to_copub():
    """ """
    path = os.path.join(mapping_dir, 'doid-to-copub-mappings-manual.txt')
    mapping_file = open(path)
    reader = csv.DictReader(mapping_file, delimiter='\t')
    doid_to_copub = dict()
    for row in reader:
        doid_id = row['doid_id']
        bi_id = int(row['bi_id'])
        doid_to_copub[doid_id] = bi_id
    return doid_to_copub

def read_bto_to_copub():
    """ """
    path = os.path.join(mapping_dir, 'bto-to-copub-mappings-gnf-manual.txt')
    mapping_file = open(path)
    reader = csv.DictReader(mapping_file, delimiter='\t')
    bto_to_copub = dict()
    for row in reader:
        bto_id = row['bto_id']
        bi_id = row['bi_id']
        if not bi_id:
            continue
        bi_id = int(bi_id)
        bto_to_copub[bto_id] = bi_id
    return bto_to_copub


def compute_doid_bto_cooccurrence():
    """ """
    bioparser.copub.load_client()
    
    bto_graph = data.bto.get_graph()
    doid_graph = data.doid.get_graph()

    bto_to_copub = read_bto_to_copub()
    doid_to_copub = read_doid_to_copub()
    copub_to_bto = {v: k for k, v in bto_to_copub.items()}
    
    copub_to_doid = {v: k for k, v in doid_to_copub.items()}
    
    copub_to_doids = dict()
    for doid_id, copub_id in doid_to_copub.items():
        copub_to_doids.setdefault(copub_id, list()).append(doid_id)
    
    
    copub_tissues = copub_to_bto.keys()
    copub_diseases = copub_to_doids.keys()
    cooccurrence_gen = bioparser.copub.term_set_cooccurrences(
        copub_diseases, copub_tissues)
    rows = list()
    for copub_disease, copub_tissue, r_scaled in cooccurrence_gen:
        doid_ids = copub_to_doids[copub_disease]
        bto_id = copub_to_bto[copub_tissue]
        bto_name = bto_graph.node[bto_id]['name']
        for doid_id in doid_ids:
            doid_name = doid_graph.node[doid_id]['name']
            row = doid_id, doid_name, bto_id, bto_name, r_scaled
            rows.append(row)
            print row

    rows.sort(key=lambda x: (x[1], x[3]))
    fieldnames = ['doid_id', 'doid_name', 'bto_id', 'bto_name', 'r_scaled']
    path = os.path.join(copub_dir, 'doid-bto-cooccurrences.txt')
    write_file = open(path, 'w')
    writer = csv.writer(write_file, delimiter='\t')
    writer.writerow(fieldnames)
    writer.writerows(rows)
    write_file.close()


def doid_bto_cooccurrence_generator():
    """ """
    path = os.path.join(copub_dir, 'doid-bto-cooccurrences-gnf-tissues.txt')
    read_file = open(path)
    reader = csv.DictReader(read_file, delimiter='\t')
    for row in reader:
        row['r_scaled'] = float(row['r_scaled'])
        yield row
    read_file.close()





def read_efo_to_copub():
    """ """
    path = os.path.join(project_dir, 'efo-to-copub-mappings-manual.txt')
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
    path = os.path.join(project_dir, 'tiger-tissue-to-copub-mappings-auto.txt')
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
    path = os.path.join(project_dir, 'tiger-tissue-to-copub-mappings-manual.txt')
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
    path = os.path.join(project_dir, 'tiger-tissue-to-disease-cooccurrences.txt')
    write_file = open(path, 'w')
    writer = csv.writer(write_file, delimiter='\t')
    writer.writerow(fieldnames)
    writer.writerows(rows)
    write_file.close()

def tiger_efo_cooccurrence_generator():
    """ """
    path = os.path.join(project_dir, 'tiger-tissue-to-disease-cooccurrences.txt')
    read_file = open(path)
    reader = csv.DictReader(read_file, delimiter='\t')
    for row in reader:
        row['r_scaled'] = float(row['r_scaled'])
        yield row
    read_file.close()
    


if __name__ == '__main__':
    
    compute_doid_bto_cooccurrence()
    
    #automap_bto_to_copub_gnf()
    #map_tiger_to_copub()
    #print read_tiger_to_copub()
    #compute_tiger_efo_cooccurrence()

