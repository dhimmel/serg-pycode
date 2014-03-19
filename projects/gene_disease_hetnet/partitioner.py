import os
import csv
import itertools
import random
import gzip
import collections

import bioparser.data

data = bioparser.data.Data()


random.seed(0)
percent_testing = 0.3
min_genes_per_disease = 10

associations_path = os.path.join(data.gwas_plus.directory, 'associations.txt')
associations_file = open(associations_path)


exclude_doids = {'DOID:0050589', 'DOID:2914'} # IBD and immune system disease
associations_path = os.path.join(data.gwas_plus.directory, 'associations.txt')
associations_file = open(associations_path)
associations_reader = csv.DictReader(associations_file, delimiter='\t')
positives = set()
doid_code_to_name = dict()
for association in associations_reader:
    doid_code = association['doid_code']
    if doid_code in exclude_doids:
        continue
    positives.add((doid_code, association['symbol']))
    doid_code_to_name[doid_code] = association['doid_name']
associations_file.close()

disease_counts = collections.Counter(doid_code for doid_code, gene in positives)
doid_codes = sorted(doid_code for doid_code, count in
                    disease_counts.items() if count >= min_genes_per_disease)

genes = sorted(gene.symbol for gene in data.hgnc.get_genes()
               if gene.locus_group == 'protein-coding gene')

rows = list()
rows_pos = list()
rows_neg = list()
for doid_code, gene in itertools.product(doid_codes, genes):
    status = int((doid_code, gene) in positives)
    #part = 'test' if random.random() < percent_testing else 'train'
    doid_name = doid_code_to_name[doid_code]
    row = {'doid_code': doid_code, 'doid_name': doid_name,
           'gene': gene, 'status': status, 'part': 'train'}
    rows_status = rows_pos if status else rows_neg
    rows_status.append(row)
    rows.append(row)

for rows_status in rows_neg, rows_pos:
    k = int(round(len(rows_status) * percent_testing))
    for test_row in random.sample(rows_status, k):
        test_row['part'] = 'test'

project_dir = '/home/dhimmels/Documents/serg/gene-disease-hetnet'
partition_path = os.path.join(project_dir, 'partitions.txt.gz')
partition_file = gzip.open(partition_path, 'w')
fieldnames = ['doid_code', 'doid_name', 'gene', 'status', 'part']
writer = csv.DictWriter(partition_file, delimiter='\t', fieldnames=fieldnames)
writer.writeheader()
writer.writerows(rows)
partition_file.close()