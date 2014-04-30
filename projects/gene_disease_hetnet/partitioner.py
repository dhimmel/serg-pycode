import os
import csv
import itertools
import random
import gzip
import collections

import bioparser.data

data = bioparser.data.Data()


random.seed(0)
percent_training = 0.8
min_genes_per_disease = 10

associations_path = os.path.join(data.gwas_plus.directory, 'processed', 'associations.txt')
associations_file = open(associations_path)
associations_reader = csv.DictReader(associations_file, delimiter='\t')
positives = set()
doid_code_to_name = dict()
for association in associations_reader:
    doid_code = association['doid_code']
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
    doid_name = doid_code_to_name[doid_code]
    row = {'doid_code': doid_code, 'doid_name': doid_name,
           'gene': gene, 'status': status}
    rows_status = rows_pos if status else rows_neg
    rows_status.append(row)
    rows.append(row)

for rows_status in rows_neg, rows_pos:
    n = len(rows_status)
    rindexes = range(n)
    random.shuffle(rindexes)
    for rindex, row in zip(rindexes, rows_status):
        percentile = float(rindex + 1) / n
        row['percentile'] = percentile
        row['part'] = 'train' if percentile <= percent_training else 'test'


project_dir = '/home/dhimmels/Documents/serg/gene-disease-hetnet'
partition_path = os.path.join(project_dir, 'partitions.txt.gz')
partition_file = gzip.open(partition_path, 'w')
fieldnames = ['doid_code', 'doid_name', 'gene', 'status', 'percentile', 'part']
writer = csv.DictWriter(partition_file, delimiter='\t', fieldnames=fieldnames)
writer.writeheader()
writer.writerows(rows)
partition_file.close()

partition_dir = os.path.join(project_dir, 'disease-partitions')
doid_to_writer = dict()
write_files = list()
for doid_code in doid_codes:
    path = os.path.join(partition_dir, '{}.txt.gz'.format(doid_code.replace(':', '_')))
    write_file = gzip.open(path, 'w')
    write_files.append(write_file)
    writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()
    doid_to_writer[doid_code] = writer

for row in rows:
    doid_to_writer[row['doid_code']].writerow(row)

for write_file in write_files:
    write_file.close()
