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
status_to_int = {'assoc_high': 1, 'negative': 0,
                 'assoc_low': -1, 'linked_low': -1, 'linked_high': -1}

pair_to_status = dict()
doid_code_to_name = dict()

statuses_path = os.path.join(data.gwas_plus.directory, 'processed2', 'association-statuses.txt')
statuses_file = open(statuses_path)
statuses_reader = csv.DictReader(statuses_file, delimiter='\t')
for association in statuses_reader:
    doid_code = association['doid_code']
    doid_name = association['doid_name']
    doid_code_to_name[doid_code] = doid_name
    symbol = association['gene_symbol']
    pair = doid_code, symbol
    pair_to_status[pair] = association['status']
statuses_file.close()

disease_counts = collections.Counter(doid_code for (doid_code, symbol), status in
    pair_to_status.iteritems() if status == 'assoc_high')
doid_codes = sorted(doid_code for doid_code, count in
                    disease_counts.items() if count >= min_genes_per_disease)

genes = sorted(gene.symbol for gene in data.hgnc.get_genes() if gene.coding)


status_to_rows = dict()
all_rows = list()
for pair in itertools.product(doid_codes, genes):
    doid_code, gene_symbol = pair
    status = pair_to_status.get(pair, 'negative')
    status_int = status_to_int[status]
    doid_name = doid_code_to_name[doid_code]
    row = {'doid_code': doid_code, 'doid_name': doid_name,
           'gene_symbol': gene_symbol, 'status': status, 'status_int': status_int}
    status_to_rows.setdefault(status, list()).append(row)
    all_rows.append(row)

for status, rows in status_to_rows.iteritems():
    n = len(rows)
    rindexes = range(n)
    random.shuffle(rindexes)
    for rindex, row in zip(rindexes, rows):
        percentile = float(rindex + 1) / n
        row['percentile'] = percentile
        row['part'] = 'train' if percentile <= percent_training else 'test'


project_dir = '/home/dhimmels/Documents/serg/gene-disease-hetnet'
partition_path = os.path.join(project_dir, 'partitions.txt.gz')
partition_file = gzip.open(partition_path, 'w')
fieldnames = ['doid_code', 'doid_name', 'gene_symbol', 'status', 'status_int', 'percentile', 'part']
writer = csv.DictWriter(partition_file, delimiter='\t', fieldnames=fieldnames)
writer.writeheader()
writer.writerows(all_rows)
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

for row in all_rows:
    doid_to_writer[row['doid_code']].writerow(row)

for write_file in write_files:
    write_file.close()
