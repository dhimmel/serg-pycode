import argparse
import csv
import json
import collections
import os
import bioparser.data

parser = argparse.ArgumentParser()
parser.add_argument('--network-dir', type=os.path.expanduser, default=
    '~/Documents/serg/gene-disease-hetnet/networks/140615-all-assoc')
args = parser.parse_args()

webdata_dir = os.path.join(args.network_dir, 'webdata')

gene_info_dir = os.path.join(webdata_dir, 'gene-info')
disease_info_dir = os.path.join(webdata_dir, 'disease-info')
feature_info_dir = os.path.join(webdata_dir, 'feature-info')
for info_dir in [gene_info_dir, disease_info_dir, feature_info_dir]:
    if not os.path.isdir(info_dir):
        os.mkdir(info_dir)


## Make gene info files
hgnc = bioparser.data.Data().hgnc
genes = hgnc.get_genes()
genes = [gene for gene in genes if gene.locus_group == 'protein-coding gene']

for gene in genes:
    json_gene = collections.OrderedDict()
    json_gene['hgnc_id'] = gene.hgnc_id
    json_gene['symbol'] = gene.symbol
    json_gene['name'] = gene.name
    json_gene['aliases'] = sorted(set(gene.previous_symbols + gene.synonyms))
    json_gene['entrez'] = gene.entrez_id or gene.entrez_id_ncbi_mapped
    json_gene['ensembl'] = gene.ensembl_id or gene.ensembl_id_ensembl_mapped
    json_gene['uniprot'] = gene.uniprot_id
    json_gene['chromosome'] = gene.chromosome
    filename = '{}.json'.format(gene.hgnc_id.replace(':', '_'))
    path = os.path.join(gene_info_dir, filename)
    with open(path, 'w') as write_file:
        json.dump(json_gene, write_file)



## Make disease info files
assoc_per_disease_path = '/home/dhimmels/Documents/serg/data-sources/gwas-catalog/140205/processed/associations-per-disease.txt'
with open(assoc_per_disease_path) as assoc_per_disease_file:
    reader = csv.DictReader(assoc_per_disease_file, delimiter='\t')
    assoc_per_disease_rows = list(reader)

doid_codes = [row['disease_code'] for row in assoc_per_disease_rows]
doid_graph = bioparser.data.Data().doid.get_graph()
for doid_code in doid_codes:
    json_doid = doid_graph.node[doid_code].copy()
    json_doid['disease_code'] = doid_code
    filename = '{}.json'.format(doid_code.replace(':', '_'))
    path = os.path.join(disease_info_dir, filename)
    with open(path, 'w') as write_file:
        json.dump(json_doid, write_file)


## Make feature description file
feature_description_path = '/home/dhimmels/Documents/serg/gene-disease-hetnet/data-integration/feature-descriptions.txt'
with open(feature_description_path) as feature_description_file:
    reader = csv.DictReader(feature_description_file, delimiter='\t')
    feature_rows = list(reader)

for feature_row in feature_rows:
    underscore_feature = feature_row['feature'].replace('|', '_')
    filename = '{}.json'.format(underscore_feature)
    path = os.path.join(feature_info_dir, filename)
    with open(path, 'w') as write_file:
        json.dump(feature_row, write_file)
