import csv

import bioparser.data
import mapping.bioportal

efo = bioparser.data.Data().efo
doid = bioparser.data.Data().doid

efo_ontology = efo.get_nx_ontology()
doid_ontology = doid.get_nx_ontology()

efo_graph = efo_ontology.graph
doid_graph = doid_ontology.graph

"""
biopartal_mapping_path = '/home/dhimmels/Documents/serg/data-mapping/bioportal/efo_doid/130914/mappings.rdf'
bioportal_efo_doid_pairs = mapping.bioportal.read_mapping(biopartal_mapping_path)
bioportal_pairs_rows = list()
for efo_id, doid_id in bioportal_efo_doid_pairs:
    doid_name = doid_graph.node[doid_id]['name'] if doid_id in doid_graph else None
    efo_name = efo_graph.node[efo_id]['name'] if efo_id in efo_graph else None
    row = {'efo_id': efo_id, 'doid_id': doid_id,
           'efo_name': efo_name, 'doid_name': doid_name}
    bioportal_pairs_rows.append(row)
bioportal_pairs_rows.sort()
bioportal_pairs_path = '/home/dhimmels/Documents/serg/data-mapping/bioportal/efo_doid/130914/mapping-pairs.txt'
with open(bioportal_pairs_path, 'w') as write_file:
    fieldnames = ['efo_id', 'doid_id', 'efo_name', 'doid_name']
    writer = csv.DictWriter(write_file, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(bioportal_pairs_rows)

bioportal_efo_to_doids = dict()
bioportal_doid_to_efos = dict()
for efo_id, doid_id in bioportal_efo_doid_pairs:
    bioportal_efo_to_doids.setdefault(efo_id, set()).add(doid_id)
    bioportal_doid_to_efos.setdefault(doid_id, set()).add(efo_id)

 
gwas_catalog = bioparser.data.Data().gwas_catalog
gwas_catalog.read_ebi_mappings()
efo_id_to_genes = gwas_catalog.get_efo_id_to_genes() # not using term cutoff
gwas_efo_ids = set(efo_id_to_genes)
efo_diseases = efo.get_diseases()
gwas_efo_disease_ids = gwas_efo_ids & efo_diseases

mapped_gwas_efos = gwas_efo_disease_ids & set(bioportal_efo_to_doids)
unmapped_gwas_efos = gwas_efo_disease_ids - set(bioportal_efo_to_doids)

print len(mapped_gwas_efos), 'GWAS EFO terms mapped by bioportal'
print len(unmapped_gwas_efos), 'GWAS EFO terms not mapped by bioportal'


gwas_pairs_rows = list()
for efo_id in mapped_gwas_efos:
    efo_name = efo_graph.node[efo_id]['name'] if efo_id in efo_graph else None
    doid_ids = bioportal_efo_to_doids[efo_id]
    if len(doid_ids) > 1:
        print efo_id, doid_ids
        #continue
    for doid_id in doid_ids:
        doid_name = doid_graph.node[doid_id]['name'] if doid_id in doid_graph else None
        row = {'efo_id': efo_id, 'doid_id': doid_id,
               'efo_name': efo_name, 'doid_name': doid_name,
               'method': 'bioportal'}
        gwas_pairs_rows.append(row)
gwas_pairs_rows.sort()
gwas_pairs_path = '/home/dhimmels/Documents/serg/data-mapping/manual/efo_doid/gwas-pairs-uneditted.txt'
with open(gwas_pairs_path, 'w') as write_file:
    fieldnames = ['efo_id', 'doid_id', 'efo_name', 'doid_name', 'method']
    writer = csv.DictWriter(write_file, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(gwas_pairs_rows)
"""
gwas_pairs_editted_path = '/home/dhimmels/Documents/serg/data-mapping/manual/efo_doid/gwas-pairs-editted.tsv'
with open(gwas_pairs_editted_path) as read_file:
    reader = csv.DictReader(read_file, delimiter='\t')
    manual_rows = list(reader)
for row in manual_rows:
    doid_id = row['doid_id']
    if doid_id not in doid_graph:
        print 'not found', doid_id

"""
for efo_id in unmapped_gwas_efos:
    efo_name = efo_graph.node[efo_id]['name'] if efo_id in efo_graph else None
    print efo_id + '\t' + efo_name


"""


