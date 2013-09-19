import csv
import os

import bioparser.data
from bioparser.metathesaurus import Concept
import mapping.bioportal

ictnet_dir = '/home/dhimmels/Documents/serg/ictnet/ictnet-creation2/'


"""
Some chemical-gene ctd interactions don't have an organism? All option?
Do we want approved symbol included in tb_gene_alias?
Are we excluding non-protein coding genes?
Several doid xrefs to omim use IDs not in the UMLS. Can table upload ignore those rows.
"""
class TableWriter(object):
    
    def __init__(self, table_name, fieldnames):
        """Class to write tab delimited text files encoding ictnet tables."""
        self.table_name = table_name
        self.fieldnames = fieldnames
        self.path = os.path.join(ictnet_dir, 'tables', table_name + '.txt')
        self.file = open(self.path, 'w')
        self.writer = csv.DictWriter(self.file, delimiter='\t',
            fieldnames=fieldnames, extrasaction='ignore')
        self.writer.writeheader()
    
    def writerow(self, row):
        self.writer.writerow(row)
    
    def writerows(self, rows, sort=True):
        if sort:
            rows = list(rows)
            rows.sort(key=lambda row: [row[key] for key in self.fieldnames])
        for row in rows:
            self.writerow(row)
        
    def close(self):
        self.file.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, *args, **kwargs):
        self.close()

def write_table(name, fielnames, rows):
    with TableWriter(name, fielnames) as writer:
        writer.writerows(rows)
    print 'Writing {} is complete'.format(name)


def get_row_tuple_to_rows(rows, primary_keys):
    """rows is a list of dictionaries. rows are considered unique based on
    their primary_keys"""
    row_tuple_to_rows = dict()
    for row in rows:
        row_tuple = tuple(row[key] for key in primary_keys)
        row_tuple_to_rows.setdefault(row_tuple, list()).append(row)
    return row_tuple_to_rows

def remove_duplicate_rows(rows, primary_keys):
    """rows is a list of dictionaries. rows are considered unique based on
    their primary_keys"""
    row_tuple_to_rows = get_row_tuple_to_rows(rows, primary_keys)
    return [rows[0] for rows in row_tuple_to_rows.itervalues()]

def condense_rows(rows, primary_keys, combine_keys, sep=', '):
    """Returns rows that are unique accross primary_keys. The values for 
    combine_keys are joined with the string sep."""
    condensed_rows = list()
    row_tuple_to_rows = get_row_tuple_to_rows(rows, primary_keys)
    for rows in row_tuple_to_rows.itervalues():
        condensed_row = rows[0].copy()
        for combine_key in combine_keys:
            joined_value = sep.join([row[combine_key] for row in rows])
        condensed_row[combine_key] = joined_value
        condensed_rows.append(condensed_row)
    return condensed_rows

hgnc = bioparser.data.Data().hgnc
symbol_to_gene = hgnc.get_symbol_to_gene()

"""
################################################################################
############# genes
tb_gene_fieldnames = ('gene_id', 'symbol', 'name', 'chromosome', 'group')
tb_gene_rows = list()
tb_gene_alias_fieldnames = ('gene_id', 'alias')
tb_gene_alias_rows = list()
genes = hgnc.get_genes()
for gene in genes:
    gene_id = gene.id_
    row = {'gene_id':gene_id, 'symbol': gene.symbol,
           'name': gene.name, 'chromosome': gene.chromosome,
           'group': gene.locus_group}
    tb_gene_rows.append(row)
    for alias in gene.synonyms + gene.previous_symbols:
        row = {'gene_id': gene_id, 'alias': alias}
        tb_gene_alias_rows.append(row)

# write tb_gene_rows
write_table('tb_gene', tb_gene_fieldnames, tb_gene_rows)

# write tb_gene_alias_rows
write_table('tb_gene_alias', tb_gene_alias_fieldnames, tb_gene_alias_rows)

################################################################################
############# ppi
tb_ppi_rows = list()
tb_ppi_fieldnames = ('source', 'target', 'pubmed', 'method', 'interaction_type',
                     'edge_type', 'sources', 'complex')

ppitrim = bioparser.data.Data().ppitrim
complex_interactions, binary_interactions = ppitrim.all_interactions()
complex_interactions = list() # excluded complexes
for is_complex, ppis in enumerate([binary_interactions, complex_interactions]):
    for ppi in ppis:
        row = ppi.copy()
        row['source'] = ppi['source'].id_
        row['target'] = ppi['target'].id_
        row['complex'] = is_complex
        tb_ppi_rows.append(row)

# write tb_ppi_rows
write_table('tb_ppi', tb_ppi_fieldnames, tb_ppi_rows)

################################################################################
############# diseases - Disease Ontology

tb_doid_fieldnames = 'doid_id', 'name'
tb_doid_rows = list()

tb_doid_ontology_fieldnames = ['parent', 'child']
tb_doid_ontology_rows = list()

tb_doid_omim_map_fieldnames = ['doid_id', 'omim_id']
tb_doid_omim_map_rows = list()

tb_doid_mesh_map_fieldnames = ['doid_id', 'mesh_id']
tb_doid_mesh_map_rows = list()



do_graph = bioparser.data.Data().doid.get_graph()
for node, data in do_graph.nodes(data=True):
    doid_id = data['id_']
    row = {'doid_id': doid_id,
           'name': data['name']}
    tb_doid_rows.append(row)
    for parent, child in do_graph.out_edges(node):
        parent_id = do_graph.node[parent]['id_']
        child_id = do_graph.node[child]['id_']
        tb_doid_ontology_rows.append({'parent': parent_id, 'child': child_id})
    
    # OMIM mappings
    for omim_id in data['xref'].get('OMIM', list()):
        row = {'doid_id': doid_id, 'omim_id': omim_id}
        tb_doid_omim_map_rows.append(row)

    # MESH mappings
    for mesh_id in data['xref'].get('MSH', list()):
        row = {'doid_id': doid_id, 'mesh_id': mesh_id}
        tb_doid_mesh_map_rows.append(row)

write_table('tb_doid', tb_doid_fieldnames, tb_doid_rows)
write_table('tb_doid_ontology', tb_doid_ontology_fieldnames, tb_doid_ontology_rows)
write_table('tb_doid_omim_map', tb_doid_omim_map_fieldnames, tb_doid_omim_map_rows)
write_table('tb_doid_mesh_map', tb_doid_mesh_map_fieldnames, tb_doid_mesh_map_rows)


################################################################################
############# diseases - EFO

tb_efo_fieldnames = 'doid_id', 'name'
tb_efo_rows = list()

tb_efo_ontology_fieldnames = 'parent', 'child'
tb_efo_ontology_rows = list()

efo = bioparser.data.Data().efo
efo_graph = efo.get_graph()
efo_diseases = efo.get_diseases()
efo_graph = efo_graph.subgraph(efo_diseases)
for node, data in efo_graph.nodes(data=True):
    row = {'doid_id': node,
           'name': data['name']}
    tb_efo_rows.append(row)
    for parent, child in efo_graph.out_edges(node):
        tb_efo_ontology_rows.append({'parent': parent, 'child': child})

with TableWriter('tb_efo', tb_efo_fieldnames) as writer:
    writer.writerows(tb_efo_rows)
with TableWriter('tb_efo_ontology', tb_efo_ontology_fieldnames) as writer:
    writer.writerows(tb_efo_ontology_rows)

################################################################################
############# diseases - OMIM
tb_omim_fieldnames = 'omim_id', 'name', 'types'
tb_omim_rows = list()

metathesaurus = bioparser.data.Data().metathesaurus
with metathesaurus:
    omim_concepts = metathesaurus.shelves['sources']['OMIM']
    #types_shelve = metathesaurus.shelves['types']
    retain_types = ['Disease or Syndrome', 'Anatomical Abnormality',
                    'Neoplastic Process', 'Congenital Abnormality']
    #types_shelve['Disease or Syndrome']
    #types_shelve['Anatomical Abnormality']
    #types_shelve['Neoplastic Process']
    #all_types = set()
    concepts_shelve = metathesaurus.shelves['concepts']
    for cui in omim_concepts:
        concept = concepts_shelve[cui]
        row = {'omim_id': concept.source_to_code['OMIM'],
               'name': concept.name,
               'types': concept.symantic_types}
        tb_omim_rows.append(row)
        #all_types |= concept.symantic_types

write_table('tb_omim', tb_omim_fieldnames, tb_omim_rows)

#omim_ids = set(row['omim_id'] for row in tb_omim_rows)
#mapped_omim_ids = set(row['omim_id'] for row in tb_doid_omim_map_rows)
#print mapped_omim_ids - omim_ids

################################################################################
############# EFO DOID bioportal mappings

tb_efo_doid_rows = list()
tb_efo_doid_fieldnames = 'efo_id', 'doid_id'

efo_doid_mapping_path = '/home/dhimmels/Documents/serg/data-mapping/bioportal/efo_doid/130914/mappings.rdf'


for source, target in mapping.bioportal.read_mapping(efo_doid_mapping_path):
    target = int(target.rsplit(':')[-1])
    row = {'efo_id': source, 'doid_id': target}
    tb_efo_doid_rows.append(row)
    
write_table('tb_efo_doid_map', tb_efo_doid_map_fieldnames, tb_efo_doid_rows)


################################################################################
############# Drugbank
tb_drugbank_fieldnames = ('drugbank_id', 'name', 'cas_number', 'type', 'groups')
tb_drugbank_rows = list()

tb_drugbank_alias_fieldnames = ('drugbank_id', 'alias')
tb_drugbank_alias_rows = list()


drugbank = bioparser.data.Data().drugbank
drugbank.read()
for drug in drugbank.drugs:
    drug['int_id'] = int(drug['drugbank_id'][2:])
    row = drug.copy()
    row['drugbank_id'] = row['int_id']
    row['cas_number'] = row.get('cas_number')
    row['groups'] = ', '.join(row.get('groups', list()))
    tb_drugbank_rows.append(row)
    aliases = drug.get('synonyms', list()) + drug.get('brands', list())
    for alias in aliases:
        row = {'drugbank_id': drug['int_id'], 'alias': alias}
        tb_drugbank_alias_rows.append(row)

write_table('tb_drugbank', tb_drugbank_fieldnames, tb_drugbank_rows)
write_table('tb_drugbank_alias', tb_drugbank_alias_fieldnames, tb_drugbank_alias_rows)
"""


################################################################################
############# CTD
tb_ctd_fieldnames = ('mesh_id', 'name')
tb_ctd_rows = list()

tb_ctd_alias_fieldnames = ('mesh_id', 'alias')
tb_ctd_alias_rows = list()

tb_ctd_drugbank_map_fieldnames = ('mesh_id', 'drugbank_id')
tb_ctd_drugbank_map_rows = list()

ctd = bioparser.data.Data().ctd
"""
for chemical in ctd.read_chemicals():
    mesh_id = chemical['ChemicalID']
    row = {'mesh_id': mesh_id,
           'name': chemical['ChemicalName'],
           'definition': chemical['Definition']}
    tb_ctd_rows.append(row)

    for alias in chemical['Synonyms']:
        row = {'mesh_id': mesh_id, 'alias': alias}
        tb_ctd_alias_rows.append(row)
    
    for drugbank_id in chemical['DrugBankIDs']:
        drugbank_int_id = int(drugbank_id[2:])
        row = {'mesh_id': mesh_id, 'drugbank_id': drugbank_int_id}
        tb_ctd_drugbank_map_rows.append(row)
    
write_table('tb_ctd', tb_ctd_fieldnames, tb_ctd_rows)
write_table('tb_ctd_alias', tb_ctd_alias_fieldnames, tb_ctd_alias_rows)
write_table('tb_ctd_drugbank_map', tb_ctd_drugbank_map_fieldnames, tb_ctd_drugbank_map_rows)



tb_medic_fieldnames = ('disease_id', 'name')
tb_medic_rows = list()

for disease in ctd.read_diseases():
    row = {'disease_id': disease['DiseaseID'],
           'name': disease['DiseaseName']}
    tb_medic_rows.append(row)

write_table('tb_medic', tb_medic_fieldnames, tb_medic_rows)


tb_ctd_medic_therapy_fieldnames = ('mesh_id', 'medic_id', 'pubmeds')
tb_ctd_medic_therapy_rows = list()

for therapy in ctd.read_chemical2diseases():
    if 'therapeutic' not in therapy['DirectEvidence']:
        continue
    pubmeds = ', '.join(therapy['PubMedIDs'])
    row = {'mesh_id': therapy['ChemicalID'],
           'medic_id': therapy['DiseaseID'],
           'pubmeds': pubmeds}
    tb_ctd_medic_therapy_rows.append(row)

write_table('tb_ctd_medic_therapy', tb_ctd_medic_therapy_fieldnames, tb_ctd_medic_therapy_rows)
"""

tb_ctd_gene_ixn_fieldnames = ('mesh_id', 'symbol', 'organism_id', 'pubmeds')
tb_ctd_gene_ixn_rows = list()

id_to_organism = dict()
for ixn in ctd.read_chemical2genes():
    mesh_id = ixn['ChemicalID']
    symbol = ixn['GeneSymbol']
    gene = symbol_to_gene.get(symbol)
    if not gene:
        continue
    pubmeds = ', '.join(ixn['PubMedIDs'])
    organism_id = ixn['OrganismID']
    organism = ixn['Organism']
    
    # Some chemical-gene ctd interactions don't have an organism 
    if not organism_id:
        organism_id = 0
        organism = 'Unkown'
    id_to_organism[organism_id] = organism
    #row_tuple = mesh_id, gene.symbol, pubmeds, organism_id
    row = {'mesh_id': mesh_id, 'symbol': gene.symbol, 'pubmeds': pubmeds,
           'organism_id': organism_id}
    tb_ctd_gene_ixn_rows.append(row)

tb_ctd_gene_ixn_rows = remove_duplicate_rows(tb_ctd_gene_ixn_rows, tb_ctd_gene_ixn_fieldnames)
tb_ctd_gene_ixn_rows = condense_rows(tb_ctd_gene_ixn_rows, ['mesh_id', 'symbol', 'organism_id'], ['pubmeds'])


write_table('tb_ctd_gene_ixn', tb_ctd_gene_ixn_fieldnames, tb_ctd_gene_ixn_rows)

tb_organism_fieldnames = 'organism_id', 'name'
tb_organism_rows = [{'organism_id': int(id_), 'name': name}
                    for id_, name in id_to_organism.iteritems()]
write_table('tb_organism', tb_organism_fieldnames, tb_organism_rows)



