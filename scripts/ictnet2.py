import csv
import os

import bioparser.data
from bioparser.metathesaurus import Concept
import mapping.bioportal

ictnet_dir = '/home/dhimmels/Documents/serg/ictnet/ictnet-creation2/'
input_dir = os.path.join(ictnet_dir, 'input')


class Table(object):
    
    def __init__(self, name, fieldnames):
        self.name = name
        self.file_name = 'tb_{}.txt'.format(name)
        self.path = os.path.join(ictnet_dir, 'tables', self.file_name)
        self.rows = list()
        self.fieldnames = fieldnames
    
    def sort_rows(self):
        self.rows.sort(key=lambda row: [row[key] for key in self.fieldnames])

    def write(self):
        self.sort_rows()
        with open(self.path, 'w') as write_file:
            writer = csv.DictWriter(write_file, delimiter='\t',
                fieldnames=self.fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(self.rows)
        print 'Writing {} is complete'.format(self.file_name)

    def append(self, row):
        self.rows.append(row)

    def get_row_tuple_to_rows(self, primary_keys):
        """rows is a list of dictionaries. rows are considered unique based on
        their primary_keys"""
        row_tuple_to_rows = dict()
        for row in self.rows:
            row_tuple = tuple(row[key] for key in primary_keys)
            row_tuple_to_rows.setdefault(row_tuple, list()).append(row)
        return row_tuple_to_rows

    def remove_duplicates(self, primary_keys):
        """Remove duplicate rows"""
        row_tuple_to_rows = self.get_row_tuple_to_rows(primary_keys)
        self.rows = [rows[0] for rows in row_tuple_to_rows.itervalues()]

    def condense(self, primary_keys, combine_keys, sep=', '):
        """Returns rows that are unique accross primary_keys. The values for 
        combine_keys are joined with the string sep."""
        condensed_rows = list()
        row_tuple_to_rows = self.get_row_tuple_to_rows(primary_keys)
        for rows in row_tuple_to_rows.itervalues():
            condensed_row = rows[0].copy()
            for combine_key in combine_keys:
                joined_value = sep.join([row[combine_key] for row in rows])
            condensed_row[combine_key] = joined_value
            condensed_rows.append(condensed_row)
        self.rows = condensed_rows

"""
Some chemical-gene ctd interactions don't have an organism? All option?
Do we want approved symbol included in tb_gene_alias?
Are we excluding non-protein coding genes?
Several doid xrefs to omim use IDs not in the UMLS. Can table upload ignore those rows.
"""


hgnc = bioparser.data.Data().hgnc
symbol_to_gene = hgnc.get_symbol_to_gene()

"""
################################################################################
############# HGNC - Genes
tb_gene = Table('gene', ['hgnc_id', 'symbol', 'name', 'chromosome', 'group'])
tb_gene_alias = Table('gene_alias', ['hgnc_id', 'alias'])

genes = hgnc.get_genes()
for gene in genes:
    gene_id = gene.id_
    row = {'hgnc_id':gene_id, 'symbol': gene.symbol,
           'name': gene.name, 'chromosome': gene.chromosome,
           'group': gene.locus_group}
    tb_gene.append(row)
    for alias in gene.synonyms + gene.previous_symbols:
        row = {'hgnc_id': gene_id, 'alias': alias}
        tb_gene_alias.append(row)

tb_gene.write()
tb_gene_alias.write()

################################################################################
############# ppi
tb_ppi = Table('ppi', ['source', 'target', 'pubmed', 'method', 'interaction_type',
              'edge_type', 'sources', 'complex'])

ppitrim = bioparser.data.Data().ppitrim
complex_interactions, binary_interactions = ppitrim.all_interactions()
complex_interactions = list() # excluded complexes
for is_complex, ppis in enumerate([binary_interactions, complex_interactions]):
    for ppi in ppis:
        row = ppi.copy()
        row['source'] = ppi['source'].id_
        row['target'] = ppi['target'].id_
        row['complex'] = is_complex
        tb_ppi.append(row)
tb_ppi.write()

################################################################################
############# diseases - Disease Ontology

tb_doid = Table('doid', ['doid_id', 'name'])
tb_doid_ontology = Table('doid_ontology', ['parent', 'child'])
tb_doid_omim_map = Table('doid_omim_map', ['doid_id', 'omim_id'])
tb_doid_medic_map = Table('doid_medic_map', ['doid_id', 'medic_id'])


do_graph = bioparser.data.Data().doid.get_graph()
for node, data in do_graph.nodes(data=True):
    doid_id = data['id_']
    row = {'doid_id': doid_id,
           'name': data['name']}
    tb_doid.append(row)
    for parent, child in do_graph.out_edges(node):
        parent_id = do_graph.node[parent]['id_']
        child_id = do_graph.node[child]['id_']
        tb_doid_ontology.append({'parent': parent_id, 'child': child_id})
    
    # OMIM mappings
    for omim_id in data['xref'].get('OMIM', list()):
        row = {'doid_id': doid_id, 'omim_id': omim_id}
        tb_doid_omim_map.append(row)
        row = {'doid_id': doid_id, 'medic_id': 'OMIM:' + omim_id}
        tb_doid_medic_map.append(row)

    # MESH mappings
    for mesh_id in data['xref'].get('MSH', list()):
        row = {'doid_id': doid_id, 'medic_id': 'MESH:' + mesh_id}
        tb_doid_medic_map.append(row)

tb_doid.write()
tb_doid_ontology.write()
tb_doid_omim_map.write()
tb_doid_medic_map.write()

################################################################################
############# diseases - EFO

tb_efo = Table('efo', ['doid_id', 'name'])
tb_efo_ontology = Table('efo_ontology', ['parent', 'child'])

efo = bioparser.data.Data().efo
efo_graph = efo.get_graph()
efo_diseases = efo.get_diseases()
efo_graph = efo_graph.subgraph(efo_diseases)
for node, data in efo_graph.nodes(data=True):
    row = {'doid_id': node,
           'name': data['name']}
    tb_efo.append(row)
    for parent, child in efo_graph.out_edges(node):
        tb_efo_ontology.append({'parent': parent, 'child': child})

tb_efo.write()
tb_efo_ontology.write()


tb_gene_efo_gwas = Table('gene_efo_gwas',
    ['hgnc_id', 'efo_id', 'p-value', 'OR or beta', 'SNPs', 'pubmed'])

gwas_catalog = bioparser.data.Data().gwas_catalog
#path = '/home/dhimmels/Documents/serg/data-sources/gwas-catalog/GWAS-EFO-Mappings092012.txt'
gwas_catalog.read_ebi_mappings()
gwas_catalog.apply_fdr()
for gwas_row in gwas_catalog.get_rows():
    efo_id = gwas_row.get('efo_id')
    if not efo_id:
        continue
    genes = set(symbol_to_gene.get(symbol) for symbol in gwas_row['genes'])
    genes.discard(None)    
    for gene in genes:
        row = {'efo_id': efo_id, 'hgnc_id': gene.id_, 
               'p-value': gwas_row['p-Value'], 'OR or beta': gwas_row['OR or beta'],
               'SNPs': gwas_row['SNPs'], 'pubmed': gwas_row['PUBMEDID']}
        tb_gene_efo_gwas.append(row)
tb_gene_efo_gwas.write()

################################################################################
############# diseases - OMIM

tb_gene_omim_morbidmap = Table('gene_omim_morbidmap', ['hgnc_id', 'omim_id'])
omim_ids = set()
morbid_map = bioparser.data.Data().morbid_map
id_to_gene_tuples = morbid_map.get_id_to_gene_tuples()
for omim_id, gene in id_to_gene_tuples:
    row = {'hgnc_id': gene.id_, 'omim_id': omim_id}
    tb_gene_omim_morbidmap.append(row)
    omim_ids.add(omim_id)
tb_gene_omim_morbidmap.write()


tb_omim = Table('omim', ['omim_id', 'name'])
# Use the UMLS metathesaurus to get the name for the OMIM ID
metathesaurus = bioparser.data.Data().metathesaurus
with metathesaurus:
    omim_to_concept = metathesaurus.get_source_code_to_concept('OMIM')
    for omim_id in omim_ids:
        concept = omim_to_concept.get(omim_id)
        omim_name = concept.source_to_name['OMIM'] if concept else None
        row = {'omim_id': omim_id, 'name': omim_name}
        tb_omim.append(row)
tb_omim.write()

################################################################################
############# EFO DOID bioportal mappings
tb_efo_doid_map = Table('efo_doid_map', ['efo_id', 'doid_id'])
efo_doid_mapping_path = '/home/dhimmels/Documents/serg/data-mapping/bioportal/efo_doid/130914/mappings.rdf'
for source, target in mapping.bioportal.read_mapping(efo_doid_mapping_path):
    target = int(target.rsplit(':')[-1])
    row = {'efo_id': source, 'doid_id': target}
    tb_efo_doid_map.append(row)
tb_efo_doid_map.write()

################################################################################
############# Drugbank
tb_drugbank = Table('drugbank', ['drugbank_id', 'name', 'cas_number', 'type', 'groups'])
tb_drugbank_alias = Table('drugbank_alias', ['drugbank_id', 'alias'])

drugank_id_to_int = lambda s: int(s[2: ])
drugbank = bioparser.data.Data().drugbank
drugbank.read()
for drug in drugbank.drugs:
    drug['int_id'] = drugank_id_to_int(drug['drugbank_id'])
    row = drug.copy()
    row['drugbank_id'] = row['int_id']
    row['cas_number'] = row.get('cas_number')
    row['groups'] = ', '.join(row.get('groups', list()))
    tb_drugbank.append(row)
    aliases = drug.get('synonyms', list()) + drug.get('brands', list())
    for alias in aliases:
        row = {'drugbank_id': drug['int_id'], 'alias': alias}
        tb_drugbank_alias.append(row)
tb_drugbank.write()
tb_drugbank_alias.write()


tb_drug_gene_drugbank = Table('drug_gene_drugbank',
    ['drugbank_id', 'hgnc_id', 'pharmacological', 'actions'])

id_to_partner = drugbank.get_id_to_partner()
for target in drugbank.targets:
    partner_id = target['partner']
    partner = id_to_partner[partner_id]
    if partner['species'] != 'Homo sapiens':
        continue
    gene = symbol_to_gene.get(partner.get('gene_name'))
    if not gene:
        continue
    row = {'drugbank_id': drugank_id_to_int(target['drugbank_id']),
           'hgnc_id': gene.id_,
           'pharmacological': int(target['known_action'] == 'yes'),
           'actions': ', '.join(target.get('actions', []))}
    tb_drug_gene_drugbank.append(row)
tb_drug_gene_drugbank.write()

################################################################################
############# CTD
tb_ctd = Table('ctd', ['mesh_id', 'name'])
tb_ctd_alias = Table('ctd_alias', ['mesh_id', 'alias'])
tb_ctd_drugbank_map = Table('ctd_drugbank_map', ['mesh_id', 'drugbank_id'])

ctd = bioparser.data.Data().ctd
for chemical in ctd.read_chemicals():
    mesh_id = chemical['ChemicalID']
    row = {'mesh_id': mesh_id,
           'name': chemical['ChemicalName'],
           'definition': chemical['Definition']}
    tb_ctd.append(row)

    for alias in chemical['Synonyms']:
        row = {'mesh_id': mesh_id, 'alias': alias}
        tb_ctd_alias.append(row)
    
    for drugbank_id in chemical['DrugBankIDs']:
        drugbank_int_id = int(drugbank_id[2:])
        row = {'mesh_id': mesh_id, 'drugbank_id': drugbank_int_id}
        tb_ctd_drugbank_map.append(row)
    
tb_ctd.write()
tb_ctd_alias.write()
tb_ctd_drugbank_map.write()

tb_medic = Table('medic', ['medic_id', 'name'])
for disease in ctd.read_diseases():
    row = {'medic_id': disease['DiseaseID'],
           'name': disease['DiseaseName']}
    tb_medic.append(row)
tb_medic.write()

#medic_ids = set(row['medic_id'] for row in tb_medic_rows)
#tb_doid_medic_map_rows = list(row for row in tb_doid_medic_map_rows
#                              if row['medic_id'] in medic_ids)
#tb_doid_medic_map.write()

tb_ctd_medic_therapy = Table('ctd_medic_therapy', ['mesh_id', 'medic_id', 'pubmeds'])

for therapy in ctd.read_chemical2diseases():
    if 'therapeutic' not in therapy['DirectEvidence']:
        continue
    pubmeds = ', '.join(therapy['PubMedIDs'])
    row = {'mesh_id': therapy['ChemicalID'],
           'medic_id': therapy['DiseaseID'],
           'pubmeds': pubmeds}
    tb_ctd_medic_therapy.append(row)
tb_ctd_medic_therapy.write()

tb_ctd_gene_ixn = Table('ctd_gene_ixn',
    ['mesh_id', 'hgnc_id', 'organism_id', 'pubmeds'])
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
    row = {'mesh_id': mesh_id, 'hgnc_id': gene.id_, 'pubmeds': pubmeds,
           'organism_id': organism_id}
    tb_ctd_gene_ixn.append(row)

tb_ctd_gene_ixn.remove_duplicates(tb_ctd_gene_ixn.fieldnames)
tb_ctd_gene_ixn.condense(['mesh_id', 'hgnc_id', 'organism_id'], ['pubmeds'])
tb_ctd_gene_ixn.write()

tb_organism = Table('organism', ['organism_id', 'name'])
for id_, name in id_to_organism.iteritems():
    tb_organism.append({'organism_id': int(id_), 'name': name})
tb_organism.write()


tb_gene_medic_ctd = Table('gene_medic_ctd', ['hgnc_id', 'medic_id', 'pubmeds'])
for ctd_row in ctd.read_gene2disease():
    if 'marker/mechanism' not in ctd_row['DirectEvidence']:
        continue
    symbol = ctd_row['GeneSymbol']
    gene = symbol_to_gene.get(symbol)
    if not gene:
        continue
    pubmeds = ', '.join(ctd_row['PubMedIDs'])
    omims = ', '.join(ctd_row['OmimIDs'])
    row = {'hgnc_id': gene.id_, 'medic_id': ctd_row['DiseaseID'], 'pubmeds': pubmeds}
    tb_gene_medic_ctd.append(row)
tb_gene_medic_ctd.write()

###############################################################################
### MicroRNA
mircat = bioparser.data.Data().mircat

tb_mircat = Table('mircat', ['mircat_name', 'sequence', 'location'])
tb_mircat.rows = list(mircat.read_mirna())
tb_mircat.write()

tb_gene_mirna = Table('gene_mirna', ['hgnc_id','mircat_name', 'pubmed'])
for row in mircat.read_targets():
    row['hgnc_id'] = row['gene'].id_
    tb_gene_mirna.append(row)
# Condensation isn't needed because each interaction has only a single pmid.
# This seems unlikely and could indicate an upstream bug.
#tb_gene_mirna.condense(['hgnc_id','mircat_name'], ['pubmed'])
tb_gene_mirna.write()

"""

def read_input(name):
    path = os.path.join(ictnet_dir, 'input', name + '.txt')
    with open(path) as read_file:
        reader = csv.DictReader(read_file, delimiter='\t')
        rows = list(reader)
    return rows

tissues = read_input('gep-tissues')
tissue_to_id = {row['gep_tissue']: i for i, row in enumerate(tissues)}
tb_tissue = Table('tissue', ['tissue_id', 'name'])
for tissue, tissue_id in tissue_to_id.items():
    row = {'tissue_id': tissue_id, 'name': tissue}
    tb_tissue.append(row)
tb_tissue.write()

doid = bioparser.data.Data().doid
do_graph = doid.get_graph()

for node, data in do_graph.nodes_iter(data=True):
    data['tissues'] = set()

tissue_doid_pairs = read_input('tissue-to-doid')
for pair in tissue_doid_pairs:
    disease_id = pair['doid']
    tissue = pair['gep_tissue']
    diseases = doid.get_descendents(disease_id)
    diseases.add(disease_id)
    for disease in diseases:
        do_graph.node[disease]['tissues'].add(tissue)

tb_doid_tissue = Table('doid_tissue', ['doid_id', 'tissue_id'])
for node, data in do_graph.nodes_iter(data=True):
    doid_id = data['id_']
    for tissue in data['tissues']:
        row = {'doid_id': doid_id, 'tissue_id': tissue_to_id[tissue]}
        tb_doid_tissue.append(row)
tb_doid_tissue.write()
"""
# Read tissues
gep_tissues_path = os.path.join(input_dir, 'gep-tissues.txt')
with open(gep_tissues_path) as read_file:
    tissues = [row.rstrip() for row in read_file]
tissue_to_id = {tissue: i for i, tissue in enumerate(tissues)}
"""




