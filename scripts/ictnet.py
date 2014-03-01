import datetime
import csv
import collections
import os
import pprint

import bioparser.data
from bioparser.metathesaurus import Concept

ictnet_dir = '/home/dhimmels/Documents/serg/ictnet/ictnet-creation/'
input_dir = os.path.join(ictnet_dir, 'input')

tables = list()

class Table(object):
    
    def __init__(self, name, fieldnames):
        self.name = name
        self.file_name = 'tb_{}.txt'.format(name)
        self.path = os.path.join(ictnet_dir, 'tables', self.file_name)
        self.rows = list()
        self.fieldnames = fieldnames
        tables.append(self)
    
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

    def __iter__(self):
        return iter(self.rows)

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

    def schema(self):
        field_to_type = collections.OrderedDict()
        for fieldname in self.fieldnames:
            values = [row[fieldname] for row in self.rows]
            # strings
            if isinstance(values[0], str):
                str_lens = {0 if s is None else len(s) for s in values}
                try:
                    length, = str_lens
                    sql_type = 'CHAR({})'.format(length)
                except ValueError:
                    length = max(str_lens)
                    sql_type = 'VARCHAR({})'.format(length)

            # integers
            elif isinstance(values[0], int):
                sql_type = 'INT({}, {})'.format(min(values), max(values))

            # floats
            elif isinstance(values[0], float):
                sql_type = 'DECIMAL({}, {})'.format(min(values), max(values))

            else:
                print type(values[0])
                raise ValueError('unknown type for mySQL conversion')

            field_to_type[fieldname] = sql_type
        return field_to_type

    def remove_invalid_references(self, other, fieldname):
        valid_refs = {row[fieldname] for row in other}
        self.rows = [row for row in self if row[fieldname] in valid_refs]


def read_input(file_name):
    path = os.path.join(ictnet_dir, 'input', file_name)
    with open(path) as read_file:
        reader = csv.DictReader(read_file, delimiter='\t')
        rows = list(reader)
    return rows

def factor_table(table, fieldname):
    """
    Create an int index for a table field. Create a new table where each row
    contains an int id and the corresponding value.
    """
    id_fieldname = '{}_id'.format(fieldname)
    name_fieldname = 'name'#'{}_name'.format(fieldname)
    table_name = '{}_{}'.format(table.name, fieldname)
    factored_table = Table(table_name, [id_fieldname, name_fieldname])
    
    value_to_id = dict()
    values = {row[fieldname] for row in table.rows}
    for i, value in enumerate(values):
        value_to_id[value] = i
        row = {id_fieldname: i, name_fieldname: value}
        factored_table.append(row)
        
    for row in table.rows:
        row[id_fieldname] = value_to_id[row[fieldname]]
    return factored_table

tb_resource_version = Table('resource_versions', ['resource', 'version'])
def add_version(resource, version):
    tb_resource_version.append({'resource': resource, 'version': version})

add_version('date', str(datetime.date.today()))

hgnc = bioparser.data.Data().hgnc
symbol_to_gene = hgnc.get_symbol_to_gene()
add_version('hgnc', hgnc.directory)
################################################################################
############# HGNC - Genes

tb_gene = Table('gene', ['gene_id', 'hgnc_id', 'symbol', 'name', 'location', 'group_id', 'type_id'])
tb_gene_alias = Table('gene_alias', ['gene_id', 'alias'])

genes = hgnc.get_genes()
for gene in genes:
    gene_id = gene.int_id
    row = {'gene_id':gene_id, 'hgnc_id': gene.hgnc_id, 'symbol': gene.symbol,
           'name': gene.name, 'location': gene.chromosome,
           'group': gene.locus_group, 'type': gene.locus_type}
    tb_gene.append(row)
    for alias in gene.synonyms + gene.previous_symbols:
        row = {'gene_id': gene_id, 'alias': alias}
        tb_gene_alias.append(row)

tb_gene_group = factor_table(tb_gene, 'group')
tb_gene_type = factor_table(tb_gene, 'type')
tb_gene_group.write()
tb_gene_type.write()

tb_gene.write()
tb_gene_alias.write()

################################################################################
############# ppi
tb_ppi = Table('ppi', ['source', 'target', 'pubmed', 'method_id', 'interaction_type_id',
              'edge_type_id', 'sources', 'complex'])

ppitrim = bioparser.data.Data().ppitrim
add_version('ppitrim', ppitrim.consolidated_path)
complex_interactions, binary_interactions = ppitrim.all_interactions()
complex_interactions = list() # excluded complexes
for is_complex, ppis in enumerate([binary_interactions, complex_interactions]):
    for ppi in ppis:
        row = ppi.copy()
        row['source'] = ppi['source'].int_id
        row['target'] = ppi['target'].int_id
        row['complex'] = is_complex
        tb_ppi.append(row)

tb_ppi_method = factor_table(tb_ppi, 'method')
tb_ppi_interaction_type = factor_table(tb_ppi, 'interaction_type')
tb_ppi_edge_type = factor_table(tb_ppi, 'edge_type')

tb_ppi_method.write()
tb_ppi_interaction_type.write()
tb_ppi_edge_type.write()

tb_ppi.write()


################################################################################
############# diseases - Disease Ontology

tb_doid = Table('doid', ['doid_id', 'doid_code', 'name'])
tb_doid_ontology = Table('doid_ontology', ['parent', 'child'])
tb_doid_omim_map = Table('doid_omim_map', ['doid_id', 'omim_id'])
tb_doid_medic_map = Table('doid_medic_map', ['doid_id', 'medic_id'])
tb_doid_efo_map = Table('doid_efo_map', ['doid_id', 'efo_id'])

doid = bioparser.data.Data().doid
do_graph = doid.get_graph()
add_version('doid', doid.directory)

for node, data in do_graph.nodes(data=True):
    doid_id = data['int_id']
    row = {'doid_id': doid_id, 'doid_code': node,
           'name': data['name']}
    tb_doid.append(row)
    for parent, child in do_graph.out_edges(node):
        parent_id = do_graph.node[parent]['int_id']
        child_id = do_graph.node[child]['int_id']
        tb_doid_ontology.append({'parent': parent_id, 'child': child_id})
    
    xref = data.get('xref', dict())
    
    # OMIM mappings
    for omim_id in xref.get('OMIM', list()):
        row = {'doid_id': doid_id, 'omim_id': omim_id}
        tb_doid_omim_map.append(row)
        row = {'doid_id': doid_id, 'medic_id': 'OMIM:' + omim_id}
        tb_doid_medic_map.append(row)

    # MESH mappings
    for mesh_id in xref.get('MSH', list()):
        row = {'doid_id': doid_id, 'medic_id': 'MESH:' + mesh_id}
        tb_doid_medic_map.append(row)

    # EFO mappings
    for efo_id in xref.get('EFO', list()):
        row = {'doid_id': doid_id, 'efo_id': 'EFO:' + efo_id}
        tb_doid_efo_map.append(row)
    for pat_id in xref.get('EFOpat_id', list()):
        row = {'doid_id': doid_id, 'efo_id': 'pat_id:' + pat_id}
        tb_doid_efo_map.append(row)

tb_doid.write()
tb_doid_ontology.write()

################################################################################
############# diseases - EFO

tb_efo = Table('efo', ['efo_id', 'name'])
tb_efo_ontology = Table('efo_ontology', ['parent', 'child'])

efo = bioparser.data.Data().efo
add_version('efo', efo.directory)
efo_diseases = efo.get_diseases()
efo_graph = efo.get_graph().subgraph(efo_diseases)
for node, data in efo_graph.nodes(data=True):
    row = {'efo_id': node,
           'name': data['name']}
    tb_efo.append(row)
    for parent, child in efo_graph.out_edges(node):
        tb_efo_ontology.append({'parent': parent, 'child': child})

tb_efo.write()
tb_efo_ontology.write()


tb_gene_efo_gwas = Table('gene_efo_gwas',
    ['gene_id', 'efo_id', 'p-value', 'OR or beta', 'SNPs', 'pubmed'])

gwas_catalog = bioparser.data.Data().gwas_catalog
add_version('gwas_catalog', gwas_catalog.gwas_dir)

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
        row = {'efo_id': efo_id, 'gene_id': gene.int_id, 
               'p-value': gwas_row['p-Value'], 'OR or beta': gwas_row['OR or beta'],
               'SNPs': gwas_row['SNPs'], 'pubmed': gwas_row['PUBMEDID']}
        tb_gene_efo_gwas.append(row)

################################################################################
############# diseases - OMIM

morbid_map = bioparser.data.Data().morbid_map
add_version('morbid_map', morbid_map.omim_dir)


tb_omim = Table('omim', ['omim_id', 'name', 'type_id'])
for disorder in morbid_map.get_disorders():
    row = {'omim_id': disorder['mim_number'],
           'name': disorder['disorder_name'],
           'type': disorder['disorder_type']}
    tb_omim.append(row)
tb_omim_type = factor_table(tb_omim, 'type')
tb_omim.write()
tb_omim_type.write()


tb_gene_omim_morbidmap_fieldnames = ['gene_id', 'omim_id', 'association_type', 'confirmed']
tb_gene_omim_morbidmap = Table('gene_omim_morbidmap', tb_gene_omim_morbidmap_fieldnames)
for association in morbid_map.get_associations():
    row = association.copy()
    row['omim_id'] = row.pop('mim_number')
    row['gene_id'] = row.pop('gene').int_id
    tb_gene_omim_morbidmap.append(row)
tb_gene_omim_morbidmap.write()



################################################################################
############# Drugbank
tb_drugbank = Table('drugbank', ['drugbank_id', 'name', 'cas_number', 'type', 'groups'])
tb_drugbank_alias = Table('drugbank_alias', ['drugbank_id', 'alias'])

drugank_id_to_int = lambda s: int(s[2: ])
drugbank = bioparser.data.Data().drugbank
add_version('drugbank', drugbank.drugbank_dir)
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
    ['drugbank_id', 'gene_id', 'pharmacological', 'actions'])

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
           'gene_id': gene.int_id,
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
add_version('ctd', ctd.directory)

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
    row = {'mesh_id': 'MESH:' + therapy['ChemicalID'],
           'medic_id': therapy['DiseaseID'],
           'pubmeds': pubmeds}
    tb_ctd_medic_therapy.append(row)
tb_ctd_medic_therapy.write()


tb_ctd_gene_ixn = Table('ctd_gene_ixn',
    ['mesh_id', 'gene_id', 'pubmeds'])
for ixn in ctd.read_chemical2genes():
    mesh_id = 'MESH:' + ixn['ChemicalID']
    symbol = ixn['GeneSymbol']
    gene = symbol_to_gene.get(symbol)
    if not gene:
        continue
    pubmeds = ', '.join(ixn['PubMedIDs'])
    organism_id = ixn['OrganismID']
    if organism_id != '9606':
        continue
    row = {'mesh_id': mesh_id, 'gene_id': gene.int_id, 'pubmeds': pubmeds}
    tb_ctd_gene_ixn.append(row)

tb_ctd_gene_ixn.remove_duplicates(tb_ctd_gene_ixn.fieldnames)
tb_ctd_gene_ixn.condense(['mesh_id', 'gene_id'], ['pubmeds'])
tb_ctd_gene_ixn.write()


tb_gene_medic_ctd = Table('gene_medic_ctd', ['gene_id', 'medic_id', 'pubmeds'])
for ctd_row in ctd.read_gene2disease_filtered():
    if 'marker/mechanism' not in ctd_row['DirectEvidence']:
        continue
    symbol = ctd_row['GeneSymbol']
    gene = symbol_to_gene.get(symbol)
    if not gene:
        continue
    pubmeds = ', '.join(ctd_row['PubMedIDs'])
    omims = ', '.join(ctd_row['OmimIDs'])
    row = {'gene_id': gene.int_id, 'medic_id': ctd_row['DiseaseID'], 'pubmeds': pubmeds}
    tb_gene_medic_ctd.append(row)
tb_gene_medic_ctd.write()

################################################################################
############# MSB PREDICT Indications
tb_drugbank_omim_predict = Table('drugbank_omim_predict', ['drugbank_id', 'omim_id'])
nature_predict = bioparser.data.Data().nature_predict
for row in nature_predict.omim_drugbank_mapper():
    row['drugbank_id'] = drugank_id_to_int(row['drugbank_id'])
    tb_drugbank_omim_predict.append(row)


###############################################################################
### MicroRNA
mircat = bioparser.data.Data().mircat
tb_mirna = Table('mirna', ['source_gene_id', 'target_gene_id', 'pubmed'])
for row in mircat.interaction_generator():
    row['source_gene_id'] = row['source'].int_id
    row['target_gene_id'] = row['target'].int_id
    tb_mirna.append(row)
# Condensation isn't needed because each interaction has only a single pmid.
# This seems unlikely and could indicate an upstream bug.
#tb_gene_mirna.condense(['gene_id','mircat_name'], ['pubmed'])
tb_mirna.write()

#######################
## BRENDA Tissue Ontology
tb_bto = Table('bto', ['bto_id', 'name'])
tb_bto_ontology = Table('bto_ontology', ['parent', 'child', 'relationship_id'])

bto = bioparser.data.Data().bto
add_version('bto', bto.directory)
bto_graph = bto.get_animal_tissue_subgraph()
for node, data in bto_graph.nodes(data=True):
    row = {'bto_id': node,
           'name': data['name']}
    tb_bto.append(row)
    for parent, child, key in bto_graph.out_edges(node, keys=True):
        tb_bto_ontology.append({'parent': parent, 'child': child, 
                                'relationship': key})

tb_bto_relationship = factor_table(tb_bto_ontology, 'relationship')

tb_bto.write()
tb_bto_ontology.write()
tb_bto_relationship.write()

#########################
## GNF Gene Atlas - Gene Expression by Tissue
tb_gnf = Table('gnf', ['gene_id', 'bto_id', 'expr', 'log10_expr'])
for expression in bioparser.data.Data().gnf.expression_generator():
    expression['gene_id'] = expression['gene'].int_id
    tb_gnf.append(expression)
tb_gnf.write()


#######################
## disease to tissue mappings
doid = bioparser.data.Data().doid
disease_ontology = doid.get_ontology()

disease_ontology.initialize_attribute('direct_tissue_mapping')
tissue_doid_pairs = read_input('Tissue Mappings - BTO-DOID.tsv')
for pair in tissue_doid_pairs:
    disease_id = pair['doid_code']
    tissue_id = pair['bto_id']
    disease_ontology.graph.node[disease_id]['direct_tissue_mapping'].add(tissue_id)
disease_ontology.most_specific_superior_annotations('direct_tissue_mapping', 'tissues')


tb_doid_bto = Table('doid_bto', ['doid_id', 'bto_id'])
for node, data in disease_ontology.graph.nodes_iter(data=True):
    doid_id = data['int_id']
    for bto_id in data['tissues']:
        row = {'doid_id': doid_id, 'bto_id': bto_id}
        tb_doid_bto.append(row)
tb_doid_bto.write()

#######################
## Side Effects - meddra
tb_side_effect = Table('side_effect', ['umls_id', 'umls_name'])
#tb_side_effect = Table('side_effect', ['umls_id', 'umls_name', 'meddra_code', 'meddra_name'])

meddra = bioparser.data.Data().meddra
meddra.annotate_umls(version='2011AB')
add_version('umls_for_meddra', '2011AB')
meddra_ontology = meddra.get_networkx_ontology()

for node, data in meddra_ontology.graph.nodes_iter(data=True):
    if not 'umls_id' in data:
        #print node, 'no umls match for meddra'
        continue
    row = {'meddra_code': node, 'umls_id': data['umls_id'], 
           'meddra_name': data['name'], 'umls_name': data['umls_name']}
    tb_side_effect.append(row)
tb_side_effect.write()

#######################
## Meddra Side Effect to Tissue
tb_side_effect_bto = Table('side_effect_bto', ['umls_id', 'bto_id'])

meddra_ontology.initialize_attribute('direct_tissue_mapping')
for tissue_mapping in read_input('Tissue Mappings - BTO-MedDRA.tsv'):
    bto_id = tissue_mapping['bto_id']
    meddra_code = tissue_mapping['meddra_code']
    meddra_ontology.graph.node[meddra_code]['direct_tissue_mapping'].add(bto_id)
meddra_ontology.most_specific_superior_annotations('direct_tissue_mapping', 'tissues')

for node, data in meddra_ontology.graph.nodes_iter(data=True):
    umls_id = data.get('umls_id')
    if not umls_id:
        continue
    for bto_id in data['tissues']:
        row = {'umls_id': umls_id, 'bto_id': bto_id}
        tb_side_effect_bto.append(row)
tb_side_effect_bto.write()


ctd = bioparser.data.Data().ctd
name_to_chemical_id = ctd.get_name_to_chemical_id()
chemical_names = set(name_to_chemical_id)
all_names_to_chemical_id = ctd.get_all_names_to_chemical_id()
all_chemical_names = set(all_names_to_chemical_id)

sider = bioparser.data.Data().sider
add_version('sider', sider.directory)
#name_to_sider_drug = sider.get_name_to_drug

sider_mesh_id_tuples = list()
for drug in sider.get_drugs():
    drug_name = drug.name
    chemical_id = None
    if drug_name in chemical_names:
        chemical_id = name_to_chemical_id[drug_name]
        sider_mesh_id_tuples.append((drug, chemical_id))
        continue
    if drug_name in all_chemical_names:
        chemical_id = all_names_to_chemical_id[drug_name]
        sider_mesh_id_tuples.append((drug, chemical_id))
        continue
    all_sider_names = drug.get_all_names(copy_to_lower=True)
    common_names = all_chemical_names & all_sider_names
    if not common_names:
        #print 'No ctd match for SIDER drug:', drug.name
        continue
    chemical_ids = set(all_names_to_chemical_id[name] for name in common_names)
    if len(chemical_ids) == 1:
        chemical_id = chemical_ids.pop()
        sider_mesh_id_tuples.append((drug, chemical_id))
        continue
    for chemical_id in chemical_ids:
        sider_mesh_id_tuples.append((drug, chemical_id))
    #print drug.name, 'mapped to multiple CTD chemicals:', chemical_ids

sider_drug_to_mesh_id = dict(sider_mesh_id_tuples)
tb_ctd_side_effect = Table('ctd_side_effect',
    ['mesh_id', 'umls_id', 'label', 'freq_category_id', 'freq_percentage'])

sider_rows = sider.get_side_effect_rows()
for row in sider_rows:
    mesh_id = sider_drug_to_mesh_id.get(row['drug'])
    if not mesh_id:
        continue
    row['mesh_id'] = mesh_id
    row['umls_id'] = row['umls_concept']
    tb_ctd_side_effect.append(row)

tb_ctd_side_effect_freq_category = factor_table(tb_ctd_side_effect, 'freq_category')
tb_ctd_side_effect_freq_category.write()


# Write holdout tables
tb_resource_version.write()

tb_doid_omim_map.remove_invalid_references(tb_omim, 'omim_id')
tb_doid_omim_map.write()

tb_doid_medic_map.remove_invalid_references(tb_medic, 'medic_id')
tb_doid_medic_map.write()

tb_doid_efo_map.remove_invalid_references(tb_efo, 'efo_id')
tb_doid_efo_map.write()

tb_gene_efo_gwas.remove_invalid_references(tb_efo, 'efo_id')
tb_gene_efo_gwas.write()

tb_drugbank_omim_predict.remove_invalid_references(tb_omim, 'omim_id')
tb_drugbank_omim_predict.write()

tb_ctd_drugbank_map.remove_invalid_references(tb_drugbank, 'drugbank_id')
tb_ctd_drugbank_map.write()

tb_ctd_side_effect.remove_invalid_references(tb_side_effect, 'umls_id') # meddra version should be changed to SIDER version
tb_ctd_side_effect.remove_duplicates(tb_ctd_side_effect.fieldnames)
tb_ctd_side_effect.write()

doc_path = os.path.join(ictnet_dir, 'table_documentation.txt')
doc_file = open(doc_path, 'w')
for table in tables:
    doc_file.write(table.file_name + '\n')
    for fieldname, sql_type in table.schema().items():
        line = '{}: {}\n'.format(fieldname, sql_type)
        doc_file.write(line)
    doc_file.write('\n')
doc_file.close()