"""
Daniel Himmelstein
daniel.himmelstein@gmail.com

This script creates tab delimited text files representing ictnet tables. Three
ictnet tables (tb_gene, tb_gene_alias, tb_tissue) are required as inputs. These
input tables are used to map terms to the ictnet id's for genes and tissues.
Input tables should be comma separated text files with " as the quoting
character. This matches myPhpAdmin's default csv export format. 

Set variable ictnet_dir to the desired directory for input and output.
ictnet_dir should contain three directories:

input-tables - contains tb_gene, tb_gene_alias, and tb_tissue
tissue-mappings - contains tissue to nci_name pairs. (tab delimited)
output-tables - directory to write output tables into

Tissue mappings were done manually. A tissue is mapped to a nci thesaurus node
if that node and all of its descendents are active in the tissue. Tissues 
should be mapped to the highest level nci node possible to make mapping easier.
"""

import os
import csv

import omictools
import data
import metamap

from metathesaurus import Concept
from sider import Drug

ictnet_dir = '/home/dhimmels/Documents/serg/ictnet/ictnet-creation/'

class TableWriter(object):
    
    def __init__(self, table_name):
        """Class to write tab delimited text files encoding ictnet tables."""
        self.table_name = table_name
        self.path = os.path.join(ictnet_dir, 'output-tables', table_name + '.txt')
        self.file = open(self.path, 'wb')
        self.writer = csv.writer(self.file, delimiter='\t')
    
    def writerow(self, row):
        self.writer.writerow(row)
    
    def writerows(self, rows, sort=True):
        if sort:
            rows = list(rows)
            rows.sort()
        for row in rows:
            self.writerow(row)
        
    def close(self):
        print 'Writing', self.table_name, 'is complete.'
        self.file.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, *args, **kwargs):
        self.close()

class TableReader(object):
    
    def __init__(self, table_name):
        self.table_name = table_name
        self.path = os.path.join(ictnet_dir, 'input-tables', table_name + '.csv')
    
    def __enter__(self):
        self.f = open(self.path)
        reader = csv.reader(self.f, quotechar='"')
        return reader
    
    def __exit__(self, *args, **kwargs):
        self.f.close()

# Store neccessary data sources
nci = data.Data().nci
nci.read_nci_thesaurus()
sider = data.Data().sider
metathesaurus = data.Data().metathesaurus
omim = data.Data().omim
omim.read()

With = omictools.With

################################################################################
################################ Side Effects ##################################

# Get sets of all SIDER side effects and indications. Encoded by concept_ids
with With(sider, metathesaurus):
    all_indications = set(sider.shelves['type_to_concept_ids']['indications'])
    all_side_effects = set(sider.shelves['type_to_concept_ids']['adverse_effects'])
    # Remove indications or side_effects that are not in the metathesaurus
    concept_id_to_concept = metathesaurus.shelves['concepts']
    for set_ in all_indications, all_side_effects:
        for concept_id in list(set_):
            if concept_id not in concept_id_to_concept:
                #print concept_id, 'not found in thesaurus'
                set_.remove(concept_id)

# Create and write tb_side_effect - side effect concept_ids and their names
with With(metathesaurus, TableWriter('tb_side_effect')) as (metathesaurus, table):
    concept_id_to_concept = metathesaurus.shelves['concepts']
    for side_effect in list(all_side_effects):
        concept = concept_id_to_concept[side_effect]
        name = concept.name
        table.writerow([side_effect, name])

# Create and write tb_drug_side_effect_net - drug mesh id and side effect pairs
with With(sider, TableWriter('tb_drug_side_effect_net')) as (sider, table):
    for drug in sider.shelves['drugs'].itervalues():
        mesh_id = drug.get_mesh()
        if not mesh_id:
            continue
        
        # Extract Adverse Effects
        concept_ids = drug.get_adverse_effects()
        for concept_id in concept_ids:
            if concept_id in all_side_effects:
                table.writerow([concept_id, mesh_id])

################################################################################
###### Map side effects and indications to NCI thesaurus terms #################
sider_map_dir = os.path.join(sider.sider_dir, 'nci-mappings')
sider_map_dir_known = os.path.join(sider.sider_dir, 'nci-mappings-known')

indication_mapper = metamap.MetaMap(sider_map_dir, sider_map_dir_known,
    fname_suffix='indications', restrict_to_sources=['NCI'],
    restrict_to_types=['dsyn','neop','anab'])
side_effect_mapper = metamap.MetaMap(sider_map_dir, sider_map_dir_known,
    fname_suffix='side_effects', restrict_to_sources=['NCI'],
    restrict_to_types=['fndg', 'sosy'])

# Find mappings from UMLS to treat as previousely manually added mappings
for kind in ['indications', 'side_effects']:
    if kind == 'indications':
        mapper = indication_mapper
        all_concept_ids = all_indications
        valid_types = {'Disease or Syndrome',  'Anatomical Abnormality', 'Neoplastic Process'}

    if kind == 'side_effects':
        mapper = side_effect_mapper
        all_concept_ids = all_side_effects
        valid_types = {'Finding',  'Sign or Symptom'}
    
    with metathesaurus:
        nci_concept_ids = set(metathesaurus.shelves['sources']['NCI'])
        concept_id_to_concept = metathesaurus.shelves['concepts']
        mappings = set()
        for concept_id in all_concept_ids & nci_concept_ids:
            concept = concept_id_to_concept[concept_id]
            if not set(concept.symantic_types) & valid_types:
                continue
                print concept.name
            mapping = metamap.Mapping(input=concept.name, 
                                      concept_name=concept.source_to_name['NCI'],
                                      concept_id=concept_id, score=1001)
            mappings.add(mapping)
        mapper.write_mapping_file(mappings, mapper.paths['previous']['added'])
        concept_ids_to_metamap = all_concept_ids - nci_concept_ids
        terms = set(concept_id_to_concept[concept_id].name
                    for concept_id in concept_ids_to_metamap)
    
    if not mapper.maps['current']['redacted']:
        mapper.map_terms(terms)
        mapper.pause_before_reading_mappings()

# Write tb_side_effect_finding_map
with With(sider, metathesaurus, TableWriter('tb_side_effect_finding_map')) as (sider, metathesaurus, table):
    side_effect_name_to_concept_id = side_effect_mapper.get_mapping_dict()
    concept_id_to_concept = metathesaurus.shelves['concepts']
    for side_effect in all_side_effects:
        name = concept_id_to_concept[side_effect].name
        nci_concept_id = side_effect_name_to_concept_id.get(name)
        if not nci_concept_id:
            continue
        nci_code = concept_id_to_concept[nci_concept_id].source_to_code['NCI']
        table.writerow([side_effect, nci_code])

with With(sider, metathesaurus, TableWriter('tb_disease_drug_net')) as (sider, metathesaurus, table):    
    # Create dicitonary of indication cuis (from sider) to nci_code
    indication_to_nci_code = dict()
    concept_id_to_concept = metathesaurus.shelves['concepts']    
    indication_to_concept_id = indication_mapper.get_mapping_dict()
    for indication in all_indications:
        term = concept_id_to_concept[indication].name
        concept_id = indication_to_concept_id.get(term)
        if not concept_id:
            continue
        nci_code = concept_id_to_concept[concept_id].source_to_code['NCI']
        indication_to_nci_code[indication] = nci_code

    # Write tb_disease_drug_net
    indication_to_concept_id = indication_mapper.get_mapping_dict()
    rows = set()
    for drug in sider.shelves['drugs'].itervalues():
        mesh_id = drug.get_mesh()
        if not mesh_id:
            continue
        # Extract Adverse Effects
        concept_ids = drug.get_indications()
        for concept_id in concept_ids:
            nci_code = indication_to_nci_code.get(concept_id)
            if nci_code:
                rows.add((mesh_id, nci_code))
    table.writerows(rows)
            

################################################################################
############################# Tissue Mappings ##################################

def tissue_to_nci_thesaurus_mapping(tissue_to_nci_fname, root_node_name,
                                    naming, tb_net_name):
    #Function to create tissue mappings and NCI Thesaurus tables.
    
    root_node = nci.concept_name_dict[root_node_name]
    nodes = nci.get_all_descendants(root_node, include_self=True)
    node_names = set(node.concept_name for node in nodes)

    tb_term_ontology_naming = TableWriter('tb_term_ontology_' + naming)
    tb_ontology_naming = TableWriter('tb_ontology_' + naming)
    tb_naming_alias = TableWriter('tb_' + naming + '_alias')

    for node in nodes:
        tb_term_ontology_naming.writerow([node.code, node.concept_name,
                                          node.definition])
        for child in node.children:
            tb_ontology_naming.writerow([child.code, node.code])
        for alias in node.synonyms:
            tb_naming_alias.writerow([node.code, alias])

    tb_term_ontology_naming.close()
    tb_ontology_naming.close()
    tb_naming_alias.close()

    # Read file with mappings of nci_symptom super nodes to associated tissues
    nci_name_to_tissues = dict()
    path = os.path.join(ictnet_dir, 'tissue-mappings', tissue_to_nci_fname)
    with open(path) as f:
        reader = csv.reader(f, delimiter='\t')
        for tissue, nci_name in reader:
            nci_name_to_tissues.setdefault(nci_name, set()).add(tissue)

    nci.add_annotation_type('tissues')
    for nci_name, tissues in nci_name_to_tissues.items():
        node = nci.concept_name_dict[nci_name]
        node.direct_annots['tissues'] |= set(tissues)
        descendants = nci.get_all_descendants(node, include_self=True)
        for descendant in descendants:
            descendant.prop_annots['tissues'] |= set(tissues)

    with TableWriter(tb_net_name) as tb_net:
        for node in nodes:
            for tissue in node.prop_annots['tissues']:
                tissue_id = tissue_name_to_id[tissue]
                tb_net.writerow([tissue_id, node.code])


# Read tb_tissue: table with tissue_id and tissue_name
with TableReader('tb_tissue') as reader:
    tissue_name_to_id = {row[1]: row[0] for row in reader}

tissue_to_nci_thesaurus_mapping(tissue_to_nci_fname = 'tissue_to_nci_finding.txt',
                                root_node_name = 'Sign_or_Symptom',
                                naming = 'finding',
                                tb_net_name = 'tb_tissue_finding_net')

tissue_to_nci_thesaurus_mapping(tissue_to_nci_fname = 'tissue_to_nci_disease.txt',
                                root_node_name = 'Diseases_and_Disorders',
                                naming = 'disease',
                                tb_net_name = 'tb_disease_tissue_net')

################################################################################
######################################## OMIM ##################################

# Table of omim disorder number's and their corresponding names
with TableWriter('tb_omim_disease') as table:
    rows = omim.disorder_number_to_name.iteritems()
    table.writerows(rows)

# Table of omim disorder number's to implicated genes
gene_symbol_to_ictnet_id = dict()
for table_name in ['tb_gene', 'tb_gene_alias']:
    with TableReader(table_name) as reader:
        for row in reader:
            gene_symbol_to_ictnet_id[row[1]] = row[0]

with TableWriter('tb_omim_report') as table:
    for disorder_num, gene_symbols in omim.disorder_num_to_genes.iteritems():
        for gene_symbol in gene_symbols:
            ictnet_gene_id = gene_symbol_to_ictnet_id.get(gene_symbol)
            if ictnet_gene_id:
                table.writerow([disorder_num, ictnet_gene_id])

# Mapping between OMIM terms and NCI Thesaurus
omim_map_dir = os.path.join(omim.omim_dir, 'nci_mappings')
omim_mapper = metamap.MetaMap(omim_map_dir,
    restrict_to_sources=['NCI'],
    restrict_to_types=['dsyn','neop','anab'])

if not omim_mapper.maps['current']['redacted']:
    print 'MetaMapping OMIM to NCI.'
    concepts_ids_to_remove = {'C0162429', 'C1458156', 'C0027627', 'C0014236',
                              'C0027651', 'C0012634', 'C0039082', 'C1333600',
                              'C0334044', 'C0151514'}
    omim_mapper.map_terms(omim.disorder_name_to_number.keys(), concepts_ids_to_remove=concepts_ids_to_remove)
    omim_mapper.pause_before_reading_mappings()

omim_name_to_concept_id = omim_mapper.get_mapping_dict()
with With(metathesaurus, TableWriter('tb_omim_disease_ontology_map')) as (metathesaurus, table):
    concept_id_to_concept = metathesaurus.shelves['concepts']
    for omim_name, concept_id in omim_name_to_concept_id.iteritems():
        omim_number = omim.disorder_name_to_number[omim_name]
        concept = concept_id_to_concept[concept_id]
        nci_code = concept.source_to_code['NCI']
        table.writerow([nci_code, omim_number])
    
