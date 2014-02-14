import os
import csv
import gzip


import data
import metathesaurus

from metathesaurus import Concept

class Drug():
        
    def __init__(self, pubchem):
        self.pubchem = pubchem
        self.labels = set()
        self.brand_names = set()
        self.generic_names = set()
        self.stereo_ids = set()
        self.name = None
        #self.side_effect_names = set()
        self.meddra_concepts = set()
        self.meddra_frequent_concepts = set()

    def __hash__(self):
        return hash(self.pubchem)
    
    def __eq__(self, other):
        return self.pubchem == other.pubchem
        
    def __repr__(self):
        return str(self.__dict__)

    def get_all_names(self, copy_to_lower=False):
        """Returns a set with all names for the drug available in SIDER."""
        names = self.brand_names | self.generic_names
        if self.name:
            names.add(self.name)
        if copy_to_lower:
            for name in list(names):
                names.add(name.lower())
        return names

class Label(object):
    
    def __init__(self, label_id):
        self.label_id = label_id
        

class SIDER(object):
    
    def __init__(self, directory=None):
        if not directory:
            directory = data.current_path('sider')
        self.directory = directory
    
    @staticmethod
    def stitch_to_pubchem(flat_id):
        """Convert a stitch flat ID to a pubchem identifier"""
        if not flat_id:
            return None
        return str(abs(int(flat_id)) - 100000000)
    
    def read_file(self, name, fieldnames):
        path = os.path.join(self.directory, '{}.tsv.gz'.format(name))
        read_file = gzip.open(path)
        reader = csv.DictReader(read_file, fieldnames=fieldnames, delimiter='\t')
        for row in reader:
            yield row
        read_file.close()

    def read_label_mapping(self):
        """
        Return a generator of the rows in label_mappings.tsv. Plural fields
        are split. Each row represents a drug label. Most drugs have
        multiple labels. Some labels represent combination drugs that could be
        problematic and should be excluded. The nature of each drug is best
        inferred from the stitch_map_type.
        """
        fieldnames = ['brand_names', 'generic_names', 'stitch_map_type', 
                      'stitch_id_flat', 'stitch_id_stereo', 'url', 'label_id']
        for row in self.read_file('label_mapping', fieldnames):
            for key in ['brand_names', 'generic_names']:
                value = row[key]
                row[key] = value.split(';') if value else list()
            stitch = row['stitch_id_flat']
            row['pubchem'] = SIDER.stitch_to_pubchem(stitch)
            yield row

    def read_adverse_effects_raw(self):
        """Unused"""
        fieldnames = ['label_id', 'concept_id', 'label_side_effect']
        return self.read_file('adverse_effects_raw', fieldnames)

    def read_meddra_adverse_effects(self):
        fieldnames = ['stitch_id_flat', 'stitch_id_stereo', 'label_concept_id',
              'drug_name', 'label_side_effect', 'meddra_level',
              'meddra_concept', 'meddra_name']
        return self.read_file('meddra_adverse_effects', fieldnames)

    def read_meddra_freq_parsed(self):
        fieldnames = ['stitch_id_flat', 'stitch_id_stereo', 'label_id',
                      'umls_concept_id', 'concept_name', 'placebo',
                      'freq', 'freq_lower', 'freq_upper',
                      'meddra_level', 'meddra_concept', 'meddra_name']
        return self.read_file('meddra_freq_parsed', fieldnames)


    def create_drugs(self):
        """
        
        """
        pubchem_to_drug = dict()
        
        for label in self.read_label_mapping():
            pubchem = label['pubchem']
            # excludes combination, unmatched, double matched, and template labels
            if not pubchem:
                continue
            drug = pubchem_to_drug.get(pubchem)
            drug = drug or Drug(pubchem)
            pubchem_to_drug[pubchem] = drug
            drug.labels.add(label['label_id'])
            drug.generic_names |= set(label['generic_names'])
            drug.brand_names |= set(label['brand_names'])
            drug.flat_id = label['stitch_id_flat']
            drug.stereo_ids.add(label['stitch_id_stereo'])
        self.pubchem_to_drug = pubchem_to_drug

        for adverse_effect in self.read_meddra_adverse_effects():
            flat_id = adverse_effect['stitch_id_flat']
            pubchem = SIDER.stitch_to_pubchem(flat_id)
            if not pubchem:
                continue
            drug = pubchem_to_drug[pubchem]
            drug_name = adverse_effect['drug_name']
            if not drug.name:
                drug.name = drug_name
            else:
                assert drug.name == drug_name
            drug.name = drug_name
            if adverse_effect['meddra_level'] != 'PT':
                continue
            drug.meddra_concepts.add(adverse_effect['meddra_concept'])

        for adverse_effect in self.read_meddra_freq_parsed():
            flat_id = adverse_effect['stitch_id_flat']
            pubchem = SIDER.stitch_to_pubchem(flat_id)
            if not pubchem:
                continue
            drug = pubchem_to_drug.get(pubchem)
            if drug is None:
                print 'In Freq Parse: pubchem drug', pubchem, 'not found'
                continue
            if adverse_effect['meddra_level'] != 'PT':
                continue
            if adverse_effect['placebo'] == 'placebo':
                continue
            freq = adverse_effect['freq']
            if freq in ['postmarketing', 'rare', 'infrequent', 'potential']:
                continue
            elif freq == 'frequent':
                drug.meddra_frequent_concepts.add(adverse_effect['meddra_concept'])
            elif '%' in freq:
                percentage = float(freq.replace('%', ''))
                if percentage > 5.0:
                    drug.meddra_frequent_concepts.add(adverse_effect['meddra_concept'])
            else:
                print 'Unparsable frequency', freq




    
    def annotate_meddra_codes(self, umls_version='2011AB', frequency=False):
        drugs = self.pubchem_to_drug.values()
        umls_path = data.version_dir('umls', umls_version)
        meta = metathesaurus.Metathesaurus(umls_path)
        with meta:
            id_to_concept = meta.shelves['concepts']
            for drug in drugs:
                concept_ids = drug.meddra_frequent_concepts if frequency else drug.meddra_concepts
                concepts = set(id_to_concept.get(cui) for cui in concept_ids)
                concepts.discard(None)
                codes = set(concept.source_to_code.get('MDR')
                            for concept in concepts)
                codes.discard(None)
                drug.meddra_codes = codes
        
    def get_drugs(self, exclude_nameless=True):
        drugs = self.pubchem_to_drug.values()
        if exclude_nameless:
            drugs = [drug for drug in drugs if drug.name]
        return drugs
    
    def get_meddra_annotated_drugs(self):
        self.create_drugs()
        self.annotate_meddra_codes()
        return self.get_drugs()

    def get_frequent_meddra_annotated_drugs(self):
        self.create_drugs()
        self.annotate_meddra_codes(frequency=True)
        return self.get_drugs()


    def get_name_to_drug(self):
        name_to_drug = dict()
        for drug in self.get_drugs():
            name_to_drug.update(dict.fromkeys(drug.get_all_names(), drug))
        return name_to_drug()
        
if __name__ == '__main__':
    sider = SIDER()
    sider.create_drugs()
    sider.annotate_meddra_codes(frequency=True)
    for drug in sider.get_drugs():
        if len(drug.meddra_concepts) != len(drug.meddra_codes):
            print drug.name, len(drug.meddra_concepts), len(drug.meddra_codes)
