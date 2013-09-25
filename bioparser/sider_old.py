import os
import csv
import collections

import utilities.omictools
import bioparser.data
import bioparser.shelved

class SIDER(bioparser.shelved.Shelved):
    
    def __init__(self, sider_dir=None):
        """
        indications_raw: The exctacted medical terms from each drug label indications section. CUI included.
        adverse_effects_raw: The exctacted medical terms from each drug label adverse effects section. CUI included.
        label_mapping: Each row represents a drug label.
        meddra_adverse_effects: Maps a drug to meddra adverse effects. Includes separate lines for LLT and PT meddra_types.
        meddra_freq_parsed: Adverse effect frequency information.
        """
        
        if not sider_dir:
            sider_dir = data.current_path('sider')
        self.sider_dir = sider_dir
        
        # Initiate shelves
        shelve_names = ['drugs', 'type_to_concept_ids',
                        'drug_to_indications', 'drug_to_adverse_effects',
                        'indication_to_drugs', 'adverse_effect_to_drugs']
        shelve_dir = os.path.join(self.sider_dir, 'shelves')
        super(SIDER, self).__init__(shelve_dir, shelve_names)
        
        # build if all shelves don't exist
        if not self.shelve_files_exist():
            print 'building sider'
            self.delete_all_shelves()
            self.build()
            print 'sider building complete'
            
    def build(self):
        """
        Parse entiriety of sider
        """
        self.read_labels()
        self.read_raw_file('indications_raw')
        self.read_raw_file('adverse_effects_raw')
        self.remove_empty_labels()
        self.create_drugs()
        self.read_meddra()
        with self:
            all_indications = set()
            all_adverse_effects = set()

            drug_to_indications = self.shelves['drug_to_indications']
            drug_to_adverse_effects = self.shelves['drug_to_adverse_effects']
            
            for drug in self.drugs:
                self.shelves['drugs'][drug.name] = drug
                indications = drug.get_indications()
                adverse_effects = drug.get_adverse_effects()
                all_indications |= set(indications)
                all_adverse_effects |= set(adverse_effects)
                drug_to_indications[drug.name] = indications
                drug_to_adverse_effects[drug.name] = adverse_effects
            
            type_to_concept_ids = self.shelves['type_to_concept_ids']
            type_to_concept_ids['indications'] = all_indications
            type_to_concept_ids['adverse_effects'] = all_adverse_effects
            
            indication_to_drugs = {key: set() for key in all_indications}
            adverse_effect_to_drugs = {key: set() for key in all_adverse_effects}
            for drug in self.drugs:
                drug_name = drug.name
                for concept_id in drug_to_indications[drug_name]:
                    indication_to_drugs[concept_id].add(drug_name)
                for concept_id in drug_to_adverse_effects[drug_name]:
                    adverse_effect_to_drugs[concept_id].add(drug_name)
            
            self.shelves['indication_to_drugs'].update(indication_to_drugs)
            self.shelves['adverse_effect_to_drugs'].update(adverse_effect_to_drugs)           

    def read_file(self, name, plural_headers = ['generic_drug_names', 'brand_drug_names'], plural_sep = ';'):
        """
        Returns a generator of dictionaries each encoding a single row. Column names
        are received from header_dict.
        
        name - name of the file with .tsv. Same as dict key.
        plural_headers - column names that require splitting. Fields that should be lists
        plural_sep - the str separating values in plural fields.
        """
        
        header_dict = {
            'indications_raw':        ['label', 'cui', 'name'],
            'adverse_effects_raw':    ['label', 'cui', 'name'],
            'label_mapping':          ['brand_drug_names', 'generic_drug_names', 'stitch_mapped_successful', 'stitch_flat', 'stitch_stereo', 'url', 'label'],
            'meddra_adverse_effects': ['stitch_flat', 'stitch_stereo', 'label_cui', 'drug_name', 'side_effect_name', 'meddra_type', 'meddra_cui', 'meddra_name'],
            'meddra_freq_parsed':     ['stitch_flat', 'stitch_stereo', 'label', 'label_cui', 'concept_name', 'placebo', 'freq', 'freq_lb', 'freq_up', 'meddra_type', 'meddra_cui', 'meddra_name']
            }

        path = os.path.join(self.sider_dir, name + '.tsv')
        headers = header_dict[name]        
        plural_headers = set(headers) & set(plural_headers)
        singular_headers = set(headers) - set(plural_headers)
        with open(path) as f:
            reader = csv.DictReader(f, fieldnames = headers, delimiter='\t')
            for row_dict in reader:
                for header in plural_headers:
                    unsplit_field = row_dict[header]
                    split_field = unsplit_field.split(plural_sep) if unsplit_field else list()
                    row_dict[header] = split_field
                for header in singular_headers:
                    if row_dict[header] == '':
                        row_dict[header] = None
                yield row_dict

    def read_raw_file(self, name):
        """
        Read either indications adverse_effects_raw.tsv or indications_raw.tsv.
        Returns a tuple of dictionaries.
        name - either 'indications_raw' or 'adverse_effects_raw'
        """
        assert name in ['indications_raw', 'adverse_effects_raw']
        
        label_attr = 'raw_adverse_cuis' if name == 'adverse_effects_raw' else 'raw_indicated_cuis'
        
        findings_raw = self.read_file(name)
        for row_dict in findings_raw:
            label_id = row_dict['label']
            cui = row_dict['cui']
            if label_id not in self.id_to_label:
                continue
            label = self.id_to_label[label_id]
            cuis = getattr(label, label_attr)
            cuis.add(cui)

    def read_labels(self):
        """
        Read labels to create a dictionary representing each label.
        Creates the following object attributes: 
        self.labels      - a list of labels
        self.id_to_label - a dictionary of label_id to label dictionaries
        """
        
        labels = list()

        label_mappings = self.read_file('label_mapping')
        for row_dict in label_mappings:
            if row_dict['stitch_mapped_successful'] is not None:
                continue
            
            label_dict = {
                'id': row_dict['label'],
                'generic_names': row_dict['generic_drug_names'],
                'brand_names': row_dict['brand_drug_names'],
                'stitch_flat': row_dict['stitch_flat'],
                'pubchem': str(abs(int(row_dict['stitch_flat'])) - 100000000),
                'stitch_stereo': row_dict['stitch_stereo']}
            
            label = Label(**label_dict)
            labels.append(label)
                
        self.labels = labels
        self.id_to_label = {label.id : label for label in labels}
            
        print 'Label reading complete:', len(labels), 'labels'

    def create_drugs(self):
        """
        Create drugs from labels
        """
        
        pubchem_to_labels = dict()
        for label in self.labels:
            pubchem_to_labels.setdefault(label.pubchem, set()).add(label)

        drugs = list()
        for pubchem, labels in pubchem_to_labels.iteritems():
            drug = Drug(pubchem)
            drug.generic_names = set()
            drug.brand_names = set()
            drug.stitch_stereo_ids = set()
            drug.indication_counter = collections.Counter()
            drug.adverse_effect_counter = collections.Counter()
            drug.total_labels = len(labels)
            for label in labels:
                drug.generic_names.update(set(label.generic_names))
                drug.brand_names.update(set(label.brand_names))
                drug.stitch_stereo_ids.update(set([label.stitch_stereo]))
                drug.indication_counter.update(label.raw_indicated_cuis)
                drug.adverse_effect_counter.update(label.raw_adverse_cuis)
            
            # Delete adverse effects that are indications
            for finding in drug.indication_counter.keys():
                del drug.adverse_effect_counter[finding]
            
            drugs.append(drug)
        
        self.drugs = drugs
        self.pubchem_to_drug = {drug.pubchem: drug for drug in drugs}
        print len(drugs), 'drugs have been created from labels.'
    
    
    def remove_empty_labels(self):
        """Remove labels that have neither an adverse_effect or indication."""
        labels = self.labels
        nonempty_labels = filter(lambda l: l.raw_adverse_cuis or l.raw_indicated_cuis, labels)
        print len(nonempty_labels), 'out of', len(labels), 'labels have findings'
        self.labels = nonempty_labels

    def remove_empty_drugs(self):
        """Remove labels that have neither an adverse_effect or indication."""
        drugs = self.drugs
        nonempty_drugs = filter(lambda d: d.adverse_effect_counter or d.indication_counter, drugs)
        print len(nonempty_drugs), 'out of', len(drugs), 'drugs have findings'
        self.drugs = nonempty_drugs
    
    def read_meddra(self):
        """ """
        cui_to_name = dict()
        meddra_gen = self.read_file('meddra_adverse_effects')
        for row_dict in meddra_gen:
            if row_dict['meddra_type'] != 'PT':
                continue
            
            #cui_to_name[row_dict['label_cui']] = row_dict['side_effect_name']
            cui_to_name[row_dict['meddra_cui']] = row_dict['meddra_name']
            
            pubchem = str(abs(int(row_dict['stitch_flat'])) - 100000000)
            drug = self.pubchem_to_drug[pubchem]
            drug.name = row_dict['drug_name']
    
        self.meddra_cui_to_name = cui_to_name        
    
class Label(object):
    
    def __init__(self, id, **kwargs):
        self.id = id
        for key, val in kwargs.items():
            setattr(self, key, val)
        self.raw_indicated_cuis = set()
        self.raw_adverse_cuis = set()

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return self.id == other.id

    def __repr__(self):
        return self.id

    def __str__(self):
        return 'Label: id = %s, PubChem = %s, Generic Names = %s' % (self.id, self.pubchem, self.generic_names)


class Drug(object):
    
    def __init__(self, pubchem, **kwargs):
        self.pubchem = pubchem
        self.name = None
        self.mesh = None
        for key, val in kwargs.items():
            setattr(self, key, val)
        
    def filter(self):
        """Delete entries in the counters which do not pass frequency filter."""
        if self.total_labels > 1:
            for counter in self.indication_counter, self.adverse_effect_counter:
                for key in counter.keys():
                    if counter[key] < 2:
                        del counter[key]
    
    def filtered_findings(self, counter):
        """Returns a list of filtered concept_ids from findings with more
        than one annotated labels if the drug has multiple labels.
        """
        if self.total_labels < 2:
            return counter.keys()
        return [concept for concept, count in counter.items() if count > 1]

    def get_indications(self):
        return self.filtered_findings(self.indication_counter)

    def get_adverse_effects(self):
        return self.filtered_findings(self.adverse_effect_counter)
    
    def get_mesh(self):
        """Get mesh id using ctd chemical vocabulary."""
        if self.mesh:
            return self.mesh
        name_list = list()
        if self.name:
            name_list.append(self.name)
        name_list.extend(self.generic_names)
        name_list.extend(self.brand_names)
        for name in name_list:
            mesh = data.Data().ctd.mesh_from_name(name)
            if mesh:
                self.mesh = mesh
                return mesh
    
    def __hash__(self):
        return hash(self.pubchem)
    
    def __eq__(self, other):
        return self.pubchem == other.pubchem
    
    def __str__(self):
        return 'Drug: name = %s, pubchem = %s, total_labels = %s' % (self.name, self.pubchem, self.total_labels)
    
    def __repr__(self):
        return str(self.__dict__)
    

if __name__ == '__main__':

    sider = SIDER()
    with sider:
        pass
        #print sider.shelves['all_terms']['indications']
        #print 'tipranavir'
        #print sider.shelves['drug_to_indications']['tipranavir']
        #print sider.shelves['drug_to_adverse_effects']['tipranavir']
        #print 'C0019158'
        #print sider.shelves['adverse_effect_to_drugs']['C0019158']
        

