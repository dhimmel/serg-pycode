import os
import csv
import gzip

import omictools



class CTD(object):
    
    def __init__(self, ctd_dir):
        self.ctd_dir = ctd_dir
        
        self.name_to_mesh = dict() # Populated by self.chemical_vocab_analyzer()
        self.synonym_to_mesh = dict() # Populated by self.chemical_vocab_analyzer()
    
    def chemical_vocab_generator(self):
        """Read CTD_chemicals.tsv.gz"""
        headings = ['ChemicalName', 'ChemicalID', 'CasRN', 'Definition', 
                    'ParentIDs', 'TreeNumbers', 'ParentTreeNumbers', 'Synonyms']
        compound_columns = ['ParentIDs', 'TreeNumbers', 'ParentTreeNumbers', 'Synonyms']
        
        path = os.path.join(self.ctd_dir, 'CTD_chemicals.tsv.gz')
        
        with gzip.GzipFile(path) as f:
            reader = csv.DictReader(omictools.CommentedFile(f), headings, delimiter='\t')
            for row_dict in reader:
                for heading in compound_columns:
                    compound_value = row_dict[heading]
                    row_dict[heading] = compound_value.split('|') if compound_value else list()
                row_dict['ChemicalID'] = row_dict['ChemicalID'][5:] # Remove 'MESH:'
                yield row_dict
    
    def chemical_vocab_analyzer(self):
        chemicals = self.chemical_vocab_generator()
        #keeper_headings = set(['ChemicalName', 'ChemicalID', 'CasRN', 'Definition', 'Synonyms'])
        #chemicals = ({key: value for key, value in chemical.items() if key in keeper_headings} for chemical in chemicals)
        for chemical in chemicals:
            mesh = chemical['ChemicalID']
            name = chemical['ChemicalName']
            name = name.lower()
            self.name_to_mesh[name] = mesh
            for syn in chemical['Synonyms']:
                syn = syn.lower()
                self.synonym_to_mesh[syn] = mesh
    
    def mesh_from_name(self, name):
        if not (self.name_to_mesh and self.synonym_to_mesh):
            self.chemical_vocab_analyzer()
        
        name = name.lower()
        
        if name in self.name_to_mesh:
            return self.name_to_mesh[name]
        
        if name in self.synonym_to_mesh:
            return self.synonym_to_mesh[name]

    
if __name__ == '__main__':
    ctd = CTD('/home/dhimmels/Documents/serg/omicnet/input-datasets/ctd/190912/')
    ctd.chemical_vocab_analyzer()
    
    names = ['rapamycin', 'isoxsuprine', '5-azacytidine', 'fluorescein', 'troglitazone', 'unoprostone', 'depsipeptide', 'pramlintide', 'triazolam', '17-hydroxyprogesterone', 'trovafloxacin', 'pyridostigmine', 'olanzapine', 'piroxicam', 'dobutamine', 'econazole']
    names = []
    for name in names:
        print ctd.mesh_from_name(name)
    
    
    
    
        