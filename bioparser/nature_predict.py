import xlrd
import os

import data

class NaturePredict(object):
    """http://www.nature.com/msb/journal/v7/n1/full/msb201126.html"""
    
    def __init__(self, directory=None):
        if not directory:
            directory = data.source_data_dir('nature-predict')
        self.directory = directory
        self.indications_path = os.path.join(directory, 'msb201126-s1.xls')
        self.disease_mappings_path = os.path.join(directory, 'msb201126-s4.xls')
    
    def read_indications(self):
        """
        Read indications excel file and return an omim disease name to 
        drugbank name dictionary.
        """
        
        wb = xlrd.open_workbook(self.indications_path)
        sheet = wb.sheet_by_name(u'Drug indications')
        
        drug_to_diseases, disease_to_drugs = dict(), dict()
        
        # Seven drug names had the first letter of their second word incorrectly capitalized.
        # Two drug names appear to have omitted their second word.
        # One drug contained the brand instead of generic name.
        drug_conversion_dict = {'Arsenic Trioxide': 'Arsenic trioxide',
                                'Meclofenamic Acid': 'Meclofenamic acid',
                                'Ipratropium': 'Ipratropium bromide',
                                'Salicyclic Acid': 'Salicyclic acid',
                                'Adenosine Monophosphate': 'Adenosine monophosphate',
                                'Ethacrynic Acid': 'Ethacrynic acid',
                                'Divalproex Sodium': 'Valproic Acid', # Valporic acid is the generic
                                'Methyl Aminolevulinate': 'Methyl aminolevulinate',
                                'Fondaparinux Sodium': 'Fondaparinux sodium',
                                'Beclomethasone': 'Beclometasone dipropionate'}
        
        for row_num in range(1, sheet.nrows):
            
            drug, disease = sheet.row_values(row_num)
            
            if drug in drug_conversion_dict:
                drug = drug_conversion_dict[drug]
            
            drug_to_diseases.setdefault(drug, set()).add(disease)
            disease_to_drugs.setdefault(disease, set()).add(drug)
        
        self.drug_to_diseases = drug_to_diseases
        self.disease_to_drugs = disease_to_drugs

    def read_disease_mappings(self):
        wb = xlrd.open_workbook(self.disease_mappings_path)
        sheet = wb.sheet_by_name(u'OMIM to UMLS mapping')
        
        column_names = ['omim_id', 'omim_name', 'concept_id', 'concept_name']
        rows = list()
        for row_num in range(1, sheet.nrows):
            row = sheet.row_values(row_num)
            row_dict = dict(zip(column_names, row))
            if row_dict['omim_name'] == "Neuropathy, Hereditary Sensory And Autonomic, Type I, With Cough And":
                row_dict['omim_name'] = "Neuropathy, Hereditary Sensory And Autonomic, Type I, With Cough And Gastroesophageal Reflux"
            row_dict['omim_id'] = str(int(row_dict['omim_id']))
            rows.append(row_dict)
        self.disease_mappings = rows

    def omim_drugbank_mapper(self):
        self.read_disease_mappings()
        self.read_indications()
        data.Data().drugbank.read(True, False, False)
        name_to_drugbank = data.Data().drugbank.get_name_to_drug()
        rows = list()
        omim_id_to_name = {disease_mapping['omim_id']: disease_mapping['omim_name']
            for disease_mapping in self.disease_mappings}
        for omim_id, omim_name in omim_id_to_name.items():
            drugs = self.disease_to_drugs[omim_name]
            drugbanks = filter(lambda x: x, (name_to_drugbank.get(drug) for drug in drugs))
            #print omim_name, len(drugs), len(drugbanks)
            for drugbank in drugbanks:
                row = {'omim_id': omim_id, 'omim_name': omim_name,
                 'drugbank_id': drugbank['drugbank_id'], 'drugbank_name': drugbank['name']}
                rows.append(row)
        return rows

    def read(self):
        print 'Reading Nature PREDICT indications:'
        self.read_indications()
        self.read_disease_mappings()

        
        concept_id_to_indicated_drugs = dict()
        for disease_mapping in self.disease_mappings:
            concept_id = disease_mapping['concept_id']
            omim_name = disease_mapping['omim_name']
            drugs = self.disease_to_drugs[omim_name]
            concept_id_to_indicated_drugs.setdefault(concept_id, set()).update(drugs)
        
        self.concept_id_to_indicated_drugs = concept_id_to_indicated_drugs
        print len(concept_id_to_indicated_drugs), 'diseases with indicated drugs.'
        print sum(map(len, concept_id_to_indicated_drugs.values())), 'total indications.'
        return concept_id_to_indicated_drugs

if __name__ == '__main__':
    np = NaturePredict()
    
    #np.read_disease_mappings()
    rows = np.omim_drugbank_mapper()
    import pprint
    #pprint.pprint(np.disease_mappings)

    pprint.pprint(rows)

    
    #np.read()
    
    
    
    