import os
import csv

import data

class Etiome(object):
    
    def __init__(self, etiome_dir=None):
        if etiome_dir is None:
            etiome_dir = data.source_data_dir('etiome')

        self.etiome_dir = etiome_dir
    
    def get_disease_to_factors(self, threshold=1):
        """ """
        path = os.path.join(self.etiome_dir, 'flattened_disease_environment.txt')
        disease_to_factors = dict()
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row['value'] < threshold:
                    continue
                disease = row['disease']
                factor = row['environmental_factor']
                disease_to_factors.setdefault(disease, set()).add(factor)
        self.disease_to_factors = disease_to_factors
        return self.disease_to_factors

    def get_factors(self):
        """ """
        path = os.path.join(self.etiome_dir, 'environmental_factors.txt')
        with open(path) as f:
            factors = [line.strip() for line in f]
        return factors
        
        
if __name__ == '__main__':
    etiome = Etiome()
    etiome.get_disease_to_factors()