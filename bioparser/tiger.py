import os
import csv

class TIGER(object):
    
    def __init__(self, tiger_dir=None):
        if tiger_dir is None:
            tiger_dir = data.current_path('tiger')
        self.tiger_dir = tiger_dir

    def get_tissues(self):
        if not hasattr(self, 'tissues'):
            path = os.path.join(self.tiger_dir, 'tissue.txt')
            with open(path) as f:
                reader = csv.DictReader(f, delimiter='\t')
                tissues = {row['Tissue'] for row in reader}
            self.tissues = tissues
        return self.tissues
    
    def get_hs_to_symbol(self):
        if not hasattr(self, 'hs_to_symbol'):
            path = os.path.join(self.tiger_dir, 'symbol2hs-Table.txt')
            unigene_to_symbol = dict()
            with open(path) as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    row['Gene_symbol']
                    row['UniGene(s)']

            self.tissues = tissues
