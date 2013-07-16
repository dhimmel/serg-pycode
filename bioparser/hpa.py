# Human Protein Atlas

import csv
import os

import data

class HPA(object):
    
    def __init__(self, hpa_dir=None):
        if hpa_dir is None:
            hpa_dir = data.current_path('hpa')
        self.hpa_dir = hpa_dir

    def get_row_generator(self, file_name):
        """Returns a generator of dictionaries each representing a row of a
        comma separated, double-quote quoted, header included text file.
        """
        path = os.path.join(self.hpa_dir, file_name)
        csv_file = open(path)
        reader = csv.DictReader(csv_file)
        for row in reader:
            yield row
        csv_file.close()

    def read_normal_tissue(self):
        ensembl_to_tissues = dict()
        rows = self.get_row_generator('normal_tissue.csv')
        valid_levels = {'Strong', 'Moderate'}
        valid_reliabilities = {'Supportive', 'Uncertain'}
        for row in rows:
            #row['Gene']
            #row['Reliability']
            #row['Level']
            #row['Tissue']
            #row['Cell type']
            
            ensembl = row['Gene']
            tissue = row['Tissue']
            # data includes 'soft tissue 1' and 'soft tissue 2'
            if tissue.startswith('soft tissue'):
                tissue = 'soft tissue'
            if row['Level'] not in valid_levels:
                continue
            if row['Reliability'] not in valid_reliabilities:
                continue
            ensembl_to_tissues.setdefault(ensembl, set()).add(tissue)
        return ensembl_to_tissues
        
if __name__ == '__main__':
    hpa = HPA()
    ensembl_to_tissues = hpa.read_normal_tissue()
    tissue_set = set()
    for set_ in ensembl_to_tissues.values():
        tissue_set |= set_
    print tissue_set
    print ensembl_to_tissues['ENSG00000162692'] #CVAM1
