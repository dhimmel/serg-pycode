import os
import csv

import omictools
import data

class MedDRA(object):
    
    def __init__(self, release_name):
        """Read MedDRA"""
        self.release_path = os.path.join(data.data_dir, 'meddra', release_name)
        self.name_to_code = None
        self.code_to_name = None
        self.codes = None
        
    def read_pt(self):
        """Read pt.asc in MedAscii directory."""
        pt_headers = ['hlt_code', 'hlt_name', 'hlt_whoart_code', 'hlt_harts_code', 
         'hlt_costart_sym', 'hlt_icd9_code', 'hlt_icd9cm_code', 
         'hlt_icd10_code', 'hlt_jart_code']
        path = os.path.join(self.release_path, 'MedAscii', 'pt.asc')
        with open(path) as f:
            reader = csv.DictReader(f, fieldnames=pt_headers, delimiter='$')
            code_to_name = {row['hlt_code']: row['hlt_name'] for row in reader}

        self.codes = set(code_to_name)        
        self.code_to_name = code_to_name
        self.name_to_code = omictools.inverse_dict(code_to_name)
    
    def code_from_name(self, name):
        if self.name_to_code is None:
            self.read_pt()
        try:
            return self.name_to_code[name]
        except KeyError:
            return None

    def name_from_code(self, code):
        if self.code_to_name is None:
            self.read_pt()
        try:
            return self.code_to_name[code]
        except KeyError:
            return None
    
    def get_codes(self):
        if self.codes is None:
            self.read_pt()
        return self.codes
    
    def write_code_to_name_table(self, path, header=False):
        with open(path, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            if header:
                writer.writerow(['code', 'name'])
            if not self.code_to_name:
                self.read_pt()
            code_name_tuples = self.code_to_name.items()
            code_name_tuples.sort()
            writer.writerows(code_name_tuples)
        

    
if __name__ == '__main__':
    meddra = MedDRA('MedDRA_15_1_English')
    meddra.read_pt()
    print meddra.code_to_name