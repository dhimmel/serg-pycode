import os
import csv
import gzip

import data


def open_ext(path, *args, **kwargs):
    open_fxn = gzip.open if path.endswith('.gz') else open
    return open_fxn(path, *args, **kwargs)


class CTDReader(csv.DictReader):
    
    def __init__(self, path, plural_fields=list()):
        self.path = path
        self.plural_fields = set(plural_fields)
        read_file = open_ext(path)
        self.read_file = read_file
        # ignore comments until the line preceeding fieldname declarations
        while read_file.next().rstrip() != '# Fields:':
            pass
        row = read_file.next().rstrip()
        fieldnames = row[2:].split('\t')
        read_file.next()
        csv.DictReader.__init__(self, read_file, delimiter='\t', fieldnames=fieldnames)
    
    def next(self):
        row = csv.DictReader.next(self)
        for field, value in row.items():
            if field in self.plural_fields:
                row[field] = value.split('|') if value else list()
            else:
                if not value:
                    row[field] = None
        return row
        
    def close(self):
        self.read_file.close()




class CTD(object):
    
    def __init__(self, directory=None):
        """
        http://ctdbase.org/downloads/
        """
        if not directory:
            directory = data.current_path('ctd')
        self.directory = directory
    
    
    def read_file(self, file_name, plural_fields):
        path = os.path.join(self.directory, file_name)
        reader = CTDReader(path, plural_fields)
        for row in reader:
            yield row
        reader.close()
    
    def read_chemicals(self):
        plural_fields = ['ParentIDs', 'TreeNumbers', 'ParentTreeNumbers', 
                         'Synonyms', 'DrugBankIDs']
        return self.read_file('CTD_chemicals.tsv.gz', plural_fields)

    def get_all_names_to_chemical_id(self):
        name_to_chemical = dict()
        for chemical in self.read_chemicals():
            names = set(chemical['Synonyms'] + [chemical['ChemicalName']])
            for name in list(names):
                names.add(name.lower())
            for name in names:
                name_to_chemical[name] = chemical['ChemicalID']
        return name_to_chemical

    def get_name_to_chemical_id(self):
        name_to_chemical = dict()
        for chemical in self.read_chemicals():
            name = chemical['ChemicalName']
            for name in [name, name.lower()]:
                name_to_chemical[name] = chemical['ChemicalID']
        return name_to_chemical

    def read_diseases(self):
        plural_fields = ['AltDiseaseIDs', 'ParentIDs', 'TreeNumbers', 
                         'ParentTreeNumbers', 'Synonyms', 'SlimMappings']
        return self.read_file('CTD_diseases.tsv.gz', plural_fields)

    def read_chemical2genes(self):
        plural_fields = ['GeneForms', 'InteractionActions', 'PubMedIDs']
        return self.read_file('CTD_chem_gene_ixns.tsv.gz', plural_fields)

    def read_chemical2diseases(self):
        plural_fields = ['DirectEvidence', 'OmimIDs', 'PubMedIDs']
        return self.read_file('CTD_chemicals_diseases.tsv.gz', plural_fields)

    def read_gene2disease(self):
        plural_fields = ['DirectEvidence', 'OmimIDs', 'PubMedIDs']
        return self.read_file('CTD_genes_diseases.tsv.gz', plural_fields)

    
if __name__ == '__main__':
    ctd = CTD()
    therapy_rows = ctd.read_chemical2diseases()
    disease_ids = set()
    for i, row in enumerate(therapy_rows):
        if 'therapeutic' not in row['DirectEvidence']:
            continue
        disease_ids.add(row['DiseaseID'])
        """
        print '{ChemicalName}\t{DiseaseName}\t{DiseaseID}\t{OmimIDs}'.format(**row)
        #print row['ChemicalName'], row['DiseaseName'], row['OmimIDs'], row['DirectEvidence']
        if i > 1000:
            break
        """
    
    for disease_id in disease_ids:
        if not disease_id.startswith('MESH:'):
            print disease_id
    
    
    
    
        