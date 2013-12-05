import os
import csv
import gzip

import data


def open_ext(path, *args, **kwargs):
    open_fxn = gzip.open if path.endswith('.gz') else open
    return open_fxn(path, *args, **kwargs)


class CTDReader(csv.DictReader):
    
    def __init__(self, path, plural_fields=list()):
        self.heading_lines = list()
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
        try:
            row = csv.DictReader.next(self)
            for field, value in row.items():
                if field in self.plural_fields:
                    row[field] = value.split('|') if value else list()
                else:
                    if not value:
                        row[field] = None
            return row
        except StopIteration:
            self.close()
            raise StopIteration


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

    def filter_ctd(self, ctd_reader, output_name, keep_row):
        """
        keep_row is a function
        """
        output_path = os.path.join(self.directory, output_name)
        plural_fields = ctd_reader.plural_fields
        write_file = gzip.open(output_path, 'w')
        write_file.write('# Fields:\n# ')
        writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=ctd_reader.fieldnames)
        writer.writeheader()
        for row in ctd_reader:
            if not keep_row(row):
                continue
            row = {k: '|'.join(v) if k in plural_fields else v for k, v in row.items()}
            writer.writerow(row)
        write_file.close()

    def read_file(self, file_name, plural_fields):
        path = os.path.join(self.directory, file_name)
        return CTDReader(path, plural_fields)
        #for row in reader:
        #    yield row
        #reader.close()
    
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

    def read_gene2disease(self, filename='CTD_genes_diseases.tsv.gz'):
        plural_fields = ['DirectEvidence', 'OmimIDs', 'PubMedIDs']
        return self.read_file(filename, plural_fields)

    def read_gene2disease_filtered(self):
        filtered_filename = 'CTD_genes_diseases-filtered.tsv.gz'
        return self.read_gene2disease(filtered_filename)


    def write_filtered_datasets(self):
        self.filter_gene2disease()

    def filter_gene2disease(self):
        filtered_filename = 'CTD_genes_diseases-filtered.tsv.gz'
        ctd_reader = self.read_gene2disease()
        keep_row = lambda row: 'marker/mechanism' in row['DirectEvidence']
        self.filter_ctd(ctd_reader, filtered_filename, keep_row)

if __name__ == '__main__':
    ctd = CTD()
    #ctd.write_filtered_datasets()