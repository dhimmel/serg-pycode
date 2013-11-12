import os
import csv

import data

class MirBase(object):
    
    def __init__(self, directory=None):
        if not directory:
            directory = data.current_path('mirbase')
        self.directory = directory
    
    
    def alias_generator(self):
        path = os.path.join(self.directory, 'aliases.txt')
        read_file = open(path)
        fieldnames = ['accession', 'identifiers']
        reader = csv.DictReader(read_file, delimiter='\t', fieldnames=fieldnames)
        for row in reader:
            row['identifiers'] = row['identifiers'].split(';')
            yield row
        read_file.close()

    def mirna_generator(self):
        path = os.path.join(self.directory, 'miRNA.tdt')
        read_file = open(path)
        reader = csv.DictReader(read_file, delimiter='\t')
        for row in reader:
            yield row
        read_file.close()

    def load_aliases(self):
        accession_to_identifiers = dict()
        identifier_to_accession = dict()
        for row in self.alias_generator():
            accession = row['accession']
            identifiers = row['identifiers']
            accession_to_identifiers[accession] = identifiers
            for identifier in identifiers:
                identifier_to_accession[identifier] = accession
        self.identifier_to_accession = identifier_to_accession
        self.accession_to_identifiers = accession_to_identifiers
    
    def load_miRNA(self, load_hgnc=True):
        accession_to_mirnas = dict()
        mirnas = list()
        for mirna in self.mirna_generator():
            mirnas.append(mirna)
            accession = mirna['Accession']
            accession_to_mirnas.setdefault(accession, list()).append(mirna)
            
            accession_mature1 = mirna['Mature1_Acc']
            accession_mature2 = mirna['Mature2_Acc']
            for accesion_mature in (accession_mature1, accession_mature2):
                if accesion_mature:
                    accession_to_mirnas.setdefault(accesion_mature, list()).append(mirna)
        self.accession_to_mirnas = accession_to_mirnas
        
        if load_hgnc:
            symbol_to_gene = data.Data().hgnc.get_symbol_to_gene()
            for mirna in mirnas:
                mirna['hgnc'] = symbol_to_gene.get(mirna['ID'])

    def identifier_to_genes(self, identifier):
        if not hasattr(self, 'identifier_to_accession'):
            self.load_aliases()
        if not hasattr(self, 'accession_to_mirnas'):
            self.load_miRNA(load_hgnc=True)
        accession = self.identifier_to_accession.get(identifier)
        if not accession:
            return None
        mirnas = self.accession_to_mirnas.get(accession)
        # this only occurs when accession is in aliases that is not in mirna.tdt
        if not mirnas:
            return None
        genes = {mirna['hgnc'] for mirna in mirnas}
        genes.discard(None)
        return genes
        
        
            

if __name__ == '__main__':
    mirbase = MirBase()
    mirbase.load_miRNA()
    import pprint
    pprint.pprint('')
    


