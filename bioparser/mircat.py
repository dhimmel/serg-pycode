import os
import csv

import data


class MirCat(object):
    
    def __init__(self, directory=None):
        """
        http://www.mirrna.org/
        """
        if not directory:
            directory = data.source_data_dir('mircat')
        self.directory = directory
    
    def read_file(self, name):
        path = os.path.join(self.directory, name)
        with open(path) as read_file:
            reader = csv.DictReader(read_file, delimiter='\t')
            for row in reader:
                yield row

    def read_mirna(self):
        return self.read_file('mirna.txt')

    def read_targets(self):
        hgnc = data.Data().hgnc
        symbol_to_gene = hgnc.get_symbol_to_gene()
        for row in self.read_file('target.txt'):
            gene = symbol_to_gene.get(row['gene_symbol'])
            if not gene:
                continue
            row['gene'] = gene
            yield row

        
if __name__ == '__main__':
    mircat = MirCat()
