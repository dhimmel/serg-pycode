import os
import csv

import data

class TIGER(object):
    
    def __init__(self, tiger_dir=None):
        if tiger_dir is None:
            tiger_dir = data.source_data_dir('tiger')
        self.tiger_dir = tiger_dir

    def get_tissues(self):
        if not hasattr(self, 'tissues'):
            path = os.path.join(self.tiger_dir, 'tissue.txt')
            with open(path) as f:
                reader = csv.DictReader(f, delimiter='\t')
                tissues = {row['Tissue'] for row in reader}
            self.tissues = tissues
        return self.tissues
    
    def get_unigene_to_symbol(self):
        if not hasattr(self, 'unigene_to_symbol'):
            path = os.path.join(self.tiger_dir, 'symbol2hs-Table.txt')
            unigene_to_symbol = dict()
            with open(path) as f:
                reader = csv.reader(f, delimiter='\t')
                reader.next() # header
                for row in reader:
                    symbol = row[0]
                    unigenes = row[1:]
                    for unigene in unigenes:
                        unigene_to_symbol[unigene] = symbol
            self.unigene_to_symbol = unigene_to_symbol
        return self.unigene_to_symbol

    def get_gene_to_tissues(self):
        """ """
        if not hasattr(self, 'gene_to_tissues'):
            unigene_to_symbol = self.get_unigene_to_symbol()
            path = os.path.join(self.tiger_dir, 'hs2tissue-Table.txt')
            with open(path) as f:
                reader = csv.reader(f, delimiter='\t')
                reader.next() # header
                gene_to_tissues = dict()
                for row in reader:
                    unigene = row[0]
                    symbol = unigene_to_symbol.get(unigene)
                    if symbol is None:
                        continue
                    tissues = row[1:]
                    for tissue in tissues:
                        gene_to_tissues.setdefault(symbol, set()).add(tissue)
            self.gene_to_tissues = gene_to_tissues
        return self.gene_to_tissues

    

tiger = TIGER()
tiger.get_unigene_to_symbol()
print tiger.get_gene_to_tissues().items()[:40]
