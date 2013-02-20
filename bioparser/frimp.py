import os
import csv
import gzip

import data

class IMP(object):
    
    def __init__(self, imp_dir=None):
        if imp_dir is None:
            imp_dir = data.current_path('imp')

        self.imp_dir = imp_dir
        
        
    def read(self, prob_cutoff = 0.0, symbol=False):
        self.entrez_to_hgnc = data.Data().hgnc.get_entrez_to_gene()  
        self.dat_path = os.path.join(self.imp_dir, 'global_average_prior.dat.gz')
        dat_file = gzip.open(self.dat_path)
        reader = csv.reader(dat_file, delimiter='\t')
        for gene_0, gene_1, prob in reader:
            prob = float(prob)
            if prob < prob_cutoff:
                continue
            if symbol:
                gene_0 = self.entrez_to_hgnc.get(gene_0)
                gene_1 = self.entrez_to_hgnc.get(gene_1)
                if gene_0 is None or gene_1 is None:
                    print 'no symbol for entrez'
                    continue
                gene_0 = gene_0.symbol
                gene_1 = gene_1.symbol
            
            yield gene_0, gene_1, prob
        
   
if __name__ == '__main__':
    
    imp = IMP()
    imp_gen = imp.read(0.8, True)
    for row in imp_gen:
        print row
    