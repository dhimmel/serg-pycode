import os
import csv
import gzip
import zipfile

import matplotlib.pyplot as plt

class IMP(object):
    
    def __init__(self, imp_dir):
        self.imp_dir = imp_dir
        
    def read(self, prob_cutoff = 0.0):
        #self.dat_path = os.path.join(self.imp_dir, 'human_global.zip')
        #zip = zipfile.ZipFile(self.dat_path)
        #dat_file = zip.open('global_average_prior.dat')
        #reader = csv.reader(dat_file, delimiter='\t')
        
        self.prob_cutoff = prob_cutoff
        
        self.dat_path = os.path.join(self.imp_dir, 'global_average_prior.dat.gz')
        dat_file = gzip.open(self.dat_path)
        reader = csv.reader(dat_file, delimiter='\t')
        for gene_0, gene_1, prob in reader:
            prob = float(prob)
            if prob < prob_cutoff:
                continue
            yield gene_0, gene_1, prob

if __name__ == '__main__':
    imp = IMP('/home/dhimmels/Documents/serg/omicnet/input-datasets/imp/082812/')
    imp_gen = imp.read(0.8)
    for row in imp_gen:
        print row
    