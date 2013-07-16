import csv
import os


import data

class PaGenBase(object):
    
    def __init__(self, pagen_dir=None):
        if pagen_dir is None:
            pagen_dir = data.source_data_dir('pagenbase')
        self.pagen_dir = pagen_dir
    
    def get_tissue_to_genes(self, spm_min=0.9, include_literature=True):
        """ """
        tissue_to_genes = dict()
        
        path = os.path.join(self.pagen_dir, 'hotisp.txt')
        hotisp_file = open(path)
        # skip the first five lines
        for i in range(5):
            hotisp_file.readline()
        reader = csv.DictReader(hotisp_file, delimiter='\t')
        for row in reader:
            gene = row['Gene Symbol']
            tissue = row['Sample']
            spm = row['SPM']
            if spm == '---':
                if not include_literature:
                    continue
            else:
                spm = float(row['SPM'])
                if spm < spm_min:
                    continue
            tissue_to_genes.setdefault(tissue, set()).add(gene)
        hotisp_file.close()
        return tissue_to_genes
    

if __name__ == '__main__':
    pagen = PaGenBase()
    tissue_to_genes = pagen.get_tissue_to_genes()
    print tissue_to_genes.keys()
    print len(tissue_to_genes)