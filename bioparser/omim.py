import csv
import os
import re

import data
import omictools

class OMIM(object):
    
    def __init__(self, omim_dir=None):
        if not omim_dir:
            omim_dir = data.current_path('omim')
        self.omim_dir = omim_dir
        
        self.mim_to_gene_symbol = dict()
        self.disorder_to_genes = dict()
    
    def read(self):
        self.read_mim2gene()
        self.compute_disorder_to_genes()
    
    def read_mim2gene(self):
        """A field-delimited file linking MIM numbers to NCBI Gene IDs and HGNC
        Approved Gene Symbols.
        """
        path = os.path.join(self.omim_dir, 'mim2gene.txt')
        mim_to_gene_symbol = dict()
        with open(path) as f:
            f = omictools.CommentedFile(f)
            fieldnames = ['mim_number', 'type', 'gene_id', 'gene_symbol']
            reader = csv.DictReader(f, delimiter='\t', fieldnames=fieldnames)
            for row in reader:
                if row['gene_symbol'] != '-' and 'move' not in row['type']:
                    mim_to_gene_symbol[row['mim_number']] = row['gene_symbol']
        self.mim_to_gene_symbol = mim_to_gene_symbol
                
    
    def read_morbidmap(self):
        path = os.path.join(self.omim_dir, 'morbidmap')
        with open(path) as f:
            fieldnames = ['disorder', 'gene_symbols', 'gene_mim_number', 'location']
            reader = csv.DictReader(f, delimiter='|', fieldnames=fieldnames)
            mapping_key_pattern = re.compile("\([1-4]\)\s*$")
            mim_number_pattern = re.compile(",*\s*[0-9]{6}$") # assumes mim numbers are 6 digits
            
            remove_phrases = ['susceptibility to', 'susceptitbility to',
                              'resistance to', 'modification of',
                              'susceptibility', 'progression of',
                              'association with', 'protection against',
                              'modifier of severity of', 'susceptiblity',
                              'modifier of', 'poor response to', 'hereditary',
                              'reduced risk of', 'decreased', 'infection by',
                              'disease progression']
            
            for row in reader:
                row['gene_symbols'] = row['gene_symbols'].split(', ')
                disorder = row['disorder']
                if 'qtl' in disorder.lower():
                    continue

                occurances = mapping_key_pattern.findall(disorder)
                phene_mapping_key = occurances[0][1] if occurances else None
                disorder = mapping_key_pattern.sub('', disorder).rstrip()
                occurances = mim_number_pattern.findall(disorder)
                mim_disorder_number = occurances[0].strip().strip(', ').strip() if occurances else None
                disorder = mim_number_pattern.sub('', disorder).rstrip()
                disorder = disorder.strip('[{}]?')
                disorder = disorder.strip('0123456789')
                for phrase in remove_phrases:
                    disorder = disorder.replace(phrase, '')
                disorder = disorder.rstrip(' ,-/')
                yield disorder, mim_disorder_number, row['gene_mim_number']
    
    def compute_disorder_to_genes(self):
        disorder_name_to_mim_genes = dict()
        disorder_name_to_number = dict()
        for disorder_name, disorder_number, mim_gene in self.read_morbidmap():
            disorder_name_to_mim_genes.setdefault(disorder_name, set()).add(mim_gene)
            if disorder_number:
                disorder_name_to_number[disorder_name] = disorder_number
        
        self.disorder_name_to_number = disorder_name_to_number
        self.disorder_number_to_name = omictools.inverse_dict(disorder_name_to_number)


        disorder_num_to_genes = dict()
        for disorder_name, mim_genes in disorder_name_to_mim_genes.iteritems():
            disorder_num = disorder_name_to_number.get(disorder_name)
            if not disorder_num:
                continue
            genes = set()
            for mim_gene in mim_genes:
                gene = self.mim_to_gene_symbol.get(mim_gene)
                if gene:
                    genes.add(gene)
            if genes:
                disorder_num_to_genes.setdefault(disorder_num, set()).update(genes)
        self.disorder_num_to_genes = disorder_num_to_genes


if __name__ == '__main__':
    omim = OMIM()
    omim.read()
    #print omim.disorder_name_to_number
    #print len(omim.disorder_name_to_number)
    #print sum(bool(x) for x in omim.disorder_name_to_number.values())
    
    
    