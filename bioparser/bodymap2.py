import os
import csv
import collections

import data

class BodyMap2(object):
    
    def __init__(self, directory=None):
        if not directory:
            directory = data.source_data_dir('bodymap2')
        self.directory = directory
        self.path = os.path.join(directory, 'E-MTAB-513-query-results.tsv')
    
    def read(self):
        read_file = open(self.path)
        line = '#'
        while line.startswith('#'):
            line = read_file.next().rstrip()
        fieldnames = line.split('\t')
        reader = csv.DictReader(read_file, delimiter='\t', fieldnames=fieldnames)
        tissues = list(fieldnames)
        for non_tissue_key in ['Gene name', 'Gene Id']:
            tissues.remove(non_tissue_key)
        rows = list()
        for row in reader:
            for tissue in tissues:
                fpkm = row[tissue]
                row[tissue] = float(fpkm) if fpkm != '' else None
            rows.append(row)
        read_file.close()
        return tissues, rows
    
    def process(self):
        tissues, rows = self.read()
        #symbol_to_gene = data.Data().hgnc.get_symbol_to_gene()
        ensembl_to_gene = data.Data().hgnc.get_ensembl_to_gene()
        gene_to_rows = dict()
        for row in rows:
            #symbol = row['Gene name']
            ensembl = row['Gene Id']
            #gene = symbol_to_gene.get(symbol)
            gene = ensembl_to_gene.get(ensembl)
            if not gene:
                continue
            gene_to_rows.setdefault(gene, list()).append(row)
        processed_rows = list()
        for gene, rows in gene_to_rows.iteritems():
            processed_row = collections.OrderedDict()
            processed_row['symbol'] = gene.symbol
            for tissue in tissues:
                fpkms = [row[tissue] for row in rows if row[tissue] is not None]
                fpkm = sum(fpkms) / len(fpkms) if fpkms else None
                processed_row[tissue] = fpkm
            processed_rows.append(processed_row)
        
        processed_rows.sort(key=lambda row: row['symbol'])
        
        path = os.path.join(self.directory, 'processed.txt')
        with open(path, 'w') as write_file:
            writer = csv.writer(write_file, delimiter='\t')
            writer.writerow(processed_rows[0].keys())
            for row in processed_rows:
                writer.writerow(row.values())
        
        return processed_rows
    
    def read_processed(self):
        path = os.path.join(self.directory, 'processed.txt')
        read_file = open(path)
        reader = csv.DictReader(read_file, delimiter='\t')
        tissues = list(reader.fieldnames)
        tissues.remove('symbol')
        rows = list()
        for row in reader:
            for tissue in tissues:
                fpkm = row[tissue]
                row[tissue] = float(fpkm) if fpkm != '' else None
            rows.append(row)
        read_file.close()
        return tissues, rows

    def get_edges(self, fpkm_cutoff = 0):
        tissues, rows = self.read_processed()
        for row in rows:
            symbol = row['symbol']
            for tissue in tissues:
                fpkm = row[tissue]
                if fpkm is None:
                    continue
                if fpkm < fpkm_cutoff:
                    continue
                edge = symbol, tissue, fpkm
                yield edge
        
if __name__ == '__main__':
    bodymap = BodyMap2()
    #bodymap.process()
    edges = list(bodymap.get_edges(100))
    print edges[:100]