# Human Protein Atlas

import csv
import os

import data

class HPA(object):
    
    def __init__(self, hpa_dir=None):
        if hpa_dir is None:
            hpa_dir = data.current_path('hpa')
        self.hpa_dir = hpa_dir

    def get_row_generator(self, file_name):
        """Returns a generator of dictionaries each representing a row of a
        comma separated, double-quote quoted, header included text file.
        """
        path = os.path.join(self.hpa_dir, file_name)
        csv_file = open(path)
        reader = csv.DictReader(csv_file)
        for row in reader:
            yield row
        csv_file.close()

    def read_normal_tissue(self, gene_symbols=True):
        gene_to_tissues = dict()
        rows = self.get_row_generator('normal_tissue.csv')
        valid_levels = {'Strong', 'Moderate'}
        valid_reliabilities = {'Supportive', 'Uncertain'}
        if gene_symbols:
            ensembl_to_gene = data.Data().hgnc.get_ensembl_to_gene()
            unmapped_ensembl_ids = set()

        for row in rows:
            #cell_type = row['Cell type']            
            gene = row['Gene']
            if gene_symbols:
                gene = ensembl_to_gene.get(gene)
                if not gene:
                    unmapped_ensembl_ids.add(row['Gene'])
                    continue
                gene = gene.symbol
            tissue = row['Tissue']
            
            # data includes 'soft tissue 1' and 'soft tissue 2'
            if tissue.startswith('soft tissue'):
                tissue = 'soft tissue'
            if row['Level'] not in valid_levels:
                continue
            if row['Reliability'] not in valid_reliabilities:
                continue
            gene_to_tissues.setdefault(gene, set()).add(tissue)
        if gene_symbols:
            print len(unmapped_ensembl_ids), 'ensembl ids in human protein atlas did not map to HGNC genes.'
        return gene_to_tissues
        
if __name__ == '__main__':
    hpa = HPA()
    gene_to_tissues = hpa.read_normal_tissue()
    tissue_set = set()
    for set_ in gene_to_tissues.values():
        tissue_set |= set_
    print tissue_set
    print gene_to_tissues['HHEX']
