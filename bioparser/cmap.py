import csv
import os

class CMAP(object):
    
    def __init__(self, current_dir, previous_dir=None):
        self.current_dir = current_dir
        self.previous_dir = previous_dir
    
    @staticmethod
    def read_gmt(path):
        """Read a Gene Matrix Transposed file format."""
        set_name_to_genes = dict()
        with open(path) as gmt_file:
            reader = csv.reader(gmt_file, delimiter = '\t')
            for row in reader:
                name, desc = row[:2]
                genes = set(row[2:])
                set_name_to_genes[name] = genes
        return set_name_to_genes

    def read_cmap(self):
        expr = dict()
        for direction in 'up', 'dn':
            path = os.path.join(self.current_dir, 'msigdb_gene_sets', 'msigdb_%s_mapped_to_HG_U133A.gmt' % direction)
            set_name_to_genes = self.read_gmt(path)
            set_name_to_genes = {key[:-3]: value for key, value in set_name_to_genes.iteritems()}
            for set_name, genes in set_name_to_genes.iteritems():
                expr.setdefault(set_name, dict())[direction] = genes
        self.expr = expr
        return expr
        
if __name__ == '__main__':
    cmapper = CMAP('/home/dhimmels/Documents/serg/omicnet/input-datasets/cmap/120824')
    cmapper.read_cmap()
    print len(cmapper.expr['ALZHEIMERS_DISEASE']['dn'])