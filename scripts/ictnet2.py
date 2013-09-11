import csv
import os

import bioparser.data

ictnet_dir = '/home/dhimmels/Documents/serg/ictnet/ictnet-creation2/'


"""
Do we want approved symbol included in tb_gene_alias?
Are we excluding non-protein coding genes?
"""
class TableWriter(object):
    
    def __init__(self, table_name):
        """Class to write tab delimited text files encoding ictnet tables."""
        self.table_name = table_name
        self.path = os.path.join(ictnet_dir, 'tables', table_name + '.txt')
        self.file = open(self.path, 'w')
        self.writer = csv.writer(self.file, delimiter='\t')
    
    def writerow(self, row):
        self.writer.writerow(row)
    
    def writerows(self, rows, sort=True):
        if sort:
            rows = list(rows)
            rows.sort()
        for row in rows:
            self.writerow(row)
        
    def close(self):
        print 'Writing', self.table_name, 'is complete.'
        self.file.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, *args, **kwargs):
        self.close()


################################################################################
############# genes
hgnc = bioparser.data.Data().hgnc
symbol_to_gene = hgnc.get_symbol_to_gene()
tb_gene_rows = list() # (gene_id, symbol, name, chromosome)
tb_gene_alias_rows = list() # (gene_id, alias)
genes = hgnc.get_genes()
for gene in genes:
    gene_id = gene.id_
    row = (gene_id, gene.symbol, gene.name, gene.chromosome)
    tb_gene_rows.append(row)
    for alias in gene.synonyms + gene.previous_symbols:
        tb_gene_alias_rows.append((gene_id, alias))

# write tb_gene_rows
with TableWriter('tb_gene') as writer:
    writer.writerow(('gene_id', 'symbol', 'name', 'chromosome'))
    writer.writerows(tb_gene_rows)

# write tb_gene_alias_rows
with TableWriter('tb_gene_alias') as writer:
    writer.writerow(('gene_id', 'alias'))
    writer.writerows(tb_gene_alias_rows)

################################################################################
############# ppi
tb_ppi_rows = list() #(source, taget, pubmed, method, sources)

ppitrim = bioparser.data.Data().ppitrim
complex_to_rows, binary_rows = ppitrim.separate_complexes()
interactions = ppitrim.binary_rows_to_interactions(binary_rows)
interactions = ppitrim.interaction_hgnc_convert(interactions)
for ppi in interactions:
    source_id = ppi['source'].id_
    target_id = ppi['target'].id_
    row = (source_id, target_id, ppi['pubmed'], ppi['method'], ppi['sources'])
    tb_ppi_rows.append(row)

# write tb_ppi_rows
with TableWriter('tb_ppi') as writer:
    writer.writerow(('source', 'taget', 'pubmed', 'method', 'sources'))
    writer.writerows(tb_ppi_rows)

################################################################################
############# ppi
