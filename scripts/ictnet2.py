import csv
import os

import bioparser.data

ictnet_dir = '/home/dhimmels/Documents/serg/ictnet/ictnet-creation2/'


"""
Do we want approved symbol included in tb_gene_alias?
Are we excluding non-protein coding genes?
"""
class TableWriter(object):
    
    def __init__(self, table_name, fieldnames):
        """Class to write tab delimited text files encoding ictnet tables."""
        self.table_name = table_name
        self.fieldnames = fieldnames
        self.path = os.path.join(ictnet_dir, 'tables', table_name + '.txt')
        self.file = open(self.path, 'w')
        self.writer = csv.DictWriter(self.file, delimiter='\t', fieldnames=fieldnames)
        self.writer.writeheader()
    
    def writerow(self, row):
        self.writer.writerow(row)
    
    def writerows(self, rows, sort=True):
        if sort:
            rows = list(rows)
            rows.sort(key=lambda row: [row[key] for key in self.fieldnames])
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
tb_gene_fieldnames = ('gene_id', 'symbol', 'name', 'chromosome')
tb_gene_rows = list()
tb_gene_alias_fieldnames = ('gene_id', 'alias')
tb_gene_alias_rows = list()
genes = hgnc.get_genes()
for gene in genes:
    gene_id = gene.id_
    row = {'gene_id':gene_id, 'symbol': gene.symbol,
           'name': gene.name, 'chromosome': gene.chromosome}
    tb_gene_rows.append(row)
    for alias in gene.synonyms + gene.previous_symbols:
        row = {'gene_id': gene_id, 'alias': alias}
        tb_gene_alias_rows.append(row)

# write tb_gene_rows
with TableWriter('tb_gene', tb_gene_fieldnames) as writer:
    writer.writerows(tb_gene_rows)

# write tb_gene_alias_rows
with TableWriter('tb_gene_alias', tb_gene_alias_fieldnames) as writer:
    writer.writerows(tb_gene_alias_rows)

################################################################################
############# ppi
tb_ppi_rows = list() #(source, taget, pubmed, method, sources, complex)
tb_ppi_fieldnames = ('source', 'taget', 'pubmed', 'method', 'interaction_type',
                     'edge_type', 'sources', 'complex')

ppitrim = bioparser.data.Data().ppitrim
complex_interactions, binary_interactions = ppitrim.all_interactions()
for is_complex, ppis in enumerate([binary_interactions, complex_interactions]):
    for ppi in ppis:
        ppi['source'] = ppi['source'].id_
        ppi['target'] = ppi['target'].id_
        ppi['complex'] = is_complex
        tb_ppi_rows.append(row)

# write tb_ppi_rows
with TableWriter('tb_ppi') as writer:
    writer.writerow(('source', 'taget', 'pubmed', 'method', 'sources', 'complex'))
    writer.writerows(tb_ppi_rows)

################################################################################
############# ppi
