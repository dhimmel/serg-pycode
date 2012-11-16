import os

import data
import utilities.omictools

class Gene(object):
    
    def __init__(self, symbol):
        self.symbol = symbol
    
    def __hash__(self):
        return hash(self.symbol)
    
    def __eq__(self, other):
        return self.symbol == other.symbol
    
class HGNC(object):
    
    def __init__(self, hgnc_dir=None):
        if hgnc_dir is None:
            hgnc_dir = data.current_path('hgnc')
        self.hgnc_dir = hgnc_dir
        self.hgnc_path = os.path.join(hgnc_dir, 'hgnc-complete.txt')
        self.genes = None
        self.symbol_to_gene = dict()
        self.entrez_to_gene = dict()
        #self.prop_to_dict = dict()
        
    def gene_generator(self):
        """
        Produce a generator where each item is a gene.
        
        Website: http://www.genenames.org/cgi-bin/hgnc_stats.pl
        Download the "complete HGNC dataset"
        
        Genes with Status 'Entry Withdrawn' or 'Symbol Withdrawn' are excluded.
        """
        key_conversions = {'HGNC ID': 'hgnc_id',
                           'Approved Symbol': 'symbol',
                           'Approved Name': 'name',
                           'Locus Type': 'locus_type',
                           'Locus Group': 'locus_group',
                           'Previous Symbols': 'previous_symbols',
                           'Synonyms': 'synonyms',
                           'Chromosome': 'chromosome',
                           'Entrez Gene ID': 'entrez_id',
                           'Ensembl Gene ID': 'ensembl_id',
                           'RefSeq IDs': 'refseq_ids', 
                           'OMIM ID (mapped data supplied by NCBI)': 'omim_id'}
        keys_requiring_splitting = ['previous_symbols', 'synonyms', 'refseq_ids']
        for line_dict in utilities.omictools.read_tdt(self.hgnc_path):
            if line_dict['Status'] != 'Approved':
                continue
            line_dict = {converted_key: line_dict[key] for key, converted_key in key_conversions.items()}
            for key, value in line_dict.items():
                if value == '':
                    line_dict[key] = None
            for key in keys_requiring_splitting:
                value = line_dict[key]
                if value is None:
                    value = list()
                else:
                    value = value.split(', ')
                line_dict[key] = value
            gene = Gene(line_dict['symbol'])
            gene.__dict__.update(line_dict)
            yield gene
    
    def get_genes(self):
        """Return a list of genes. Each gene is a dictionary describing the gene."""
        if self.genes is None:
            gene_generator = self.gene_generator()
            self.genes = list(gene_generator)
        return self.genes
    
    def get_symbol_to_gene(self):
        if self.symbol_to_gene:
            return self.symbol_to_gene
        genes = self.get_genes()
        for gene in genes:
            keys = [gene.symbol]
            keys.extend(gene.synonyms)
            keys.extend(gene.previous_symbols)
            #keys.append(gene.name)
            self.symbol_to_gene.update(dict.fromkeys(keys, gene))
        return self.symbol_to_gene

    def get_entrez_to_gene(self):
        if self.entrez_to_gene:
            return self.entrez_to_gene
        genes = self.get_genes()
        entrez_to_gene = {gene.entrez_id: gene for gene in genes}
        return self.entrez_to_gene        
        
if __name__ == '__main__':
    hugu = HGNC()
    gene_number = len(hugu.get_genes())
    symbol_number = len(hugu.get_symbol_to_gene())
    print symbol_number, 'symbols representing', gene_number, 'genes'
    