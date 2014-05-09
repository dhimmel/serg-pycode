import os
import csv
import gzip
import collections

import data
import utilities.omictools

class Gene(object):
    
    def __init__(self, symbol):
        self.symbol = symbol
    
    def __hash__(self):
        return hash(self.symbol)
    
    def __eq__(self, other):
        return self.symbol == other.symbol
    
    def __str__(self):
        return self.symbol

    def __repr__(self):
        return self.symbol    

class HGNC(object):
    
    def __init__(self, directory=None):
        if directory is None:
            directory = data.current_path('hgnc')
        self.directory = directory
        self.hgnc_path = os.path.join(directory, 'hgnc_complete_set.txt.gz')
        self.genes = None
        self.symbol_to_gene = dict()
        self.entrez_to_gene = dict()
        self.ensembl_to_gene = dict()
        #self.prop_to_dict = dict()

    def gene_generator(self):
        """
        Produce a generator where each item is a gene.
        
        Website: http://www.genenames.org/cgi-bin/hgnc_stats
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
                           'UniProt ID (supplied by UniProt)': 'uniprot_id',
                           'Entrez Gene ID (supplied by NCBI)': 'entrez_id_ncbi_mapped',
                           'Ensembl ID (supplied by Ensembl)': 'ensembl_id_ensembl_mapped',
                           'OMIM ID (supplied by NCBI)': 'omim_id'}
        keys_requiring_splitting = ['previous_symbols', 'synonyms', 'refseq_ids']
        read_file = gzip.open(self.hgnc_path)
        reader = csv.DictReader(read_file, delimiter='\t')
        for row in reader:
            if row['Status'] != 'Approved':
                continue
            row = {converted_key: row[key] for key, converted_key in key_conversions.items()}
            for key, value in row.items():
                if value == '':
                    row[key] = None
            for key in keys_requiring_splitting:
                value = row[key]
                if value is None:
                    value = list()
                else:
                    value = value.split(', ')
                row[key] = value
            gene = Gene(row['symbol'])
            gene.__dict__.update(row)
            gene.int_id = int(gene.hgnc_id.split(':')[1])
            gene.coding = gene.locus_group == 'protein-coding gene'
            yield gene
        read_file.close()
    
    def get_genes(self):
        """Return a list of genes."""
        if self.genes is None:
            gene_generator = self.gene_generator()
            self.genes = list(set(gene_generator))
        return self.genes



    def get_symbol_to_gene(self):
        """
        The approved symbol for a gene points to that gene. Previous symbols point to their
        gene if the previous symbol is not an approved symbol and symbol is not a previous
        symbol for multiple genes. Synonyms point to their gene if the synonym is not an approved
        or previous symbol and if the symbol is not a synonym for multiple genes.
        """
        if self.symbol_to_gene:
            return self.symbol_to_gene
        genes = self.get_genes()

        approved_symbol_to_gene = dict()
        previous_symbol_to_genes = dict()
        synonyms_to_genes = dict()
        for gene in genes:
            approved_symbol_to_gene[gene.symbol] = gene
            for previous_symbol in gene.previous_symbols:
                previous_symbol_to_genes.setdefault(previous_symbol, set()).add(gene)
            for synonym in gene.synonyms:
                synonyms_to_genes.setdefault(synonym, set()).add(gene)

        symbol_to_gene = approved_symbol_to_gene.copy()

        previous_symbols = set(previous_symbol_to_genes) - set(approved_symbol_to_gene)
        for previous_symbol in previous_symbols:
            genes = previous_symbol_to_genes[previous_symbol]
            try:
                gene, = genes
                symbol_to_gene[previous_symbol] = gene
            except ValueError:
                pass

        synonyms = set(synonyms_to_genes) - (set(approved_symbol_to_gene) | set(previous_symbol_to_genes))
        for synonym in synonyms:
            genes = synonyms_to_genes[synonym]
            try:
                gene, = genes
                symbol_to_gene[synonym] = gene
            except ValueError:
                pass

        self.symbol_to_gene = symbol_to_gene
        return self.symbol_to_gene

    def get_entrez_to_gene(self):
        """ """
        if self.entrez_to_gene:
            return self.entrez_to_gene
        genes = self.get_genes()
        self.entrez_to_gene = {gene.entrez_id: gene for gene in genes}
        self.entrez_to_gene.update({gene.entrez_id_ncbi_mapped: gene for gene in genes})
        return self.entrez_to_gene        

    def get_ensembl_to_gene(self):
        """ """
        if self.ensembl_to_gene:
            return self.ensembl_to_gene
        genes = self.get_genes()
        self.ensembl_to_gene = {gene.ensembl_id: gene for gene in genes}
        self.ensembl_to_gene.update({gene.ensembl_id_ensembl_mapped: gene for gene in genes})
        return self.ensembl_to_gene        

    def identifiers_to_genes(self, identifiers, id_type='symbol', coding=False):
        if id_type == 'symbol':
            id_to_gene = self.get_symbol_to_gene()
        elif id_type == 'entrez':
            id_to_gene = self.get_entrez_to_gene()
        elif id_type == 'ensembl':
            id_to_gene = self.get_ensembl_to_gene()
        else:
            raise ValueError('id_type unkown')
        if coding:
            id_to_gene = {id_: gene for id_, gene in id_to_gene.iteritems()
                          if gene.locus_group == 'protein-coding gene'}
        genes = {id_to_gene.get(identifier) for identifier in identifiers}
        genes.discard(None)
        return genes

    def write_as_table(self):
        """Write as tab delimited table of protein coding genes.
        """
        genes = self.get_genes()
        genes = [gene for gene in genes if gene.locus_group == 'protein-coding gene']
        rows = list()
        for gene in genes:
            row = collections.OrderedDict()
            rows.append(row)
            row['hgnc_id'] = gene.hgnc_id
            row['symbol'] = gene.symbol
            row['name'] = gene.name
            row['aliases'] = '|'.join(sorted(set(gene.previous_symbols + gene.synonyms)))
            row['entrez'] = gene.entrez_id or gene.entrez_id_ncbi_mapped
            row['ensembl'] = gene.ensembl_id or gene.ensembl_id_ensembl_mapped
            row['uniprot'] = gene.uniprot_id
            row['chromosome'] = gene.chromosome
        rows.sort(key=lambda x: x['symbol'])
        path = os.path.join(self.directory, 'protein-coding.txt')
        write_file = open(path, 'w')
        writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
        write_file.close()


if __name__ == '__main__':
    hgnc = HGNC()
    #hgnc.write_as_table()
    entrez_to_gene = hgnc.get_entrez_to_gene()
    print entrez_to_gene['']
    #gene_number = len(hgnc.get_genes())
    #symbol_number = len(hgnc.get_symbol_to_gene())
    #print symbol_number, 'symbols representing', gene_number, 'genes'
    