import os
import csv
import logging
import itertools

import data


class MirCat(object):
    
    def __init__(self, directory=None):
        """
        http://www.mirrna.org/
        """
        if not directory:
            directory = data.source_data_dir('mircat')
        self.directory = directory
    
    def read_file(self, name):
        path = os.path.join(self.directory, name)
        with open(path) as read_file:
            reader = csv.DictReader(read_file, delimiter='\t')
            for row in reader:
                yield row

    def read_mirna(self):
        name_to_mircat = dict()
        identifier_to_genes = data.Data().mirbase.identifier_to_genes
        matched = 0
        for row in self.read_file('mirna.txt'):
            mircat_name = row['mircat_name']
            hgncs = identifier_to_genes(mircat_name)
            if hgncs is None:
                replaced = mircat_name.replace('-miR-', '-mir-')
                hgncs = identifier_to_genes(replaced)
            if hgncs:
                matched += 1
            row['hgncs'] = hgncs
            name_to_mircat[mircat_name] = row
        logging.info('MirCat miRNAs matched to HGNC {} out of {}'.format(
            matched, len(name_to_mircat)))
        return name_to_mircat
    
    def read_targets(self):
        symbol_to_gene = data.Data().hgnc.get_symbol_to_gene()
        name_to_mircat = self.read_mirna()
        mirna_targets = list()
        for row in self.read_file('target.txt'):
            target = symbol_to_gene.get(row['gene_symbol'])
            mircat_row = name_to_mircat[row['mircat_name']]
            sources = mircat_row['hgncs']
            row['target'] = target
            row['sources'] = sources
            
            if not target or not sources:
                continue
            yield row
    
    def interaction_generator(self):
        for row in self.read_targets():
            target = row['target']
            for source in row['sources']:
                interaction = {'source': source, 'target': target, 'pubmed': row['pubmed']}
                yield interaction
        
        
if __name__ == '__main__':
    mircat = MirCat()
    #mircat.read_mirna()
    #mircat.read_targets()
    import pprint
    ixns = list(mircat.interaction_generator())
    pprint.pprint(ixns)
