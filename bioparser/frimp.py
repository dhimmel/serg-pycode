import os
import csv
import logging
import gzip

import data

class IMP(object):
    
    def __init__(self, imp_dir=None):
        """ """
        if imp_dir is None:
            imp_dir = data.current_path('imp')

        self.imp_dir = imp_dir

    
    def entrez_to_symbols(self, *args):
        """ """
        if not hasattr(self, 'entrez_to_hgnc'):
            self.entrez_to_hgnc = data.Data().hgnc.get_entrez_to_gene()
        if not hasattr(self, 'unmatched_entrez'):
            self.unmatched_entrez = set()

        symbols = list()
        for entrez in args:
            gene = self.entrez_to_hgnc.get(entrez)
            if gene is None:
                self.unmatched_entrez.add(entrez)
                return None
            symbol = gene.symbol
            symbols.append(symbol)
        
        if len(set(symbols)) < len(args):
            #logging.warning('Functional self relations: ' + str(args) + 'translated to' +  str(symbols))
            return None
        
        return symbols

    def read_dat(self, path):
        """ """
        dat_file = gzip.open(path)
        reader = csv.reader(dat_file, delimiter='\t')
        for row in reader:
            row[2] = float(row[2])
            row = row[:3] # in case extra columns or trailings tabs included
            row = tuple(row)
            yield row
        dat_file.close()

    def write_dat(self, path, relationships):
        """ """
        with gzip.open(path, 'w') as wf:
            writer = csv.writer(wf, delimiter='\t')
            writer.writerows(relationships)

    def read_dat_subset(self, path, prob_cutoff = 0.0, symbol=False):
        """ """
        rows = self.read_dat(path)
        for gene_0, gene_1, prob in rows:
            if prob < prob_cutoff:
                continue
            if symbol:
                symbols = self.entrez_to_symbols(gene_0, gene_1)
                if symbols is None:
                    continue
                gene_0, gene_1 = symbols       
            yield gene_0, gene_1, prob

    def read_predictions(self, prob_cutoff = 0.0, symbol=False):
        """ """
        dat_path = os.path.join(self.imp_dir, 'originals', 'global_average_prior.dat.gz')
        generator = self.read_dat_subset(dat_path, prob_cutoff=prob_cutoff, symbol=symbol)
        return generator

    def read_gold_standard_positives(self, symbol=False):
        """ """
        dat_path = os.path.join(self.imp_dir, 'originals', 'human_positives.dat.gz')
        generator = self.read_dat_subset(dat_path, prob_cutoff=0.5, symbol=symbol)
        return generator
    
    def get_relationships(self, predictions=True, positives=True, prob_cutoff=0.5, symbol=True):
        """
        Returns a list of tuples which each represent a functional
        relationship between two genes. The tuples are formatted 
        (gene_0, gene_1, probabiltiy). Predictions of functional relationship
        status are included when predictions=True and are filtered with the
        prob_cutoff argument. Gold standard positives are included when
        positives=True. To convert entrez ids to gene symbol use symbol=True. If
        symbol=True, unmatched symbols are written to the log file and pairs
        with an unmatched symbol are omitted.
        """
        messages = list()
        messages.append('IMP functional relationships')
        messages.append('Include gold standard positives: ' + str(positives))
        messages.append('Include predictions: ' + str(positives) + '. Filter probabilities below ' + str(prob_cutoff))
        messages.append('Convert to symbols: ' + str(symbol))
        message = '\n'.join(messages)
        logging.info(message)
        genes_to_prob = dict()
        
        if predictions:
            fr_tuple = self.read_predictions(prob_cutoff, symbol)
            for gene_0, gene_1, prob in fr_tuple:
                genes = frozenset([gene_0, gene_1])
                genes_to_prob[genes] = prob

        if positives:
            fr_tuple = self.read_gold_standard_positives(symbol)
            for gene_0, gene_1, prob in fr_tuple:
                genes = frozenset([gene_0, gene_1])
                genes_to_prob[genes] = prob

        if symbol:
            for entrez in self.unmatched_entrez:
                logging.warning('No symbol found for IMP entrez ' + entrez)
        
        relation_tuples = [tuple(sorted(genes) + [prob]) 
                           for genes, prob in genes_to_prob.items()]
        relation_tuples.sort()
        return relation_tuples

    def process(self):
        processed_dir = os.path.join(self.imp_dir, 'processed')

        path = os.path.join(processed_dir, 'positives.txt.gz')
        relationships = self.get_relationships(predictions=False, positives=True, symbol=True)
        self.write_dat(path, relationships)

        prob_cutoffs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        for prob_cutoff in prob_cutoffs:
            path = os.path.join(processed_dir, 'predictions_%s.txt.gz' % prob_cutoff)
            relationships = self.get_relationships(predictions=True, positives=False, prob_cutoff=prob_cutoff, symbol=True)
            self.write_dat(path, relationships)

            path = os.path.join(processed_dir, 'positives_and_predictions_%s.txt.gz' % prob_cutoff)
            relationships = self.get_relationships(predictions=True, positives=True, prob_cutoff=prob_cutoff, symbol=True)
            self.write_dat(path, relationships)
    
    def read_processed_relationships(self, name):
        path = os.path.join(self.imp_dir, 'processed', name + '.txt.gz')
        relationships = self.read_dat(path)
        return relationships

        
if __name__ == '__main__':
    
    imp = IMP()
    imp.process()
    #imp.get_pairs(predictions=True, positives=False)
    