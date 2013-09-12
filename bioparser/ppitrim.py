import os
import csv
import gzip
import itertools
import collections

import data

class PPITrim(object):
    
    def __init__(self, consolidated_path=None):
        """
        website: http://www.ncbi.nlm.nih.gov/CBBresearch/Yu/downloads/ppiTrim.html
        dataset download: ftp://ftp.ncbi.nlm.nih.gov/pub/qmbpmn/ppiTrim/datasets/
        publication: 10.1093/database/bar036
        """
        if not consolidated_path:
            consolidated_path = os.path.join(data.source_data_dir('iref'),
                'release_12.0', '9606.mitab.06062013.txt.ppiTrim.txt.gz')
        self.consolidated_path = consolidated_path
    
    def row_generator(self):
        """
        
        http://database.oxfordjournals.org/content/suppl/2011/07/14/bar036.DC1/ppidb-supp.pdf
        """
        read_file = gzip.open(self.consolidated_path)
        fieldnames = read_file.readline().rstrip().lstrip('#').split('\t')
        if len(fieldnames) == 35:
            fieldnames.append('negative')
        reader = csv.DictReader(read_file, delimiter='\t', fieldnames=fieldnames)
        plural_fields = ['uidA', 'uidB', 'confidence']
        multiplex_fields = ['altA', 'altB', 'aliasA', 'aliasB', 'pmids', 'interactionIdentifier']
        for row in reader:
                        
            # split all fields into elements
            for fieldname, field in row.items():
                if field != '-':
                    row[fieldname] = field.split('|')
                        
            # fields where each element contains a unique id
            for fieldname in plural_fields:
                field = row[fieldname]
                if field == '-':
                    row[fieldname] = dict()
                    continue
                field_dict = dict()
                for elem in field:
                    key, value = elem.split(':')
                    field_dict[key] = value
                row[fieldname] = field_dict

            # fields with multiple elements sharing the same key
            for fieldname in multiplex_fields:
                field = row[fieldname]
                if field == '-':
                    row[fieldname] = dict()
                    continue
                field_dict = dict()
                for elem in field:
                    if elem == '-':
                        continue
                    key, value = elem.split(':', 1) # go:GO:0005783 requires split(':', 1)
                    field_dict.setdefault(key, list()).append(value)
                row[fieldname] = field_dict
            
            yield row
        read_file.close()
    
    def separate_complexes(self):
        complex_to_rows = dict()
        binary_rows = list()
        for row in self.row_generator():
            complex = row['uidA'].get('complex')
            if complex:
                complex_to_rows.setdefault(complex, list()).append(row)
            else:
                binary_rows.append(row)
        return complex_to_rows, binary_rows
    
    
    def get_interaction_info(self, row):
        interaction_info = collections.OrderedDict()
        interaction_info['pubmed'] = ', '.join(row['pmids'].get('pubmed', ''))
        interaction_info['method'] = ', '.join(row['method'])
        interaction_info['interaction_type'] = ', '.join(row['interactionType'])
        interaction_info['edge_type'] = ', '.join(row['interactionIdentifier'].get('edgetype', list()))
        interaction_info['sources'] = ', '.join(row['sourcedb'])
        # all binary interactions are 'none' and complexes are 'bipartite' for expansion
        #interaction_info['expansion'] = ', '.join(row['expansion'])
        return interaction_info

    def binary_rows_to_interactions(self, binary_rows):
        for row in binary_rows:
            interaction_info = self.get_interaction_info(row)
            aliases_a = row['aliasA']['entrezgene/locuslink']
            aliases_b = row['aliasB']['entrezgene/locuslink']
            symbol_pairs = itertools.product(aliases_a, aliases_b)
            for symbol_a, symbol_b in symbol_pairs:
                interaction = collections.OrderedDict()
                interaction['source'] = symbol_a
                interaction['target'] = symbol_b
                interaction.update(interaction_info)
                yield interaction
    
    def complexes_to_interactions(self, complex_to_rows):
        """ """
        for complex, rows in complex_to_rows.iteritems():
            interaction_info = self.get_interaction_info(rows[0])
            aliases = set()
            for row in rows:
                aliases_b = row['aliasB']['entrezgene/locuslink']
                aliases |= set(aliases_b)
            symbol_pairs = itertools.combinations(aliases, 2) # matrix model
            for symbol_a, symbol_b in symbol_pairs:
                #interaction_info_temp = get_interaction_info(row)
                #assert interaction_info_temp.values() == interaction_info.values()
                interaction = collections.OrderedDict()
                interaction['source'] = symbol_a
                interaction['target'] = symbol_b
                interaction.update(interaction_info)
                yield interaction

    def interaction_hgnc_convert(self, interactions, exclude_self_references=True):
        hgnc = data.Data().hgnc
        symbol_to_gene = hgnc.get_symbol_to_gene()
        for interaction in interactions:
            source = symbol_to_gene.get(interaction['source'])
            target = symbol_to_gene.get(interaction['target'])
            if not source or not target:
                continue
            if exclude_self_references and source == target:
                continue
            interaction = interaction.copy()
            interaction['source'] = source
            interaction['target'] = target
            yield interaction
            
    def all_interactions(self):
        complex_to_rows, binary_rows = self.separate_complexes()

        complex_interactions = self.complexes_to_interactions(complex_to_rows)
        complex_interactions = self.interaction_hgnc_convert(complex_interactions)
        
        binary_interactions = self.binary_rows_to_interactions(binary_rows)
        binary_interactions = self.interaction_hgnc_convert(binary_interactions)
        
        complex_interactions = list(complex_interactions)
        binary_interactions = list(binary_interactions)
        
        return complex_interactions, binary_interactions

    
if __name__ =='__main__':
    ppi = PPITrim()
    complex_interactions, binary_interactions = ppi.all_interactions()
    print 'complex', len(complex_interactions)
    print 'binary', len(binary_interactions)