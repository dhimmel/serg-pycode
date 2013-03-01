import csv
import os

import mapping.metamap

import data


class DvD(object):

    def __init__(self, dvd_dir=None):
        if dvd_dir is None:
            dvd_dir = data.current_path('dvd', require_dated_format=True)

        self.dvd_dir = dvd_dir
    
    def read_expression_file(self, name):
        """ """
        path = os.path.join(self.dvd_dir, name)
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                yield row
    
    def get_drug_expression(self):
        """Returns tuples a tuple describing an instance of a differently
        expressed gene for a specific drug."""
        hgnc = data.Data().hgnc
        symbol_to_hgnc = hgnc.get_symbol_to_gene()

        rows = self.read_expression_file('drug-expression.txt')
        drug_to_concept_id = self.get_drug_to_concept_id()
        
        for row in rows:
            concept_id = drug_to_concept_id.get(row['condition'])
            gene = symbol_to_hgnc.get(row['gene'])
            if not concept_id or not gene:
                continue
            gene = gene.symbol
            yield concept_id, gene, row['direction']
            
    def get_disease_expression(self):
        """Returns tuples a tuple describing an instance of a differently
        expressed gene for a specific disease."""
        hgnc = data.Data().hgnc
        symbol_to_hgnc = hgnc.get_symbol_to_gene()
        
        rows = self.read_expression_file('disease-expression.txt')
        disease_to_concept_id = self.get_disease_to_concept_id()

        for row in rows:
            concept_id = disease_to_concept_id.get(row['condition'])
            gene = symbol_to_hgnc.get(row['gene'])
            if not concept_id or not gene:
                continue
            gene = gene.symbol
            yield concept_id, gene, row['direction']
    
    def read_each_line(self, fname):
        path = os.path.join(self.dvd_dir, fname)
        with open(path) as f:
            list_of_lines = [line.rstrip() for line in f]
        return list_of_lines

    def get_disease_to_concept_id(self):
        diseases = self.read_each_line('disease-names.txt')
        map_dir = os.path.join(self.dvd_dir, 'disease-mappings')
        mapper = mapping.metamap.MetaMap(map_dir, restrict_to_sources=['NCI'])
        
        if not mapper.maps['current']['redacted']:
            mapper.map_terms(diseases)
            mapper.pause_before_reading_mappings()
        
        return mapper.get_mapping_dict()

    def get_drug_to_concept_id(self):
        drugs = self.read_each_line('drug-names.txt')
        
        map_dir = os.path.join(self.dvd_dir, 'drug-mappings')
        mapper = mapping.metamap.MetaMap(map_dir, restrict_to_sources=['MSH'], restrict_to_types=['phsu'])
        
        if not mapper.maps['current']['redacted']:
            mapper.map_terms(drugs)
            mapper.pause_before_reading_mappings()
        
        return mapper.get_mapping_dict()


if __name__ == '__main__':
    dvd = DvD()
    print dvd.metamap_drugs()