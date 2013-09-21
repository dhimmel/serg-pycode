import collections
import csv
import os

import utilities.omictools
import data
import utilities.shelved


class Concept(object):

    def __init__(self, concept_id):
        """A Concept object represents a UMLS metathesaurus concept."""
        self.concept_id = concept_id
        self.semantic_types = set()
        self.name = None
        self.source_to_rank = dict()
        self.source_to_code = dict()
        self.source_to_name = dict()
        self.name_rank = 0
    
    def __hash__(self):
        return hash(self.concept_id)

    def __eq__(self, other):
        return self.concept_id == other.concept_id
    
    def __str__(self):
        s = 'concept_id: ' + self.concept_id
        for attribute in ['name', 'semantic_types']:
            if not getattr(self, attribute):
                continue
            s += '\n  %s: %s' % (attribute, getattr(self, attribute))
        if self.source_to_code:
            s += '\n  Source Vocabularies'
            for key, value in self.source_to_code.items():
                s += '\n    %s: %s' % (key, value)
        return s

class Metathesaurus(utilities.shelved.Shelved):

    def __init__(self, umls_dir=None):
        """
        Sources: http://www.nlm.nih.gov/research/umls/sourcereleasedocs/index.html#
        """
        if not umls_dir:
            umls_dir = os.path.join(data.data_dir, 'umls', '2013AA')
        self.umls_dir = umls_dir
        self.meta_dir = os.path.join(self.umls_dir, 'META')
        
        self.pymeta_dir = os.path.join(self.umls_dir, 'pyMETA')
        if not os.path.exists(self.pymeta_dir):
            os.mkdir(self.pymeta_dir)
        
        shelve_names = ['sources', 'concepts', 'types']
        super(Metathesaurus, self).__init__(self.pymeta_dir, shelve_names)

        if not self.shelve_files_exist():
            self.extract_source_vocabs = {'MSH', 'NCI', 'OMIM', 'MDR'}
            self.build()
        
    def build(self, printing=True):
        """
        Reference Manual: http://www.ncbi.nlm.nih.gov/books/NBK9685/
        Abbreviations: http://www.nlm.nih.gov/research/umls/knowledge_sources/metathesaurus/release/abbreviations.html
        Column descriptions: http://www.nlm.nih.gov/research/umls/knowledge_sources/metathesaurus/release/columns_data_elements.html
        """
        if printing:
            print 'Building UMLS metathesaurus.'
            print 'Extracting:', ', '.join(self.extract_source_vocabs)
        self.fname_to_fieldnames = dict()
        id_to_concept = dict()
        
        # Gather concept_ids and their semantic types.
        for row in self.read_rrf('MRSTY.RRF'):
            # Each row represents a concept to semantic_type pair
            semantic_type = row['STY']
            concept_id = row['CUI']
            concept = id_to_concept.setdefault(concept_id, Concept(concept_id))
            concept.semantic_types.add(semantic_type)
        
        source_and_tty_to_rank = dict()
        source_to_tty_ranks = dict()
        for row in self.read_rrf('MRRANK.RRF'):
            rank = row['RANK']
            source = row['SAB']
            term_type = row['TTY']
            source_and_tty_to_rank[(source, term_type)] = rank
            source_to_tty_ranks.setdefault(source, list()).append(term_type)
            
        
        # Iterate through each atom (each occurrence of each unique string 
        # or concept name within each source vocabulary.
        for row in self.read_rrf('MRCONSO.RRF'):
            
            concept_id = row['CUI']
            concept = id_to_concept.get(concept_id)
            if concept is None:
                continue
            
            name = row['STR']
            term_type = row['TTY']
            name = name.decode('utf-8').encode('ascii', 'replace')
            # Extract source name and source vocabulary.
            source = row['SAB']
            row_rank = source_and_tty_to_rank[(source, term_type)]
            if (self.extract_source_vocabs is None or 
                source in self.extract_source_vocabs):
                
                source_and_tty_to_rank[(source, term_type)]
                current_rank = concept.source_to_rank.setdefault(source, 0)
                #if term_type == source_to_tty_ranks[source][0]:                
                if row_rank > current_rank:
                    concept.source_to_rank[source] = row_rank
                    concept.source_to_code[source] = row['CODE']
                    concept.source_to_name[source] = name
            
            if row_rank > concept.name_rank:
                concept.name = name
                concept.name_rank = row_rank
        
        # check all concepts have names.
        for concept_id, concept in id_to_concept.items():
            if not concept.name:
                print 'no concept name\n', concept
                
        # Shelve creation
        type_to_ids, source_to_ids = dict(), dict()
        for concept_id, concept in id_to_concept.iteritems():
            del concept.name_rank
            del concept.source_to_rank
            for semantic_type in concept.semantic_types:
                type_to_ids.setdefault(semantic_type, set()).add(concept_id)
            for source in concept.source_to_code:
                source_to_ids.setdefault(source, set()).add(concept_id)
        
        if printing:
            print 'Writing metathesaurus to shelves.'
        with self:
            self.shelves['concepts'].update(id_to_concept)
            self.shelves['types'].update(type_to_ids)
            self.shelves['sources'].update(source_to_ids)
        if printing:
            print 'Metathesaurus building complete.'
        
    def fieldnames(self, fname):
        """Get the fieldnames for a rrf file from MRFILES.RRF.
        Source: http://www.ncbi.nlm.nih.gov/books/NBK9685/
        """
        if not self.fname_to_fieldnames:
            fieldnames = ['FIL', 'DES', 'FMT', 'CLS', 'RWS', 'BTS']
            rows = self.read_rrf('MRFILES.RRF', fieldnames=fieldnames)
            for row in rows:
                row['FMT'] = row['FMT'].split(',')
                self.fname_to_fieldnames[row['FIL']] = row['FMT']
        
        return self.fname_to_fieldnames[fname]
    
    def read_rrf(self, fname, fieldnames=None):
        """
        Read an rrf file returning a generator of dictionaries. Each
        dictionary represents a row with fieldnames as keys and row values
        as values.
        """
        if not fieldnames:
            fieldnames = self.fieldnames(fname)
        path = os.path.join(self.meta_dir, fname)
        with open(path) as f:
            reader = csv.DictReader(f, fieldnames=fieldnames, delimiter='|')
            for row in reader:
                yield row
    
    def get_source_code_to_concept_id(self, source):
        """ """
        with self:
            id_to_concept = self.shelves['concepts']
            concept_ids = self.shelves['sources'][source]
            source_code_to_concept_id = {id_to_concept[cid].source_to_code[source]: cid for cid in concept_ids}
        return source_code_to_concept_id
    
    def get_source_code_to_concept(self, source):
        """ """
        source_code_to_concept = dict()
        with self:
            id_to_concept = self.shelves['concepts']
            concept_ids = self.shelves['sources'][source]
            for cid in concept_ids:
                concept = id_to_concept[cid]
                source_code = concept.source_to_code[source]
                source_code_to_concept[source_code] = concept
        return source_code_to_concept

if __name__ == '__main__':
    meta = Metathesaurus()
    with meta:
        concepts = meta.shelves['concepts']
        print concepts['C0155555']
  
    