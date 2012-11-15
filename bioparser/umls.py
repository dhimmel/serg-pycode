import ast
import collections
import csv
import os
import shelve

import omictools

class Concept(object):

    def __init__(self, concept_id):
        """A Concept object represents a UMLS metathesaurus concept."""
        self.concept_id = concept_id
        self.symantic_types = set()
        self.name = None
        self.vocabs = dict()
    
    def __hash__(self):
        return hash(self.concept_id)
    
    def __str__(self):
        s = 'concept_id: ' + self.concept_id
        for attribute in ['name', 'symantic_types']:
            if not getattr(self, attribute):
                continue
            s += '\n  %s: %s' % (attribute, getattr(self, attribute))
        if self.vocabs:
            s += '\n  Source Vocabularies'
            for key, value in self.vocabs.items():
                s += '\n    %s: %s' % (key, value)
        return s


class Metathesaurus(object):
    
    def __init__(self, meta_dir=None):
        """
        A Metathesaurus object parses and stores UMLS Metathesaurus concepts.
        
        Parsing is performed by the read method. The meta_dir parameter
        specifies the location of the UMLS Metathesaurus rrf files which
        serve as the data source to parse. The set_profile method 
        
        Storage exists in two object attributes:
        self.concepts      - a set of Concept objects
        self.id_to_concept - a concept_id to Concept object dictionary
        
        # http://www.nlm.nih.gov/research/umls/META3_current_semantic_types.html
        # Source vocab info: http://www.nlm.nih.gov/research/umls/knowledge_sources/metathesaurus/release/source_vocabularies.html
        # RRF description: http://www.ncbi.nlm.nih.gov/books/NBK9685/
        """

        if not meta_dir:
            meta_dir = os.path.join(omictools.data_dir, 'umls', '2012AA', 'META')
            
        self.meta_dir = meta_dir
        self.fname_to_fieldnames = dict()
                
        self.profile_attributes = ['keep_symantic_types',
                                   'extract_source_vocabs',
                                   'require_one_source',
                                   'require_all_sources']
        self.set_profile('empty')
        
    def set_profile(self, profile):
        """A profile directs the parsing and filtering while reading the
        metathesaurus. profile is either a string defining a predefined profile
        defined in this function or a dictionary with keys matching
        self.profile_attributes and user set values. Values of None for
        profile_attributes specify that no restriction or filtering
        should occur based on that attribute.
        """
                
        if isinstance(profile, dict):
            
            assert set(profile) == set(self.profile_attributes)
            for key, value in profile.items():
                if value is not None:
                    value = set(value)
                setattr(self, key, value)
            
        else:
            
            assert profile in ['full', 'all', 'omicnet', 'empty']
            self.profile = profile
            
            if profile == 'empty':
                # No concepts are read
                self.keep_symantic_types = set()
                self.extract_source_vocabs = set()
                self.require_one_source = set()
                self.require_all_sources = set()
    
            if profile == 'full':
                # All concepts and sources are read.
                self.keep_symantic_types = None
                self.extract_source_vocabs = None
                self.require_one_source = None
                self.require_all_sources = None

            if profile == 'all':
                # All concepts are read but no sources are retained.
                self.keep_symantic_types = None
                self.extract_source_vocabs = set()
                self.require_one_source = None
                self.require_all_sources = None

            if profile == 'omicnet':
                self.keep_symantic_types = set(['Disease or Syndrome', 'Neoplastic Process', 'Anatomical Abnormality'])
                self.extract_source_vocabs = set(['ICD10CM', 'ICD9CM', 'MTHICD9', 'SNOMEDCT', 'MSH', 'NCI', 'OMIM'])
                self.require_one_source = set(['MSH'])
                self.require_all_sources = set(['MSH'])
            
    
    def read(self):
        """
        Parse UMLS metathesaurus files and return the resulting set of
        concepts. Parsing and filtering is performed in accordance the
        attributes named in self.profile_attributes. This function populates
        self.concepts and self.id_to_concept.
        """
        
        id_to_concept = dict()
        
        # Gather concept_ids and their symantic types.
        for row in self.read_rrf('MRSTY.RRF'):
            # Each row represents a concept to symantic_type pair
            symantic_type = row['STY']
            if self.keep_symantic_types is None or symantic_type in self.keep_symantic_types:
                concept_id = row['CUI']
                concept = id_to_concept.setdefault(concept_id, Concept(concept_id))
                concept.symantic_types.add(symantic_type)
        
        # Iterate through each atom (each occurrence of each unique string 
        # or concept name within each source vocabulary.
        for row in self.read_rrf('MRCONSO.RRF'):
            
            concept_id = row['CUI']
            concept = id_to_concept.get(concept_id)
            if concept is None:
                continue
            
            # Extract source name and source vocabulary.
            source_vocab = row['SAB']
            if (self.extract_source_vocabs is None or 
                source_vocab in self.extract_source_vocabs):
                concept.vocabs[source_vocab] = row['CODE']
            
            # Check if the source designates the preferred concept name.
            if (row['TS'] == 'P' and
                row['STT'] == 'PF' and
                row['ISPREF'] == 'Y' and 
                row['SUPPRESS'] == 'N' and
                row['LAT'] == 'ENG'):
                # Check that no other preferred concept names have been added.
                if concept.name is not None:
                    raise Exception('Previosely set Concept name: Multiple Preferred Concept Names')
                name = row['STR']
                name = name.decode('utf-8').encode('ascii', 'replace')
                concept.name = name

        # Filter concepts which do not pass source vocabulary constraints.
        for concept_id, concept in id_to_concept.items():
            vocab_names = set(concept.vocabs.keys())
            
            if (self.require_one_source is not None and 
                not self.require_one_source & vocab_names):
                # Delets CUIs which do not have one source matching self.require_one_source
                del id_to_concept[concept_id]
                
            elif (self.require_all_sources is not None and 
                  not self.require_all_sources <= vocab_names):
                # Delets concepts without every concept specified by self.require_all_sources
                del id_to_concept[concept_id]
                
            elif concept.name is None:
                # Delete concepts with no Preferred Cocept Name
                del id_to_concept[concept_id]
                
        self.id_to_concept = id_to_concept
        self.concepts = set(id_to_concept.values())
        return self.concepts
    
    def refresh_id_to_concept(self):
        """Recreate the mapping dict id_to_concept using self.concepts."""
        self.id_to_concept = {concept.concept_id: concept for concept in self.concepts}
    
    def get(self, concept_id):
        """Returns the concept with the specified concept_id. Returns None
        if no concept exists with that id."""
        return self.id_to_concept.get(concept_id)

    def __getitem__(self, concept_id):
        return self.id_to_concept[concept_id]


    def __str__(self):
        """Returns a string summary of the object."""
        
        s = 'UMLS Object from ' + self.meta_dir
        
        s += '\nProfile:'
        for profile_attr in self.profile_attributes:
            s += '\n  %s: %s' % (profile_attr, getattr(self, profile_attr))

        self.compute_counters()
        
        s += '\nNumber of concepts by type:'
        for key, count in sorted(self.symantic_type_counter.items()):
            s += '\n  %s: %s' % (key, count)

        s += '\nNumber of concepts by source:'
        for key, count in sorted(self.vocab_counter.items()):
            s += '\n  %s: %s' % (key, count)

        return s
    
    def compute_counters(self):
        """Sets two counters"""
        self.symantic_type_counter = collections.Counter()
        self.vocab_counter = collections.Counter()
        for concept in self.concepts:
            self.symantic_type_counter.update(concept.symantic_types)
            self.vocab_counter.update(concept.vocabs.keys())

    
    def fieldnames(self, fname):
        """Get the fieldnames for a rrf file from MRFILES.RRF.
        Source: http://www.ncbi.nlm.nih.gov/books/NBK9685/
        """
        fname_to_fieldnames = self.fname_to_fieldnames
        
        if not fname_to_fieldnames:
        
            with open(os.path.join(self.meta_dir, 'MRFILES.RRF')) as f:
                fieldnames = ['FIL', 'DES', 'FMT', 'CLS', 'RWS', 'BTS']
                dict_reader = csv.DictReader(f, fieldnames=fieldnames, delimiter='|')
                rows = list(dict_reader)
                
            for row in rows:
                row['FMT'] = row['FMT'].split(',')
                fname_to_fieldnames[row['FIL']] = row['FMT']
        
        return fname_to_fieldnames[fname]
    
    def read_rrf(self, fname):
        """
        Read an rrf file returning a generator of dictionaries. Each
        dictionary represents a row with fieldnames as keys and row values
        as values.
        """
        fieldnames = self.fieldnames(fname)
        path = os.path.join(self.meta_dir, fname)
        with open(path) as f:
            reader = csv.DictReader(f, fieldnames=fieldnames, delimiter='|')
            for row in reader:
                yield row
        
class BaseMetathesaurus(Metathesaurus):
    
    def __init__(self, **kwargs):
        """This class is designed to provide access to the concept id, name,
        and symantic types for any concept. The profile all does not filter
        based on symantic type or source vocabularies but also does not retain
        source vocabulary annotations. A tab delimited text file concepts.txt
        stores id, name, and type information for all concepts. If concepts.txt
        exists, concepts are populated by parsing this file. If concepts.txt
        does not exist, the concept information is parsed from the META rrf
        files and concepts.txt is written to the /META directory.
        """
        super(BaseMetathesaurus, self).__init__(**kwargs)
        self.set_profile('all')
        self.concepts_path = os.path.join(self.meta_dir, 'concepts.txt')
        self.shelve_path = os.path.join(self.meta_dir, 'concepts.shelve')
        self.id_to_concept_shelve = shelve.open(self.shelve_path)
        if not os.path.exists(self.concepts_path):
            print 'Creating and writing', self.concepts_path
            self.read()
            self.write_concepts_file()
            self.update_shelve()
        
    def write_concepts_file(self):
        """Write a tab delimited text file with headers where each row
        represents a concept. Fields in order are concept_id, name, and
        symantic_types. symantic_types fields are saved as the str() conversion
        or a python list.
        """
        with open(self.concepts_path, 'w') as csv_file:
            writer = csv.writer(csv_file, delimiter='\t')
            writer.writerow(['concept_id', 'name', 'symantic_types'])
            for concept in self.concepts:
                writer.writerow([concept.concept_id, concept.name,
                                 list(concept.symantic_types)])
    
    def read_concepts_file(self):
        """Read tab delimited concept.txt file. Requies proper header with 
        concept_id, name, and symantic_types as fields. The symantic_type field
        assumes a valid str representation of a python list, set, or tuple of
        strings.
        """
        with open(self.concepts_path) as csv_file:
            concepts = set()
            reader = csv.DictReader(csv_file, delimiter='\t')
            for row in reader:
                concept = Concept(row['concept_id'])
                concept.name = row['name']
                symantic_types = ast.literal_eval(row['symantic_types'])
                concept.symantic_types = set(symantic_types)
                concepts.add(concept)
        self.concepts = concepts
        self.refresh_id_to_concept()

    def update_shelve(self):
        for concept_id, concept in self.id_to_concept.iteritems():
            self.id_to_concept_shelve[concept_id] = concept
        self.id_to_concept_shelve.sync()
        #self.id_to_concept_shelve.close()
        #self.id_to_concept_shelve = shelve.open(self.shelve_path)

    def get(self, concept_id, default=None):
        """Returns the concept with the specified concept_id. Returns None
        if no concept exists with that id."""
        return self.id_to_concept_shelve.get(concept_id)

    def __getitem__(self, concept_id):
        return self.id_to_concept_shelve[concept_id]

    def close(self):
        self.id_to_concept_shelve.close()


if __name__ == '__main__':
    #metathesaurus = omictools.data.metathesaurus
    #metathesaurus.set_profile('omicnet')
    #metathesaurus.read()
    #print metathesaurus
    base_metathesaurus = omictools.data.base_metathesaurus
    print base_metathesaurus.get('C0011849')
    print base_metathesaurus.get('C00118491')
    base_metathesaurus.close()


    