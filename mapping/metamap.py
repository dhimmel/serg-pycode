import os
import csv
import shutil
import subprocess
import time

import utilities.omictools

class Mapping(object):
    
    # fieldnames for row representation of object
    fieldnames = ['score', 'input', 'concept_name', 'concept_id']
    
    def __init__(self, input, concept_name=None, concept_id=None, score=None):
        """
        """
        self.input = input
        self.concept_id = concept_id
        self.concept_name = concept_name
        self.score = score
        if self.score is not None:
            self.score = int(self.score)
        self.tuple_ = (self.input, self.concept_id)
    
    def __hash__(self):
        return hash(self.tuple_)

    def __eq__(self, other):
        return self.tuple_ == other.tuple_

    def __ne__(self, other):
        return self.tuple_ != other.tuple_

    def __lt__(self, other):
        return self.score < other.score
    
    def __str__(self):
        return self.input + ' <---> ' + str(self.concept_name)
    
    def __repr__(self):
        return self.__str__()
    
    def as_row(self):
        return [getattr(self, field) for field in Mapping.fieldnames]

class MetaMap(object):
    
    def __init__(self, map_dir, previous_map_dir=None, fname_suffix='',
                 metamap_program_dir='/opt/public_mm/',
                 restrict_to_sources=['MSH'],
                 restrict_to_types=['dsyn','neop','anab']):
        """
        pre-edit -- automatically matched mappings either by metamapper 
            or records of previous manual curation
        redacted -- pre-edit file after manual removals and additions

        removed -- cumulative mappings that have been manually removed
        added -- cumulative mappings that have been manually added
        
        Semantic Type Abbreviations
        http://mmtx.nlm.nih.gov/semanticTypes.shtml
        """
        
        self.metamap_program_dir = os.path.expanduser(metamap_program_dir)
        self.metamap_bin = os.path.join(metamap_program_dir, 'bin/metamap12')
        
        if not os.path.isdir(map_dir):
            os.mkdir(map_dir)
        
        self.map_dir = map_dir
        self.previous_map_dir = previous_map_dir
        
        self.restrict_to_sources = restrict_to_sources
        self.restrict_to_types = restrict_to_types
        
        self.paths = dict()
        self.maps = dict()
        
        categories = ['removed', 'added', 'pre-edit', 'redacted']
        version_dict = {'current': map_dir, 'previous': previous_map_dir}
        
        if fname_suffix:
            fname_suffix = '-' + fname_suffix
        
        for version in version_dict.keys():
            
            for d in self.paths, self.maps:
                d[version] = dict().fromkeys(categories)
            
            for category in categories:
                if version != 'previous' or previous_map_dir:
                    head = version_dict[version]
                    tail = 'metamappings-' + category + fname_suffix + '.txt'
                    path = os.path.join(head, tail)
                    self.paths[version][category] = path
        
        if os.path.exists(self.paths['current']['redacted']):
            mappings = self.read_mapping_file(self.paths['current']['redacted'])
            self.maps['current']['redacted'] = mappings

        if previous_map_dir:
            for category in categories:
                path = self.paths['previous'][category]
                if os.path.isfile(path):
                    mappings = self.read_mapping_file(path)
                    self.maps['previous'][category] = mappings
                else:
                    self.write_mapping_file(set(), path)
        
                
    def start_servers(self):
        """Start the requried servers for metamap."""
        for server_name in ['skrmedpostctl', 'wsdserverctl']:
            server_path = os.path.join(self.metamap_program_dir, 'bin', server_name)
            subprocess.call([server_path, 'start'], stdout=subprocess.PIPE)
            time.sleep(1)

    def stop_servers(self):
        """Stop the requried servers for metamap."""
        for server_name in ['skrmedpostctl', 'wsdserverctl']:
            server_path = os.path.join(self.metamap_program_dir, 'bin', server_name)
            subprocess.call([server_path, 'stop'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            time.sleep(1)

    def metamap_term(self, term):
        """Use metamap to map a single term to a UMLS comcept"""
        restrict_to_sources = ','.join(self.restrict_to_sources)
        restrict_to_sts = ','.join(self.restrict_to_types)
        quoted_term = '"' + term + '"'
        arg_list = ['echo', quoted_term, '|', self.metamap_bin, '--restrict_to_sts ' + restrict_to_sts,
                    #'--restrict_to_sources ICD10CM,ICD9CM,MTHICD9,MSH,NCI,OMIM',
                    '--restrict_to_sources ' + restrict_to_sources,
                    '--silent', '--show_cuis', '--term_processing']
        cmd = ' '.join(arg_list)
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.PIPE)
        result_dict = self.parse_metamap_output(output)
        return Mapping(term, **result_dict)
    
    def parse_metamap_output(self, output):
        """Parse output from a command line metamaop query."""
        pull_line = False
        for line in output.split('\n'):
            if pull_line:
                meta_mapping = line.strip()
                score, concept = meta_mapping.split('  ')
                score = int(score)
                concept_id, concept_name = concept.split(':')
                result_dict = {'concept_id': concept_id,
                               'concept_name': concept_name, 'score':score}
                return result_dict
            elif line.startswith('Meta Mapping'):
                pull_line = True
        return {'concept_id': None, 'concept_name': None, 'score': None}

    def map_terms(self, terms, score_threshold=0, concepts_ids_to_remove=None):
        """Map an iterable of terms (strings) to UMLS concepts"""
        self.start_servers()
        
        terms = set(terms)        
        mappings = set(self.metamap_term(term) for term in terms)
        if concepts_ids_to_remove is None:
            concepts_ids_to_remove = {'C0012634', 'C0039082', 'C0162429'}
        
        for mapping in list(mappings):
            if (mapping.score < score_threshold or
                mapping.concept_id in concepts_ids_to_remove):
                mappings.remove(mapping)
                mappings.add(Mapping(mapping.input))
            
        self.stop_servers()
        
        if self.previous_map_dir:
            
            # Filter out previously removed mappings
            previously_removed = self.maps['previous']['removed']
            mappings -= previously_removed
            
            # Add previous manual mappings
            # Previous manual mapping overrides current automatic mapping
            input_to_mapping = {mapping.input: mapping for mapping in mappings}
            for manual_map in self.maps['previous']['added']:
                input_to_mapping[manual_map.input] = manual_map
            
            # Add previous automatic mappings
            # Previous automatic mappings do not override current automatic mappings
            for auto_map in self.maps['previous']['redacted']:
                input = auto_map.input
                if input not in input_to_mapping:
                    input_to_mapping[input] = auto_map
            
            mappings = set(input_to_mapping.values())

        self.maps['current']['pre-edit'] = mappings
        self.write_mapping_file(mappings, self.paths['current']['redacted'])
        mappings = self.excluding_unmapped(mappings)
        self.write_mapping_file(mappings, self.paths['current']['pre-edit'])
        print 'Perform manual redaction and then run read_mappings()'
        
    def read_mappings(self):
        """Returns an input to concept_id dict. Filtered terms are excluded"""
        
        for category in ['pre-edit', 'redacted']:
            mappings = self.read_mapping_file(self.paths['current'][category])
            self.maps['current'][category] = mappings
            
        # Create removed and added mappings files
        # mappings that appear in pre-edit but not redacted are removed
        removed = self.maps['current']['pre-edit'] - self.maps['current']['redacted']
        # add mappings that appear in redacted and not pre-edit
        added = self.maps['current']['redacted'] - self.maps['current']['pre-edit']
        
        if self.previous_map_dir:
            previous_added = self.maps['previous']['added']
            previous_removed = self.maps['previous']['removed'] 
            previous_added -= removed
            previous_removed -= added
            # mappings were previousely added get added to the added list
            added |= previous_added
            # mappings were previousely removed get added to the removed list
            removed &= previous_removed        

        self.maps['current']['removed'] = removed
        self.write_mapping_file(removed, self.paths['current']['removed'])

        self.maps['current']['added'] = added
        self.write_mapping_file(added, self.paths['current']['added'])

    @staticmethod
    def read_mapping_file(path):
        mappings = set()
        f = open(path)
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            for key, value in row.items():
                if value is '':
                    row[key] = None
            mapping = Mapping(**row)
            mappings.add(mapping)
        f.close()
        mappings = MetaMap.excluding_unmapped(mappings)
        return mappings
    
    @staticmethod
    def write_mapping_file(mappings, path):
        mappings = list(mappings)
        #mappings = sorted(mappings, key = lambda x: x.score, reverse=True)
        mappings.sort(reverse=True)
        with open(path, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(Mapping.fieldnames)
            for mapping in mappings:
                writer.writerow(mapping.as_row())
        
    def pause_before_reading_mappings(self):
        while raw_input('Type <proceed> when editing of mapping-redacted.txt is complete.\n') != 'proceed':
            pass
        self.read_mappings()
        
    def get_mapping_dict(self, key='input', value='concept_id'):
        return {getattr(mapping, key): getattr(mapping, value)
                for mapping in self.maps['current']['redacted'] if mapping.concept_id}

    @staticmethod
    def excluding_unmapped(mappings):
        return set(mapping for mapping in mappings if mapping.concept_id)


if __name__ == '__main__':
    mm = MetaMap('/home/dhimmels/Desktop/mm_test/')
    terms = ['multiple sclerosis', 'ear problems', 'diarrhea']
    mm.map_terms(terms)
    mm.pause_before_reading_mappings()
    #mm.read_mappings()
    print mm.get_mapping_dict()