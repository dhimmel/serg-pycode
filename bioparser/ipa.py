import os
import csv

import networkx

import utilities.omictools
import networks.ontologies

import data

class Function(object):
    
    effect_kinds = {'affects', 'increases', 'decreases'}
    
    def __init__(self, name):
        """
        An ingenuity Function
        """
        self.name = name
        self.molecules = dict.fromkeys(self.effect_kinds, list())
        
    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, other):
        return self.name == other.name
    
    #def __str__(self):
    #    return self.name

    def __repr__(self):
        return str(self.__dict__)

class Molecule(object):

    def __init__(self, symbol, kind, synonyms, **kwargs):
        """
        set of all kinds:
        set(['chemical - kinase inhibitor', 'biologic drug', 'chemical drug', 'growth factor', 'translation regulator', 'kinase', 'chemical - other', 'enzyme', 'other', 'phosphatase', 'transcription regulator', 'ion channel', 'cytokine', 'chemical - endogenous mammalian', 'mature microRNA', 'chemical - protease inhibitor', 'transporter', 'peptidase', 'microRNA', 'chemical - endogenous non-mammalian', 'chemical reagent', 'transmembrane receptor', 'chemical toxicant', 'ligand-dependent nuclear receptor', 'G-protein coupled receptor'])
        """
        self.symbol = symbol
        self.kind = kind
        self.synonyms = synonyms
    
    def __hash__(self):
        return hash(self.symbol)
    
    def __eq__(self, other):
        return self.symbol == other.symbol
    
    def __str__(self):
        return self.symbol

    def __repr__(self):
        return self.symbol

class Drug(Molecule):
    
    def __init__(self, **kwargs):
        super(Drug, self).__init__(**kwargs)
        if 'targets' in kwargs:
            self.targets = kwargs['targets']

class Gene(Molecule):
    
    def __init__(self, entrez_gene_name, location, drugs, entrez_id_human, **kwargs):
        super(Gene, self).__init__(**kwargs)
        self.entrez_gene_name = entrez_gene_name
        self.location = location
        self.drugs = drugs
        self.entrez_id_human = entrez_id_human

class Other(Molecule):
    
    def __init__(self, **kwargs):
        super(Other, self).__init__(**kwargs)


class IPA(object):
    
    def __init__(self, ipa_dir=None):
        """
        Ingenuity IPA parser
        """
        if not ipa_dir:
            ipa_dir = data.current_path('ipa')
        self.ipa_dir = ipa_dir
        
        self.networkx = None
        
        self.drugs = set()
        self.molecules = set()
        self.genes = set()
        
    def build(self):
        """
        """
        print 'Building Ingenuity IPA'
        self.drugs = set(self.read_drugs('drugs.txt'))   
        self.molecules = set(self.read_molecules('molecules.txt'))
                
        self.functions = set()
        self.categories = set()
        function_generator = self.read_functions('functions.txt')
        for function in function_generator:
            self.categories.add(function.category)
            self.functions.add(function)
        self.name_to_function = {function.name: function for function in self.functions}
        print len(self.name_to_function), 'functions'
        self.build_ontology()
        
        # Incorporate effect on function information
        effect_generator = self.read_effects('effects.txt')
        for effect in effect_generator:
            name = effect['name']
            molecules = effect['molecules']
            function = self.name_to_function[name]
            function.molecules[effect['effect']] = molecules
        
        
        """
        # Parse drugs for each gene
        # NOT POSSIBLE UNTIL WE SURPASS THE 5000 result search limit for
        # genes and chemicals which affects the "chemical - endogenous mammal"
        # category.
        valid_symbols = {drug.symbol for drug in self.drugs}
        print len(valid_symbols), 'drugs'
        for gene in self.genes:
            drugs = gene.drugs
            drugs = self.parse_molecules(drugs, ', ', valid_symbols)
            gene.drugs = drugs
        """
        self.molecules = self.drugs | self.molecules
        self.symbol_to_molecule = {molecule.symbol: molecule for molecule in self.molecules}
        self.genes = {molecule for molecule in self.molecules if isinstance(molecule, Gene)}        
        self.others = {molecule for molecule in self.molecules if isinstance(molecule, Other)}        

    def build_ontology(self):
        """Build IPA Ontology"""
        self.ontology = IngenuityOntology()
        self.ontology.name_to_function = dict() # function name to node. does not include category, class, or root nodes
        function_category_supernode = IngenuityOntologyNode('function_category', 'root')
        function_class_supernode = IngenuityOntologyNode('function_class', 'root')
        for root_node in function_category_supernode, function_class_supernode:
            self.ontology.root_nodes.add(root_node)
            self.ontology.nodes.add(root_node)

        node_to_node = dict()
        for function in self.functions:
                        
            node_attributes = ([function.category, 'function_category'],
                               [function.function_class, 'function_class'],
                               [function.name, 'function'])
            
            ipa_kind_to_node = dict()
            for name, ipa_kind in node_attributes:
                new_node = IngenuityOntologyNode(name, ipa_kind)            
                node = node_to_node.get(new_node)
                if not node:
                    node = new_node
                    self.ontology.nodes.add(node)
                    node_to_node[node] = node
                ipa_kind_to_node[ipa_kind] = node
            
            function_category_node = ipa_kind_to_node['function_category']
            function_class_node = ipa_kind_to_node['function_class']
            function_node = ipa_kind_to_node['function']
            
            self.ontology.create_relationships([function_category_supernode], [function_category_node])
            self.ontology.create_relationships([function_class_supernode], [function_class_node])
            
            parents = [function_category_node, function_class_node]
            self.ontology.create_relationships(parents, [function_node])
            self.ontology.name_to_function[function_node] = function_node.name
        
        # Remove function_classes with a single child with the same name
        for node in list(self.ontology.nodes):
            if node.ipa_kind == 'function_category':
                children = list(node.children)
            if node.ipa_kind == 'function_class':
                children = list(node.children)
                if len(children) == 1 and node.name == children[0].name:
                    self.ontology.remove_node(node)
    
        
    def parse_molecules(self, unsplit, splitter, valid_symbols, total_number=None):
        """ """
        symbols = list()
        while unsplit:
            symbol, unsplit = self.split_one(unsplit, splitter)
            #print '-------------------------------------'
            while unsplit and symbol not in valid_symbols:
                #print symbol
                addition, unsplit = self.split_one(unsplit, splitter)
                symbol += splitter + addition
                
            if symbol not in valid_symbols:
                print symbol
                raise ValueError
            symbols.append(symbol)
        assert total_number is None or len(symbols) == total_number
        return symbols
    
    @staticmethod
    def split_one(unsplit, splitter):
        """Perform a single comma-space left string split. Returns the head 
        before the first split and the tail as the remainder of the string."""
        split = unsplit.split(', ', 1)
        if len(split) == 1:
            head = split[0]
            tail = ''
        else:
            head, tail = split
        return head, tail

    @staticmethod    
    def parse_function_annotation(function_annotation):
        try:
            name, synonyms = function_annotation.split(' [')
            synonyms = synonyms.rstrip(']')
            synonyms = synonyms.split(',')
            if '...' in synonyms:
                synonyms.remove('...')
        except ValueError:
            name = function_annotation
            synonyms = list()
        return name, synonyms
        
    def read_molecules(self, file_name='molecules.txt'):
        """
        With all results from functions query selected, click "Annotations". 
        A new window appears titled "Molecule Annotations". Click on 
        "Customize Table" and select all checkboxes to include "Synonym(s)".
        Click checkbox to select all annotations. Then click the export icon
        and save file in tab delimited format as molecules.txt. Close window.
        """
        assert self.drugs
        
        path = os.path.join(self.ipa_dir, file_name)
        fieldnames = ['symbol', 'synonyms', 'entrez_gene_name', 'location',
                      'kind', 'biomarker_applications', 'drugs',
                      'entrez_id_human', 'entrez_id_mouse', 'entrez_id_rat']
        with IPAExportReader(path, fieldnames) as dict_reader:
            for row in dict_reader:
                for key, value in row.items():
                    if value in {'', ' ', '--'}:
                        row[key] = None
                molecule = Other(**row)
                if molecule in self.drugs:
                    molecule = Drug(**row)
                if row['entrez_id_human'] and row['kind'] != 'microRNA':
                    molecule = Gene(**row)
                yield molecule
    
    def read_drugs(self, file_name='drugs.txt'):
        """
        Under "Genes and Chemicals" tab, click "Advanced Search".
        Under "Molecule Type(s) check "biologic drug" and "chemical drug".
        Export limit is 5000 items so future exports may need to do separate
        exports for biologic and chemical drugs. Click export and save file in
        tab delimited format as drugs.txt.
        
        To get a usuable symbol from ipa gene symbols run:
        symbol = gene.symbol.split('/')[0]
        symbol = symbol.split(' (')[0]
        """
        path = os.path.join(self.ipa_dir, file_name)
        fieldnames = ['number', 'symbol', 'matched_term', 'synonyms', 'entrez',
                      'location', 'kind', 'biomarker_applications', 'drugs',
                      'targets']
        plural_fields = {'targets'}
        with IPAExportReader(path, fieldnames) as dict_reader:
            for row in dict_reader:
                for key, value in row.items():
                    if value == '' or value == ' ':
                        value = None
                        row[key] = value
                    if key in plural_fields:
                        row[key] = list() if value is None else value.split(',')
                drug = Drug(**row)
                yield drug
      
    def read_functions(self, file_name='functions.txt'):
        """
        Under the "Functions and Diseases" tab query "disease". Ensure the
        display mode is set to show categories (default) as opposed to a simple
        list of functions. If a "Show categories" button is visible, it must
        be clicked. The query matches 9805 functions and diseases. 
        Check "Matching Functions & Diseases" to select all. Click export
        and save file in tab delimited format as drugs.txt.
        """
        path = os.path.join(self.ipa_dir, file_name)
        fieldnames = ['category','function_class', 'function_annotation', 'molecules',
                      'number_of_molecules']
        with IPAExportReader(path, fieldnames) as dict_reader:
            for row in dict_reader:
                
                # Function annotation and creation of name and synonyms
                name, synonyms = self.parse_function_annotation(row['function_annotation'])
                row['name'] = name
                row['synonyms'] = synonyms
                del row['function_annotation']
                
                # Molecules
                number_of_molecules = int(row['number_of_molecules'])
                row['number_of_molecules'] = number_of_molecules
                function = Function(row['name'])
                del row['molecules']
                function.__dict__.update(row)
                yield function

    def read_effects(self, file_name='effects.txt'):
        """
        With all results from a function query selected click
        "Effect on Function". A dialog box which says "Opening" appears.
        Search takes several minutes to complete. When complete a new window
        opens titled "Effect on function". Check box to select all processes
        and export as a text file named effects.txt.
        """
        
        assert self.molecules
        
        path = os.path.join(self.ipa_dir, file_name)
        fieldnames = ['function_annotation','molecules', 'number_of_molecules']
        valid_symbols = {molecule.symbol for molecule in self.molecules}
        with IPAExportReader(path, fieldnames) as dict_reader:
            for row in dict_reader:
                function_annotation = row['function_annotation']
                effect, function_annotation = function_annotation.split(' ', 1)
                assert effect in Function.effect_kinds
                row['effect'] = effect
                function_annotation, occurances = function_annotation.rsplit(' ', 1)
                #occurances = occurances.strip('()')
                #current_num_of_molecules, total_num_of_molecules = occurances.split('/')
                #current_num_of_molecules = int(current_num_of_molecules)
                #total_num_of_molecules = int(total_num_of_molecules)

                name, synonyms = self.parse_function_annotation(function_annotation)
                del row['function_annotation']
                
                number_of_molecules = int(row['number_of_molecules'])
                row['number_of_molecules'] = number_of_molecules
                unsplit = row['molecules']
                molecules = self.parse_molecules(unsplit, ', ', valid_symbols, number_of_molecules)                
                row['molecules'] = molecules
                row['name'] = name
                row['function'] = self.name_to_function[name]
                yield row
        
        
        
class IngenuityOntology(networks.ontologies.AcyclicDirectedGraph):
    
    def __init__(self):
        """
        """
        super(IngenuityOntology, self).__init__()
        

class IngenuityOntologyNode(networks.ontologies.AcyclicDirectedGraphNode):
    
    def __init__(self, name, ipa_kind):
        """
        """
        assert ipa_kind in ['root', 'function_category', 'function_class', 'function']
        super(IngenuityOntologyNode, self).__init__()
        self.name = name
        self.ipa_kind = ipa_kind
    
    def __hash__(self):
        return hash((self.name, self.ipa_kind))
    
    def __eq__(self, other):
        return self.name == other.name and self.ipa_kind == other.ipa_kind
    
    def __repr__(self):
        return self.name

class IPAExportReader(object):
    
    none_strings = {'', ' ', '--'}
    
    def __init__(self, path, fieldnames):
        self.path = path
        self.fieldnames = fieldnames
    
    def __enter__(self):
        self.f = open(self.path)
        for i in xrange(3):
            self.f.next() # skip line
        reader = csv.DictReader(self.f, fieldnames=self.fieldnames, delimiter='\t')
        for row in reader:
            if None in row:
                del row[None]
            for key, value in row.items():
                row[key] = self.to_ascii(value)
            yield row
    
    def __exit__(self, *args, **kwargs):
        self.f.close()

    @staticmethod
    def to_ascii(s):
        return unicode(s, encoding='utf-8').encode('ascii', 'replace')

if __name__ == '__main__':
    ipa = IPA()
    print ipa.ipa_dir
    ipa.build()
    
    
    g = ipa.build_networkx()
    
    #onto_structure_path = os.path.join(ipa.ipa_dir, 'ipa-ontology-structure.txt')
    #ipa.ontology.write_structure(onto_structure_path)    
    
