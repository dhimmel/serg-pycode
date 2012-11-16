import os
import csv

import utilities.omictools
import networks.ontologies

import data

class Function(object):
    
    effect_types = {'affects', 'increases', 'decreases'}
    
    def __init__(self, name):
        """
        An ingenuity Function
        """
        self.name = name
        self.molecules = dict.fromkeys(self.effect_types, list())
        
    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, other):
        return self.name == other.name
    
    #def __str__(self):
    #    return self.name

    def __repr__(self):
        return str(self.__dict__)

class Molecule(object):

    def __init__(self, symbol):
        """
        set of all types:
        set(['chemical - kinase inhibitor', 'biologic drug', 'chemical drug', 'growth factor', 'translation regulator', 'kinase', 'chemical - other', 'enzyme', 'other', 'phosphatase', 'transcription regulator', 'ion channel', 'cytokine', 'chemical - endogenous mammalian', 'mature microRNA', 'chemical - protease inhibitor', 'transporter', 'peptidase', 'microRNA', 'chemical - endogenous non-mammalian', 'chemical reagent', 'transmembrane receptor', 'chemical toxicant', 'ligand-dependent nuclear receptor', 'G-protein coupled receptor'])
        """
        self.symbol = symbol
    
    def __hash__(self):
        return hash(self.symbol)
    
    def __eq__(self, other):
        return self.symbol == other.symbol
    
    def __str__(self):
        return self.symbol

    def __repr__(self):
        return self.symbol

class Drug(Molecule):
    
    def __init__(self, symbol):
        super(Drug, self).__init__(symbol)

class IPAGene(Molecule):
    
    def __init__(self, symbol):
        super(Drug, self).__init__(symbol)

class Other(Molecule):
    
    def __init__(self, symbol):
        super(Other, self).__init__(symbol)


class IPA(object):
    
    def __init__(self, ipa_dir=None):
        """
        Ingenuity IPA parser
        """
        if not ipa_dir:
            ipa_dir = data.current_path('ipa')
        self.ipa_dir = ipa_dir
    
    def build(self):
        """
        """
        print 'Building Ingenuity IPA'
                
        molecule_generator = self.read_annotations('annotations-query_disease.txt')
        #for molecule in molecules:
        self.molecules = set(molecule_generator)
            
        
        
        self.symbol_to_molecule = {molecule.symbol: molecule for molecule in self.molecules}
        self.functions = set()
        function_generator = self.read_functions('associated_molecules-query_disease.txt')
        for function in function_generator:
            self.functions.add(function)
        self.name_to_function = {function.name: function for function in self.functions}
        print len(self.name_to_function), 'functions'
        self.build_ontology()

        # Incorporate effect on function information
        effect_generator = self.read_effect_on_function('effect_on_function-query_disease.txt')
        for effect in effect_generator:
            name = effect['name']
            molecules = effect['molecules']
            function = self.name_to_function[name]
            function.molecules[effect['effect']] = molecules
        
        self.drugs = set(self.read_drugs('drugs.txt'))
        self.genes = data.Data().hgnc.get_genes()


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
            if node.ipa_kind == 'function_class':
                children = list(node.children)
                if len(children) == 1 and node.name == children[0].name:
                    self.ontology.remove_node(node)

    def parse_molecules(self, unsplit, total_number):
        """ """
        molecules = list()
        while unsplit:
            symbol, unsplit = self.split_one(unsplit)
            while unsplit and symbol not in self.symbol_to_molecule:
                addition, unsplit = self.split_one(unsplit)
                symbol += ', ' + addition
            molecule = self.symbol_to_molecule[symbol]
            molecules.append(molecule)
        assert len(molecules) == total_number
        return molecules
    
    
    @staticmethod
    def split_one(unsplit, splitter=', '):
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
        
    
    def read_annotations(self, file_name):
        """
        """
        path = os.path.join(self.ipa_dir, file_name)
        fieldnames = ['symbol','synonyms', 'entrez_gene_name', 'location',
                      'type', 'biomarker_applications', 'drugs',
                      'entrez_id_human', 'entrez_id_mouse', 'entrez_id_rat']
        with IPAExportReader(path, fieldnames) as dict_reader:
            for row in dict_reader:
                #if row['entrez_gene_name'] == '--':
                #    row['entrez_gene_name'] == None
                for key, value in row.items():
                    if value in {'', ' ', '--'}:
                        row[key] = None
                molecule = Molecule(row['symbol'])
                molecule.__dict__.update(row)
                yield molecule
    
    def read_drugs(self, file_name):
        """
        Under "genes and chemicals" tab click advanced search. Under "molecule
        types check "biological drug" and "chemical drug". Then search and
        export.
        """
        path = os.path.join(self.ipa_dir, file_name)
        fieldnames = ['number', 'symbol', 'matched_term', 'synonyms', 'entrez',
                      'location', 'type', 'biomarker_applications', 'drugs',
                      'targets']
        delete_fields = {'number', 'matched_term', 'entrez', 'location',
                         'biomarker_applications', 'drugs'}
        plural_fields = {'targets'}
        with IPAExportReader(path, fieldnames) as dict_reader:
            for row in dict_reader:
                for fieldname in delete_fields:
                    del row[fieldname]
                for key, value in row.items():
                    if value == '':
                        row[key] = None
                    if key in plural_fields:
                        row[key] = list() if value is None else value.split(',')
                drug = Drug(row['symbol'])
                drug.__dict__.update(row)
                yield drug
        
        
    def read_functions(self, file_name):
        """
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
                """
                unsplit = row['molecules']
                molecules = self.parse_molecules(unsplit, number_of_molecules)
                row['molecules'] = molecules
                """              
                function = Function(row['name'])
                del row['molecules']
                function.__dict__.update(row)
                yield function

    def read_effect_on_function(self, file_name):
        """
        Query "disease" in functions and disease searchbar (top middle tab).
        Select all results and click "effect on function".
        Select all processes and export
        """
        path = os.path.join(self.ipa_dir, file_name)
        fieldnames = ['function_annotation','molecules', 'number_of_molecules']
        with IPAExportReader(path, fieldnames) as dict_reader:
            for row in dict_reader:
                function_annotation = row['function_annotation']
                effect, function_annotation = function_annotation.split(' ', 1)
                assert effect in Function.effect_types
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
                molecules = self.parse_molecules(unsplit, number_of_molecules)                
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
    
    def __init__(self, path, fieldnames):
        self.path = path
        self.fieldnames = fieldnames
    
    def __enter__(self):
        self.f = open(self.path)
        for i in xrange(3):
            self.f.next() # skip line
        reader = csv.DictReader(self.f, fieldnames=self.fieldnames, delimiter='\t')
        return reader
    
    def __exit__(self, *args, **kwargs):
        self.f.close()


if __name__ == '__main__':
    ipa = IPA()
    ipa.build()
    #onto_structure_path = os.path.join(ipa.ipa_dir, 'ipa-ontology-structure.txt')
    #ipa.ontology.write_structure(onto_structure_path)
    for i, molecule in enumerate(ipa.molecules):
        if molecule.entrez_id_human:
            gene = data.Data().hgnc.get_symbol_to_gene().get(molecule.symbol)
            #if not gene:
            #    print molecule.symbol, '\t', molecule.entrez_id_human
            #print molecule.symbol, '\t', molecule.type, '\t', molecule.entrez_id_human
            if gene:
                drugs = molecule.drugs
                if not drugs:
                    continue
                print molecule.symbol, drugs


    #print(ipa.name_to_node['phospholipidosis of lysosome'].children)
    #print list(ipa.functions)[1]    
    