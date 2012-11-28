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
        self.drugs = set(self.read_drugs())   
        self.molecules = set(self.read_annotations('annotations-query_disease.txt'))
                
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
            if node.ipa_kind == 'function_class':
                children = list(node.children)
                if len(children) == 1 and node.name == children[0].name:
                    self.ontology.remove_node(node)
    
    def build_networkx(self):
        """
        """
        if self.networkx:
            return self.networkx
        self.build()
        
        hgnc = data.Data().hgnc
        entrez_to_hgnc = hgnc.get_entrez_to_gene()
        symbol_to_hgnc = hgnc.get_symbol_to_gene()
                
        g = networkx.Graph(name='ipanet')
        ################################################################################
        ################################# Create Nodes #################################
        
        # Create drug nodes
        targets = set()
        for drug in self.drugs:
            g.add_node(drug.symbol, kind='drug')
            targets |= set(drug.targets)
        
        # Create disease nodes
        for disease in self.functions:
            g.add_node(disease.name, kind='disease')
        
        # Create gene nodes
        hugu_genes_added = set()
        for gene in self.genes:
            entrez_id = gene.entrez_id_human.split('|')[0]
            hgnc_gene = entrez_to_hgnc.get(entrez_id)
            if hgnc_gene:
                hugu_genes_added.add(hgnc_gene)
            hugu_symbol = hgnc_gene.symbol if hgnc_gene else None
            g.add_node(gene.symbol, hugu_symbol=hugu_symbol, kind='gene')
        for target in targets:
            if target not in g:
                hgnc_gene = symbol_to_hgnc.get(target)
                if hgnc_gene:
                    hugu_genes_added.add(hgnc_gene)
                hugu_symbol = hgnc_gene.symbol if hgnc_gene else None
                g.add_node(target, hugu_symbol=hugu_symbol, kind='gene')
        del targets
        
        missing_hugu_genes = set(hgnc.get_genes()) - hugu_genes_added
        for hgnc_gene in missing_hugu_genes:
            if hgnc_gene.symbol in g:
                raise Exception('pre-existing ipa symbol matching gene name')
            g.add_node(hgnc_gene.symbol, hugu_symbol=hgnc_gene.symbol, kind='gene')
        
        ################################################################################
        ################################# Create Edges #################################
        # Create drug-gene links from drug target annotations.
        for drug in self.drugs:
            for target in drug.targets:
                g.add_edge(drug.symbol, target, kind='target')
        
        # Create disease-gene and disease-drug links from ipa function annotations.
        for disease in self.functions:
            for effect, molecules in disease.molecules.items():
                for molecule in molecules:
                    if disease.name in g and molecule in g:
                        if g.node[molecule]['kind'] == 'drug':
                            if effect == 'decreases':
                                kind = 'indication'
                            else:
                                kind = 'disease_modifying_drug'
                        else:
                            kind = 'disease_gene'
                        g.add_edge(disease.name, molecule, effect=effect, kind=kind)
        
        self.networkx = g
        return self.networkx
    
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
        
    def read_annotations(self, file_name):
        """
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
        Under "genes and chemicals" tab click advanced search. Under "molecule
        kinds check "biological drug" and "chemical drug". Then search and
        export.
        
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

    g = ipa.build_networkx()
    print networkx.info(g)
    connected_nodes = (node for node, degree in g.degree_iter() if degree)
    g_connected = g.subgraph(connected_nodes)
    print networkx.info(g_connected)
    gml_path = '/home/dhimmels/Documents/serg/ipanet/ipanet.gml'
    networkx.write_gml(g_connected, gml_path)
    print 'IPA network written as GML'
    pkl_path = '/home/dhimmels/Documents/serg/ipanet/ipanet.pkl'
    networkx.write_gpickle(g_connected, pkl_path)
    print 'IPA network written as pickle'
    
    """
    ipa.build()
    onto_structure_path = os.path.join(ipa.ipa_dir, 'ipa-ontology-structure.txt')
    ipa.ontology.write_structure(onto_structure_path)

    """
    
    
