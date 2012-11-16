import networks.ontologies
from operator import attrgetter


class NCINode(networks.ontologies.AcyclicDirectedGraphNode):
    
    def __init__(self, code, concept_name, parents, synonyms, definition):
        """parents and synonyms are pipe delimited and should come directly
        from the flat nci text thesaurus.
        """
        super(NCINode, self).__init__()
        self.code = code
        self.concept_name = concept_name
        self.parent_concept_names = parents.split('|')
        self.synonyms = synonyms.split('|')
        self.definition = definition
    
    def __hash__(self):
        return hash(self.code)
    
    def __eq__(self, other):
        return self.code == other.code
        
    def __repr__(self):
        return ('Code: ' + self.code  + '\n' +
                'Name: ' + self.concept_name + '\n' +
                'Synonyms: ' + str(self.synonyms) + '\n' +
                'Defintion: ' + str(self.definition) + '\n' +
                'Parents: ' + str(self.parent_concept_names)  + '\n' +
                'Children: ' + str(map(attrgetter('concept_name'), self.children)))


class NCIOntology(networks.ontologies.AcyclicDirectedGraph):
    
    def __init__(self, path):
        """A DiseaseOntology object represents the ontology encoded by the
        NCI thesaurus. path is the path for the flattened text file.
        """
        super(NCIOntology, self).__init__()
        self.concept_name_dict = dict()
        self.code_name_dict = dict()
        self.NodeClass = NCINode
        self.path = path
    
    def read_nci_thesaurus(self):
        """Read a flattened nci thesaurus text file. File description:
        http://evs.nci.nih.gov/ftp1/NCI_Thesaurus/ReadMe.txt
        Download page: http://ncicb.nci.nih.gov/download/downloadevs.jsp
        """
        with open(self.path) as f:
            for line in f:
                self.nodes.add(self.NodeClass(*line.rstrip('\r\n').split('\t')))
        self.compute_dicts()
        self.construct_ontology_relations()
        print 'reading nci thesaurus complete'

    def compute_dicts(self):
        for node in self.nodes:
            self.concept_name_dict[node.concept_name] = node
            self.code_name_dict[node.code] = node

    def construct_ontology_relations(self):
        """Populate node.parents and node.children for each node and
        self.root_nodes.
        """
        for node in self.nodes:
            for parent_concept_name in node.parent_concept_names:
                if parent_concept_name == 'root_node':
                    self.root_nodes.add(node)
                else:
                    parent = self.concept_name_dict[parent_concept_name]
                    node.parents.add(parent)
                    parent.children.add(node)

