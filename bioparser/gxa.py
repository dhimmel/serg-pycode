#http://rest.bioontology.org/bioportal/virtual/ontology/1136?conceptid=efo:EFO_0000408&apikey=041eb758-d944-4493-bfe2-3dfd6830d2cd
from ontologies import AcyclicDirectedGraph, AcyclicDirectedGraphNode
import urllib2
import json
import pickle
import os
import math
import omicnet

class ExpirimentalFactorOntology(AcyclicDirectedGraph):

    def __init__(self):
        super(ExpirimentalFactorOntology, self).__init__()
        self.NodeClass = ExpirimentalFactorOntologyNode
        self.id_dict = dict()
        self.name_dict = dict()
    
    def read_efo_obo(self, path):
        with open(path) as f:
            line = f.readline()
            while line:
                if line.startswith('[Term]'):
                    id = f.readline().split(' ! ')[0][4:]
                    id = id.replace(':', '_')
                    name = f.readline().strip()[6:]
                    parents = []
                    line = f.readline().strip()
                    while line.startswith('is_a:'):
                        parents.append(line[6:].replace(':', '_'))
                        line = f.readline().strip()
                    self.nodes.add(self.NodeClass(id, name, parents))
                line = f.readline()
        for node in self.nodes:
            self.id_dict[node.id] = node
            self.name_dict[node.name] = node
        self.construct_ontology_relations()
        print 'reading efo obo complete'
    
    def construct_ontology_relations(self):
        """Populate node.parents and node.children for each node and
        self.root_nodes.
        """
        for node in self.nodes:
            if not node.parent_ids:
                self.root_nodes.add(node)
            else:
                for parent_id in node.parent_ids:
                    parent = self.id_dict[parent_id]
                    node.parents.add(parent)
                    parent.children.add(node)

class ExpirimentalFactorOntologyNode(AcyclicDirectedGraphNode):
    
    def __init__(self, id, name, parents):
        super(ExpirimentalFactorOntologyNode, self).__init__()
        self.id = id
        self.name = name
        self.parent_ids = parents
    
    def __hash__(self):
        return hash(self.id)
    
    def __repr__(self):
        return '%s\t%s' % (self.id, self.name)


class GXA:
    
    def __init__(self, 
            gxa_query_dir='/Users/dhimmels/Documents/serg/omicnet/input-datasets/gxa/queries/120724/', 
            efo_obo_path='/Users/dhimmels/Documents/serg/omicnet/input-datasets/gxa/efo/trunk/src/efoinobo/efo.obo'):
        self.gxa_query_dir = gxa_query_dir
        self.efo_obo_path = efo_obo_path
    
    def read_efo(self):
        self.efo = ExpirimentalFactorOntology()
        self.efo.read_efo_obo(self.efo_obo_path)

    def query_gxa(self, query_efo_id, rows=200):
        expr_dict = dict()
        starting_from = 0
        while True:
            url = ('http://www.ebi.ac.uk/gxa/api/v1?updownInEfo=' + query_efo_id + 
                   '&species=Homo+Sapiens&format=json&start=' + str(starting_from) + '&rows=' + str(rows))
            #print url
            try:
                usock = urllib2.urlopen(url)
                data = json.load(usock)
                usock.close()
                total_results = data['totalResults']
                starting_from = data['startingFrom']
                print "Parsing results %s to %s of %s" % (starting_from, starting_from + rows, total_results)
                for result in data['results']:
                    gene_name = result['gene']['name']
                    for expression in result['expressions']:
                        efo_term = expression['efoTerm']
                        efo_id = expression['efoId']
                        preserve_keys = set(['upPvalue', 'downPvalue', 'upExperiments', 'downExperiments', 'nonDEExperiments'])
                        gene_dict = {key: value for key, value in expression.items() if key in preserve_keys}
                        expr_dict.setdefault(efo_term, dict())[gene_name] = gene_dict
                starting_from += rows
                if starting_from >= total_results:
                    break
            except urllib2.HTTPError:
                print 'urllib2.HTTPError'
        return expr_dict

    def full_query(self):
        """Perform GXA Queries"""
        gxa_nodes = self.efo.descendants_of_depth(efo.id_dict['EFO_0000408'], 2)
        for node in gxa_nodes:
            print 'Starting gxa query with parent node:', node
            expr_dict = self.query_gxa(node.id)
            if not expr_dict:
                continue
            path = os.path.join(self.gxa_query_dir, node.id + '.json')
            with open(path, 'w') as json_file:
                json.dump(expr_dict, json_file)

    def read_gxa_json_queries(self, condition_subset=None):
        print 'Reading gxa queries:', 
        files = os.listdir(self.gxa_query_dir)
        files = filter(lambda s: s.endswith('.json'), files)
        gxa_dict = dict()
        for file in files:
            path = os.path.join(self.gxa_query_dir, file)
            with open(path) as f:
                expr_dict = json.load(f)
                for key, value in expr_dict.iteritems():
                    if (condition_subset is None
                        or key in condition_subset):
                        gxa_dict[key] = value
        
        for disease, gene_dict in gxa_dict.iteritems():
            for expr_dict in gene_dict.itervalues():
                total_experiments = float(expr_dict['downExperiments'] + 
                                          expr_dict['upExperiments'] + 
                                          expr_dict['nonDEExperiments'])
                expr_dict['totalExperiments'] = int(total_experiments)
                expr_dict['percentUp'] = expr_dict['upExperiments'] / total_experiments
                expr_dict['percentNonDE'] = expr_dict['nonDEExperiments'] / total_experiments
                expr_dict['percentDown'] = expr_dict['downExperiments'] / total_experiments
        self.gxa_dict = gxa_dict
        print 'complete.'
        return gxa_dict

    def filter_gxa_dict(self, p_value_cutoff=0.001):
        """Filters genes that do not meet an up or down p-value cutoff."""
        for disease, gene_dict in self.gxa_dict.iteritems():
            for gene, expr_dict in gene_dict.items():
                if expr_dict['upPvalue'] > p_value_cutoff and expr_dict['downPvalue'] > p_value_cutoff:
                    del gene_dict[gene]
        return self.gxa_dict
    
    def get_genes_per_condition(self):
        return {disease: len(gene_dict) for disease, gene_dict in self.gxa_dict.iteritems()} 

    def get_disease_to_genes_dict(self, entrez=False):
        disease_to_genes = {disease: gene_dict.keys() for disease, gene_dict in self.gxa_dict.iteritems()}
        if entrez:
            neo_db_dir = '/Users/dhimmels/Documents/serg/omicnet/neo4j-omicnet-db'
            g = omicnet.OmicNet(neo_db_dir)
            all_genes = set()
            for genes in disease_to_genes.itervalues():
                all_genes.update(genes)
            all_genes = set(filter(lambda x: ':' not in x, all_genes))
            symbol_to_entrez = dict()
            for symbol in all_genes:
                node = g.retrieve_gene(symbol)
                #if not node:
                #    print symbol
                if node and 'Entrez Gene ID' in node.keys():
                    entrez = node['Entrez Gene ID']
                    if entrez:
                        symbol_to_entrez[symbol] = entrez
            g.shutdown()
            disease_to_entrez_genes = dict()
            for disease, genes in disease_to_genes.iteritems():
                entrez_genes = [symbol_to_entrez[gene] for gene in genes if gene in symbol_to_entrez]
                disease_to_entrez_genes[disease] = set(entrez_genes)
            disease_to_genes = disease_to_entrez_genes                
        return disease_to_genes
        
if __name__ == '__main__':
    ashg_disease_dict = {'amyotrophic lateral sclerosis': 'ALS',
                         'breast carcinoma': 'BreastCancer',
                         'cataract': 'Cataract',
                         'Parkinson\'s disease': 'Coriell_Parkinsons',
                         "Alzheimer's disease": 'LOAD',
                         'lung carcinoma': 'LungCancer',
                         'melanoma': 'Melanoma',
                         'pancreatic carcinoma': 'PanScan',
                         'prostate carcinoma': 'ProstateCancer',
                         'psoriasis': 'Psoriasis',
                         'schizophrenia': 'Schizophrenia_AA',
                         'obesity': 'VisceralAdiposity'}
    
    gxa = GXA()
    gxa.read_gxa_json_queries(ashg_disease_dict)
    print 'Before filtering:', gxa.get_genes_per_condition()
    gxa.filter_gxa_dict()
    print 'After filtering:', gxa.get_genes_per_condition()
    disease_to_symbols = gxa.get_disease_to_genes_dict()
    disease_to_entrez = gxa.get_disease_to_genes_dict(entrez=True)
    print 'After entrez conversion:', {disease: len(genes) for disease, genes in disease_to_entrez.iteritems()}    
    
    ashg_dir = '/Users/dhimmels/Documents/ashg-gene-sets'
    for disease, entrez_ids in disease_to_entrez.items():
        ashg_disease = ashg_disease_dict[disease]
        path = os.path.join(ashg_dir, ashg_disease + '.txt')
        with open(path, 'w') as f:
            f.writelines(id + '\n' for id in entrez_ids)
        

"""
def meta_p(p_values):
    n = len(p_values)
    k = prod(p_values)
    meta_p = k * sum((-math.log(k)) ** i / math.factorial(i) for i in range(n))
    return meta_p

for disease, gene_dict in gxa_dict.iteritems():
    
    if disease not in ashg_diseases: continue
    num_consistent_05 = sum(max([expr_dict['percentUp'], expr_dict['percentDown']]) == 1.0
              for expr_dict in gene_dict.itervalues())
    num_consistent_01 = sum(max([expr_dict['percentUp'], expr_dict['percentDown']]) == 1.0 and
              (expr_dict['upPvalue'] < 0.01 or 
              expr_dict['downPvalue'] < 0.01)
              for expr_dict in gene_dict.itervalues())    
    num_consistent_001 = sum(max([expr_dict['percentUp'], expr_dict['percentDown']]) == 1.0 and
              (expr_dict['upPvalue'] < 0.001 or 
              expr_dict['downPvalue'] < 0.001)
              for expr_dict in gene_dict.itervalues())    
    num_consistent_multiple_exp = sum(max([expr_dict['percentUp'], expr_dict['percentDown']]) == 1.0 and
              (expr_dict['totalExperiments'] > 1)
              for expr_dict in gene_dict.itervalues())
    print '%30s %8s %8s %8s %8s %8s' % (disease, len(gene_dict), num_consistent_05, num_consistent_01, num_consistent_001, num_consistent_multiple_exp)

print '%30s %8s %8s %8s %8s %8s' % ('disease', 'all_0.05', 'same_0.05', 'same_0.01', 'same_0.001', 'same_>1exp')
"""   
    






