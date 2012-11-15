import os

from neo4j_graph import Graph
import hgnc
import umls
import drugbank
import metamap
import gwas_catalog
import nature_predict
import frimp
import omictools

class OmicNet(Graph):
    
    node_kinds = ['drug', 'disease', 'gene']
    rel_kinds = ['disease2gene_gwas', 'disease2gene_expr', 'drug2gene_target',
                 'drug2gene_cmap', 'disease2drug_predict', 'gene2gene_imp']
    
    def __init__(self, neo_db_dir):
        """neo_db_dir contains the directory storing the neo4j database"""
        Graph.__init__(self, neo_db_dir)
        
    def create_genes(self, hgnc_dir):
        """Create gene nodes from HGNC."""
        print ' Creating Gene Nodes '.center(79, '#')
        index = ['symbol', 'previous_symbols', 'synonyms', 'ensembl_id', 'entrez_id']
        hugu = hgnc.HGNC(hgnc_dir)
        for gene in hugu.gene_generator():
            self.create_node('gene', gene, index)
    
    def create_drugs(self, drugbank_dir):
        """Create drug nodes from a drugbank download."""
        print ' Creating Drug Nodes '.center(79, '#')
        dbank = drugbank.DrugBank(drugbank_dir)
        dbank.read(parse_partners=False, parse_targets=False)
        index = ['drugbank_id', 'name', 'brands', 'synonyms']
        for drug in dbank.drugs:
            self.create_node('drug', drug, index)
    
    def create_diseases(self, meta_dir):
        """Create disease nodes from the UMLS."""
        print ' Creating Disease Nodes '.center(79, '#')
        metathesaurus = umls.UMLS(meta_dir)
        diseases = metathesaurus.read()
        index = ['concept_id']
        for disease_dict in diseases:
            self.create_node('disease', disease_dict, index)

    def create_disease2gene_gwas(self, current_dir, previous_dir=None):
        """Create disease2gene_gwas relationships."""
        print ' Creating disease2gene_gwas Relationships '.center(79, '#')
        gwas = gwas_catalog.GwasCatalog(current_dir, previous_dir)
        concept_id_to_genes = gwas.read()
        for concept_id, genes in concept_id_to_genes.items():
            for gene in genes:
                # Create concept_id to gene relationship
                gene_node = self.retrieve_gene(gene)
                disease_node = self.retrieve_disease(concept_id)
                if gene_node and disease_node:
                    self.create_relationship(gene_node, disease_node, 'disease2gene_gwas')

    def create_drug2gene_target(self, drugbank_dir):
        """Create drug2gene_target relations."""
        print ' Creating drug2gene_target Relationships '.center(79, '#')
        dbank = drugbank.DrugBank(drugbank_dir)
        dbank.read(parse_drugs=False)
        id_to_partner = dbank.get_id_to_partner()
        for target in dbank.targets:
            # Only take targets with a known action (Pharmacological action)
            if target['known_action'] != 'yes':
                continue
            drug_node = self.retrieve_drug(target['drugbank_id'], drugbank_id=True)
            partner_id = target.pop('partner')
            del target['drugbank_id']
            partner = id_to_partner[partner_id]
            if 'gene_name' not in partner:
                continue
            if partner['species'] != 'Homo sapiens':
                continue
            gene_node = self.retrieve_gene(partner['gene_name'])
            if drug_node and gene_node:
                self.create_relationship(drug_node, gene_node, 'drug2gene_target', target)
    
    
    def create_disease2gene_expr(self, gxa_dir, previous_gxa_dir):
        """GXA disease2gene_expr"""
        print ' Creating disease2gene_expr Relationships '.center(79, '#')

    
    def create_disease2drug_predict(self, predict_dir):
        """Create disease2drug_predict relationships."""
        print ' Creating disease2drug_predict Relationships '.center(79, '#')
        np = nature_predict.NaturePredict(predict_dir)
        concept_id_to_indicated_drugs = np.read()
        for concept_id, indicated_drugs in concept_id_to_indicated_drugs.items():
            for drug in indicated_drugs:
                disease_node = self.retrieve_disease(concept_id)
                drug_node = self.retrieve_drug(drug)
                if disease_node and drug_node: 
                    self.create_relationship(disease_node, drug_node,
                                             'disease2drug_predict')

    def create_gene2gene_imp(self, imp_dir):
        """Create gene2gene_imp relationships."""
        print ' Creating gene2gene_imp Relationships '.center(79, '#')
        imp = frimp.IMP(imp_dir)
        imp_generator = imp.read(prob_cutoff = 0.5)
        for gene_0, gene_1, prob in imp_generator:
            gene_0_node = self.retrieve_gene(gene_0, entrez=True)
            gene_1_node = self.retrieve_gene(gene_1, entrez=True)
            if gene_0_node and gene_1_node: 
                self.create_relationship(gene_0_node, gene_1_node,
                                         'gene2gene_imp', {'probability': prob})

    def retrieve_disease(self, concept_id):
        """Returns the disease node with the specified UMLS concept_id."""
        return self.query_node('kind:disease AND concept_id:' + concept_id)
    
    def retrieve_gene(self, name, entrez=False):
        """Retreive a gene by hugu gene name (default) or ensembl id."""
        name = name.replace(' ', '\\ ')
        if entrez:
            return self.query_node('kind:gene AND entrez_id:' + name)
        query_list = ['kind:gene AND symbol:' + name,
                      'kind:gene AND synonyms:' + name,
                      'kind:gene AND previous_symbols:' + name]
        for query in query_list:
            node = self.query_node(query)
            if node: return node
        return None

    def retrieve_drug(self, name, drugbank_id=False):
        """Retreive a drug by a name (default) or drugbank_id."""
        if drugbank_id:
            return self.query_node('kind:drug AND drugbank_id:' + name)
        
        name = name.replace(' ', '\\ ')
        query_list = ['kind:drug AND name:' + name,
                      'kind:drug AND brands:' + name,
                      'kind:drug AND synonyms:' + name]
        for query in query_list:
            node = self.query_node(query)
            if node: return node
        return None

    def create_subgraph(self, node_kinds, rel_kinds, sub_db_dir):
        """ """
        sub_g = OmicNet(sub_db_dir)
        super(OmicNet, self).create_subgraph(node_kinds, rel_kinds, sub_g)
        return sub_g
        
        
