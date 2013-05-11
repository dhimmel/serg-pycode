import os
import csv
import collections

import data

class GwasCatalog(object):
    
    def __init__(self, gwas_dir=None):
        if gwas_dir is None:
            gwas_dir = data.current_path('gwas-catalog')
        self.gwas_dir = gwas_dir
        
    def row_generator(self, path=None):
        """ """
        invalid_genes = set(['Intergenic', 'NR'])
        if not path:
            path = os.path.join(self.gwas_dir, 'gwascatalog.txt')
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                genes = row['Reported Gene(s)'].split(',')
                genes = set(gene.strip() for gene in genes)
                genes -= invalid_genes
                row['genes'] = genes
                yield row
    
    def get_phenotypes(self):
        """Returns the set of Disease/Trait terms."""
        return {row['Disease/Trait'] for row in self.row_generator()}
    
    def efo_map(self):
        gwas_terms = self.get_phenotypes()
    
    def read_ebi_mappings(self, path):
        """
        Read the mapping file available at the EBI's GWAS Diagram Browser:
        http://www.ebi.ac.uk/fgpt/gwas/#downloadstab
        Returns a dictionary of GWAS Catalog term to EFO IDs (set).
        """
        efo_graph = data.Data().efo.get_graph()

        catalog_term_to_efo_ids = dict()
        mapped_pmids = set()
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                catalog_term = row['DISEASETRAIT']
                efo_id = row['EFOURI']
                efo_id = efo_id.rsplit('/', 1)[-1]
                efo_id = efo_id.replace('rdfns#', 'ORP_')
                efo_id = efo_id.replace('CL#', '')
                if efo_id not in efo_graph.node:
                    print efo_id, 'from ebi gwas catalog to EFO mappings not found in EFO.'
                    raise KeyError # Exception can be commented out 
                    continue
                catalog_term_to_efo_ids.setdefault(catalog_term, set()).add(efo_id)
                mapped_pmids.add(row['PUBMEDID'])
        print 'Distribution of number of EFO terms mapping to unique catalog terms in EBI mappings:'
        map_counter = collections.Counter(len(id_set) for id_set in catalog_term_to_efo_ids.values())
        print map_counter
        return catalog_term_to_efo_ids, mapped_pmids
        

if __name__ =='__main__':
    gcat = GwasCatalog()
    path = '/home/dhimmels/Documents/serg/data-sources/gwas-catalog/GWAS-EFO-Mappings092012.txt'
    catalog_term_to_efo_ids, mapped_pmids = gcat.read_ebi_mappings(path)
