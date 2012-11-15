import os
import metamap

class GwasCatalog(object):
    
    def __init__(self, current_dir, previous_dir=None):
        self.current_dir = current_dir
        self.previous_dir = previous_dir
        
        
    def read(self):
        """
        Path is the lacation of a gwascatalog.txt file.
        Website: http://www.genome.gov/gwastudies/
        Direct download link: http://www.genome.gov/admin/gwascatalog.txt
        """
        print 'initiating reading of the gwas catalog'
        mm = metamap.MetaMap(self.current_dir, self.previous_dir)

        disease_to_genes = dict()
        concept_id_to_genes = dict()
        invalid_genes = set(['Intergenic', 'NR'])
        self.catalog_path = os.path.join(self.current_dir, 'gwascatalog.txt')
        for line_dict in metamap.read_tdt(self.catalog_path):
            disease = line_dict['Disease/Trait']
            genes = set(gene.strip() for gene in
                        line_dict['Reported Gene(s)'].split(','))
            genes -= invalid_genes
            disease_to_genes.setdefault(disease, set()).update(genes)
        diseases = disease_to_genes.keys()

        if not mm.maps['current']['redacted']:
            mm.map_terms(diseases)
            mm.pause_before_reading_mappings()
        
        disease_to_concept_id = mm.get_mapping_dict()
        for disease, genes in disease_to_genes.items():
            if disease in disease_to_concept_id:
                concept_id = disease_to_concept_id[disease]
                concept_id_to_genes.setdefault(concept_id, set()).update(genes)

        self.concept_id_to_genes = concept_id_to_genes
        return concept_id_to_genes
