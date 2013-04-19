"""
TODO: Should node reference in indexes be deleted upon node deletion?


Find properties in GML file that will break the import
\n\t\t\w*?[^\w\s]+.*?\s[^"0-9]

gremlin query
g = new Neo4jGraph('/scratch/omicnet.db')
g.saveGML('/home/dhimmels/Documents/serg/omicnet/omicnet.gml')
g.shutdown()
"""


all_simple_paths


import omicnet
import omictools

try:
    execfile('omicnet.py')
    g = OmicNet('/home/dhimmels/Documents/serg/omicnet/omicnet.db')
    
    
    
    ###########################################################################
    #### disease2gene_gwas statistics and info
    """
    diseases = g.get_nodes_of_kind('disease')
    disease_to_gwas_genes = dict()
    for disease in diseases:
        genes = [rel.start for rel in getattr(disease, 'disease2gene_gwas')]
        if not genes:
            continue
        gene_names = [gene['symbol'] for gene in genes]
        disease_name = disease['name']
        disease_to_gwas_genes[disease_name] = gene_names
    
    #print disease_to_gwas_genes
    gene_dist = [len(genes) for genes in disease_to_gwas_genes.values()]
    all_genes = set()
    for genes in disease_to_gwas_genes.values():
        all_genes |= set(genes)
    print len(disease_to_gwas_genes), 'diseases'
    print len(all_genes), 'total associated genes'
    gene_counts = [(disease, len(set(genes))) for disease, genes in disease_to_gwas_genes.items()]
    gene_counts.sort(key=lambda x: x[1], reverse=True)
    print gene_counts
    print disease_to_gwas_genes['Multiple Sclerosis']
    """
    ###########################################################################

    ###########################################################################
    #### drug2gene_target statistics and info
    drugs = g.get_nodes_of_kind('drug')
    drug_to_targets = dict()
    print 'total drugs', len(drugs)
    for drug in drugs:
        genes = [rel.end for rel in getattr(drug, 'drug2gene_target')]
        if not genes:
            continue
        gene_names = [gene['symbol'] for gene in genes]
        drug_name = drug['name']
        drug_to_targets[drug_name] = gene_names
    all_genes = set()
    for genes in drug_to_targets.values():
        all_genes |= set(genes)
    print len(drug_to_targets), 'drugs with targets'
    print len(all_genes), 'total targets'
    gene_counts = [len(genes) for drug, genes in drug_to_targets.items()]
    gene_counts = [[drug, len(genes)] for drug, genes in drug_to_targets.items()]
    gene_counts.sort(key=lambda x: x[1], reverse=True)
    print gene_counts[:12]





finally:
    g.shutdown()
    omictools.data.base_metathesaurus.close()