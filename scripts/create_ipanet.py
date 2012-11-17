import networkx

import bioparser.data

def create_ipanet():
    hgnc = bioparser.data.Data().hgnc
    entrez_to_hgnc = hgnc.get_entrez_to_gene()
    symbol_to_hgnc = hgnc.get_symbol_to_gene()
    
    ipa = bioparser.data.Data().ipa
    ipa.build()
    
    g = networkx.Graph(name='ipanet')
    
    ################################################################################
    ################################# Create Nodes #################################
    
    # Create drug nodes
    targets = set()
    for drug in ipa.drugs:
        g.add_node(drug.symbol, kind='drug')
        targets |= set(drug.targets)
    
    # Create disease nodes
    for disease in ipa.functions:
        g.add_node(disease.name, kind='disease')
    
    # Create gene nodes
    hugu_genes_added = set()
    for gene in ipa.genes:
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
    for drug in ipa.drugs:
        for target in drug.targets:
            g.add_edge(drug.symbol, target)
    
    # Create disease-gene and disease-drug links from ipa function annotations.
    for disease in ipa.functions:
        for effect, molecules in disease.molecules.items():
            for molecule in molecules:
                if disease.name in g or molecule in g:
                    g.add_edge(disease.name, molecule, effect=effect)
    
    return g


################################################################################
################################# Network Stats ################################
#print 'calculating statistics'


if __name__ == '__main__':
    g = create_ipanet()
    # pagerank
    pr = networkx.pagerank(g)
    pr_sorted = sorted(pr.items(), key=lambda x: x[1])
    print pr_sorted[:5]







