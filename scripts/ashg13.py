import bioparser.data


data = bioparser.data.Data()
gcat = data.gwas_catalog

print gcat.get_efo_id_to_genes(fdr_cutoff=0.05)