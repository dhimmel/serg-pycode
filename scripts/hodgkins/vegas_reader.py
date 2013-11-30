import re
import os
import csv

import bioparser.data



def fdr_generator(p_values, num_total_tests=None):
    """
    Calculates the Benjamini-Hochberg correction for multiple hypothesis
    testing from a list of p-values *sorted in ascending order*.

    See
    http://en.wikipedia.org/wiki/False_discovery_rate#Independent_tests
    for more detail on the theory behind the correction.

    **NOTE:** This is a generator, not a function. It will yield values
    until all calculations have completed.

    :Parameters:
    - 'p_values': a list or iterable of p-values sorted in ascending
      order
    - 'num_total_tests': the total number of tests (p-values)

    """
    if num_total_tests is None:
        num_total_tests = len(p_values)
    p_values = sorted(p_values)
    prev_bh_value = 0
    for i, p_value in enumerate(p_values):
        bh_value = p_value * num_total_tests / (i + 1)
        # Sometimes this correction can give values greater than 1,
        # so we set those values at 1
        bh_value = min(bh_value, 1)

        # To preserve monotonicity in the values, we take the
        # maximum of the previous value or this one, so that we
        # don't yield a value less than the previous.
        bh_value = max(bh_value, prev_bh_value)
        prev_bh_value = bh_value
        yield bh_value


def read_vegas(path):
    symbol_to_gene = bioparser.data.Data().hgnc.get_symbol_to_gene()
    read_file = open(path)
    reader = csv.DictReader(read_file, delimiter='\t')
    gene_to_vegasp = dict()
    for row in reader:
        vegas_gene = row['Gene']
        gene = symbol_to_gene.get(vegas_gene)
        if not gene:
            continue
        #assert gene not in gene_to_vegasp
        gene_to_vegasp[gene] = float(row['Pvalue'])
    read_file.close()
    return gene_to_vegasp

def get_fdr_genes(path, cutoff=0.05):
    gene_to_vegasp = read_vegas(path)
    genes = gene_to_vegasp.keys()
    vegasp = gene_to_vegasp.values()
    fdr_p = fdr_generator(vegasp)
    fdr_genes = {gene for gene, fdr in zip(genes, fdr_p) if fdr <= cutoff}
    return fdr_genes


def genes_from_directory(vegas_dir, fdr_cutoff):
    disease_to_genes = dict()
    fnames = os.listdir(vegas_dir)
    pattern = re.compile(r"^meta_USC_UC_IARC_updated_(.*?)\.hg19_vegas_results$")
    for fname in fnames:
        hodgkin_type = re.search(pattern, fname).group(1)
        path = os.path.join(vegas_dir, fname)
        genes = get_fdr_genes(path, fdr_cutoff)
        disease_to_genes[hodgkin_type] = genes
    return disease_to_genes

if __name__ == '__main__':
    hodgkin_dir = '/home/dhimmels/Documents/serg/hodgkins/'
    vegas_dir = os.path.join(hodgkin_dir, 'vegas')
    disease_to_genes = genes_from_directory(vegas_dir, 0.2)
    import pprint
    pprint.pprint(disease_to_genes)