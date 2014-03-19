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

def get_fdr_genes(gene_to_vegasp, cutoff=0.2):
    genes = gene_to_vegasp.keys()
    vegasp = gene_to_vegasp.values()
    fdr_p = fdr_generator(vegasp)
    fdr_genes = {gene for gene, fdr in zip(genes, fdr_p) if fdr <= cutoff}
    return fdr_genes

def get_nominal_genes(gene_to_vegasp, cutoff=0.05):
    """
    Get nominally significant genes at the cutoff threshold
    """
    genes = {gene for gene, vegasp in gene_to_vegasp.iteritems() if vegasp <= cutoff}
    return genes


def name_to_path_from_dir(directory,
                          head='meta_USC_UC_IARC_updated_',
                          tail='.hg19_vegas_results'):
    """
    Given a directory containing vegas results, return a dictionary of the
    GWAS names as keys and full file paths as values. GWAS names are extracted
    from the file names as the text in between the strings specified by the
    head and tail arguments (must include start and end of string).
    File names without a match are excluded.
    """
    filenames = os.listdir(directory)
    pattern = re.compile(r"^{}(.*?){}$".format(re.escape(head), re.escape(tail)))
    name_to_path = dict()
    for filename in filenames:
        match = re.search(pattern, filename)
        if not match:
            continue
        name = match.group(1)
        name_to_path[name] = os.path.join(directory, filename)
    return name_to_path

def genes_from_directory(vegas_dir, method=get_nominal_genes, cutoff=0.05):
    """
    To FDR adjust set method equal to vegas_reader.get_fdr_genes.
    """
    disease_to_genes = dict()
    name_to_path = name_to_path_from_dir(vegas_dir)
    for name, path in name_to_path.items():
        gene_to_vegasp = read_vegas(path)
        genes = method(gene_to_vegasp, cutoff)
        disease_to_genes[name] = genes
    return disease_to_genes

if __name__ == '__main__':
    hodgkin_dir = '/home/dhimmels/Documents/serg/hodgkins/'
    vegas_dir = os.path.join(hodgkin_dir, 'input', 'vegas')
    disease_to_genes = genes_from_directory(vegas_dir, get_fdr_genes, 0.1)
    import pprint
    pprint.pprint(disease_to_genes)