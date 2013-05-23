import os
import csv
import collections
import re

import data


def fdr_generator(p_values, num_total_tests):
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
                
                # Convert the p-value to a float
                try:
                    row['p-Value'] = float(row['p-Value'])
                except ValueError:
                    row['p-Value'] = None
                
                # Extract the reported genes into a set
                genes = row['Reported Gene(s)'].split(',')
                genes = set(gene.strip() for gene in genes)
                genes -= invalid_genes
                row['genes'] = genes
                
                # Extract the number of SNPs assesed. None for Not Reported
                snp_num = row['Platform [SNPs passing QC]']
                try:
                    snp_num = re.search(r"\[.*?\]", snp_num).group(0)
                    snp_num = re.search(r"(?<=\[).*?(?=[;\]])", snp_num).group(0)
                    snp_num.replace('million', ',000,000')
                    snp_num = filter(str.isdigit, snp_num)
                    snp_num = int(snp_num)
                except (AttributeError, ValueError):
                    # AttributeError when no regex matches due to 'NR' field
                    # ValueError when 'NR' (no digits) inside of brackets.
                    snp_num = None
                row['Number of SNPs'] = snp_num
                
                yield row
    
    def get_rows(self):
        if not hasattr(self, 'rows'):
            self.rows = list(self.row_generator())
        return self.rows
    
    def get_pmid_to_rows(self):
        """ """
        pmid_to_rows = dict()
        for row in self.get_rows():
            pmid_to_rows.setdefault(row['PUBMEDID'], list()).append(row)
        self.pmid_to_rows = pmid_to_rows
        return pmid_to_rows

    def get_trait_tuple_to_rows(self):
        """Trait tuple is (pubmed_id, trait)"""
        trait_tuple_to_rows = dict()
        for row in self.get_rows():
            trait_tuple = row['PUBMEDID'], row['Disease/Trait']
            trait_tuple_to_rows.setdefault(trait_tuple, list()).append(row)
        self.trait_tuple_to_rows = trait_tuple_to_rows
        return trait_tuple_to_rows
    
    def apply_fdr(self):
        """
        """
        trait_tuple_to_rows = self.get_trait_tuple_to_rows()
        for trait_tuple, rows in trait_tuple_to_rows.iteritems():
            rows.sort(key=lambda row: row['p-Value'])
            pvals = [row['p-Value'] for row in rows]
            snp_num = list({row['Number of SNPs'] for row in rows})
            if len(snp_num) != 1 or snp_num[0] is None or any(pval is None for pval in pvals):
                fdr_pvals = [None] * len(rows)
            else:
                fdr_pvals = list(fdr_generator(pvals, snp_num[0]))
            for row, fdr_pval in zip(rows, fdr_pvals):
                row['FDR p-Value'] = fdr_pval
    
    def get_phenotypes(self):
        """Returns the set of Disease/Trait terms."""
        return {row['Disease/Trait'] for row in self.get_rows()}
    
    def efo_map(self):
        gwas_terms = self.get_phenotypes()
    
    def read_ebi_mappings(self, path=None):
        """
        Read the mapping file available at the EBI's GWAS Diagram Browser:
        http://www.ebi.ac.uk/fgpt/gwas/#downloadstab
        Returns a dictionary of GWAS Catalog term to EFO IDs (set).
        """
        if path is None:
            file_names = os.listdir(self.gwas_dir)
            file_names = filter(lambda s: re.match(r"GWAS-EFO-Mappings[0-9\-]*\.txt$", s), file_names)
            file_names.sort()
            file_name = file_names[-1]
            path = os.path.join(self.gwas_dir, file_name)

        efo_graph = data.Data().efo.get_graph()
        trait_tuple_to_efo_ids = dict()

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
                    #raise KeyError # Exception can be commented out 
                    continue
                trait_tuple = row['PUBMEDID'], catalog_term
                trait_tuple_to_efo_ids.setdefault(trait_tuple, list()).append(efo_id)
        
        trait_tuple_to_rows = self.get_trait_tuple_to_rows()
        for trait_tuple, efo_ids in trait_tuple_to_efo_ids.iteritems():
            gcat_rows = trait_tuple_to_rows.get(trait_tuple, list())
            for gcat_row in gcat_rows:
                gcat_row['efo_ids'] = efo_ids
                if len(efo_ids) == 1:
                    gcat_row['efo_id'] = efo_ids[0]
                
        
    def get_efo_id_to_genes(self, p_cutoff=None, fdr_cutoff=None, mapped_term_cutoff=None):
        """ """
        self.read_ebi_mappings()
        if fdr_cutoff is not None:
            self.apply_fdr()
        efo_id_to_genes = dict()
        for row in self.get_rows():
            pval = row['p-Value']
            if p_cutoff is not None and (pval is None or pval > p_cutoff):
                continue
            if fdr_cutoff is not None:
                fdr = row['FDR p-Value']
                if fdr is None or fdr > fdr_cutoff:
                    continue
            # exclude GWAS diseases or traits that do not map to a single EFO term.
            if not 'efo_ids' in row:
                continue
            efo_ids = row['efo_ids']
            if mapped_term_cutoff is not None and len(efo_ids) > mapped_term_cutoff:
                continue
            for efo_id in efo_ids:
                efo_id_to_genes.setdefault(efo_id, set()).update(row['genes'])
        return efo_id_to_genes

if __name__ =='__main__':
    gcat = GwasCatalog()
    path = '/home/dhimmels/Documents/serg/data-sources/gwas-catalog/GWAS-EFO-Mappings092012.txt'
    gcat.read_ebi_mappings(path)
    #catalog_term_to_efo_ids, mapped_pmids = gcat.read_ebi_mappings(path)
    gcat.get_rows()
    gcat.apply_fdr()
    rows = gcat.get_rows()
    print 'SNPs passing 0.05 P-value -', sum(row['p-Value'] <= 0.05 for row in rows if row['p-Value'])
    print 'SNPs passing 0.05 P-value and reporting genes -', sum(row['p-Value'] <= 0.05 and bool(row['genes']) for row in rows if row['FDR p-Value'])
    print 'SNPs passing 0.05 FDR P-value -', sum(row['FDR p-Value'] <= 0.05 for row in rows if row['FDR p-Value'])
    print 'SNPs passing 0.05 FDR P-value and reporting genes -', sum(row['FDR p-Value'] <= 0.05 and bool(row['genes']) for row in rows if row['FDR p-Value'])
    print 'SNPs passing 0.05 FDR P-value and reporting genes and mapping to EFO -', sum(row['FDR p-Value'] <= 0.05 and bool(row['genes']) and 'efo_id' in row for row in rows if row['FDR p-Value'])
    
    efo_id_to_genes = gcat.get_efo_id_to_genes(fdr_cutoff=0.05, mapped_term_cutoff=1)
    #print efo_id_to_genes
    print len(efo_id_to_genes), 'mapped efo IDs with genes passing FDR'
    print 'Genes per EFO term counter:'
    print collections.Counter(map(len, efo_id_to_genes.values()))
    
    print 'Top 50 EFO Diseases for the number of associated GWAS genes:'
    efo_graph = data.Data().efo.get_graph()
    disease_terms = data.Data().efo.get_diseases() - data.Data().efo.get_neoplasms()
    efo_id_counter = collections.Counter({efo_graph.node[k]['name']: len(v) for k, v in efo_id_to_genes.iteritems() if k in disease_terms})
    for name, count in efo_id_counter.most_common(100):
        print name, '\t', count




