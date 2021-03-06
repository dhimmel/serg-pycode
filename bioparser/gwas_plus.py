import os
import csv
import collections
import re
import logging
import datetime
import operator
import math

import data


def wabf_from_theta(estimate, variance, prior_variance):
    """
    From doi:10.1038/nrg2615 (page 683, formula 6)
    """
    total_variance = variance + prior_variance
    z_score = estimate / math.sqrt(variance)
    try:
        return math.sqrt(variance / total_variance) * math.exp(
            prior_variance * z_score ** 2 / (2 * total_variance))
    except OverflowError:
        return float('inf')


def theta_variance_from_pvalue(theta_hat, mlog10_pvalue):
    log10_pvalue = -mlog10_pvalue
    log_pvalue = log10_pvalue / math.log10(math.e)
    log_pvalue -= math.log(2)
    import rpy2.robjects
    qnorm = rpy2.robjects.r['qnorm']
    z_score = qnorm(p=log_pvalue, log_p=True)[0]
    theta_hat_sd = theta_hat / z_score
    return theta_hat_sd ** 2

def theta_variance_from_pvalue_old(theta_hat, pvalue):
    import scipy.stats
    z_score = abs(scipy.stats.norm.ppf(pvalue / 2.0))
    print z_score
    theta_hat_sd = theta_hat / z_score
    return theta_hat_sd ** 2


def calculate_posterior(bayes_factor, prior_prob):
    """
    Calculate the posterior probability of association from the bayes factor
    and prior probability of association.
    """
    if bayes_factor == float('inf'):
        return 1.0
    posterior_odds = bayes_factor * prior_prob / (1.0 - prior_prob)
    return posterior_odds / (1.0 + posterior_odds)


def posterior_probability(odds_ratio, mlog10_pvalue, prior_sd=0.2, prior_prob=1e-6):
    """
    Returns the posterior probability of association for a GWAS association
    calculated from the association odds_ratio and pvalue. Priors are
    required for the standard deviation of the true effect size distribution
    and prior probability of association.
    """
    estimate = math.log(odds_ratio)
    variance = theta_variance_from_pvalue(estimate, mlog10_pvalue)
    wabf = wabf_from_theta(estimate, variance, prior_sd ** 2)
    return calculate_posterior(wabf, prior_prob)


class GwasCatalog(object):

    def __init__(self, directory=None, processed_dirname='processed', ignore_pmids=set()):
        if directory is None:
            directory = data.current_path('gwas-catalog', require_dated_format=True)
        self.directory = directory
        self.processed_dir = os.path.join(directory, processed_dirname)
        if not os.path.isdir(self.processed_dir):
            os.mkdir(self.processed_dir)
        self.ignore_pmids = ignore_pmids

    def get_catalog_rows(self, path=None):
        if hasattr(self, 'catalog_rows'):
            return self.catalog_rows

        if not path:
            path = os.path.join(self.directory, 'gwascatalog.txt')

        #symbol_to_gene = data.Data().hgnc.get_symbol_to_gene()
        entrez_to_gene = data.Data().hgnc.get_entrez_to_gene()
        symbols_to_genes = data.Data().hgnc.identifiers_to_genes

        # Parsing Parameters
        rsid_pattern = re.compile(r"rs[0-9]+")
        reported_gene_split_pattern = re.compile(r"[, ]+")
        mapped_gene_split_pattern = re.compile(r" - |;")
        invalid_symbols = {'Intergenic', 'NR', 'Pending'}
        date_format = '%m/%d/%Y'
        snp_number_pattern = re.compile(r"[[({][^\d]*([\d,.]+)[^\d]*[\])};]")
        non_decimal_pattern = re.compile(r"[^\d.]+")
        odds_ratio_ci_pattern = re.compile(r"([\d.]+)-([\d.]+)")
        case_pattern = r'([\d,]+)[^\d]*?cases'
        control_pattern = r'([\d,]+)[^\d]*?controls'
        sample_size_pattern = r'[\d,]+'



        read_file = open(path)
        reader = csv.DictReader(read_file, delimiter='\t')
        catalog_rows = list()
        for row in reader:
            parsed_row = collections.OrderedDict()
            parsed_row['phenotype'] = row['Disease/Trait']
            parsed_row['rsids'] = re.findall(rsid_pattern, row['SNPs'])
            parsed_row['date'] = datetime.datetime.strptime(row['Date'], date_format)
            reported_symbols = re.split(reported_gene_split_pattern, row['Reported Gene(s)'])
            reported_symbols = set(reported_symbols) - invalid_symbols
            reported_genes = symbols_to_genes(reported_symbols, coding=False)
            parsed_row['reported_symbols'] = reported_symbols
            parsed_row['reported_genes'] = reported_genes

            #mapped_symbols_field = row['Mapped_gene']
            #mapped_symbols = re.split(mapped_gene_split_pattern, mapped_symbols_field)
            #mapped_symbols.remove('')
            #mapped_genes = symbols_to_genes(mapped_symbols, coding=False)
            #parsed_row['mapped_genes'] = mapped_genes
            #if ' - ' not in mapped_symbols_field and len(mapped_genes) == 1:
            #    parsed_row['mapped_gene'] = mapped_genes,

            mapped_entrez = set(row['Snp_gene_ids'].split(';'))
            mapped_entrez.discard('')
            mapped_genes = symbols_to_genes(mapped_entrez, id_type='entrez', coding=False)
            parsed_row['mapped_genes'] = mapped_genes

            upstream_gene = entrez_to_gene.get(row['Upstream_gene_id'])
            if upstream_gene:
                parsed_row['upstream_gene'] = upstream_gene
            downstream_gene = entrez_to_gene.get(row['Downstream_gene_id'])
            if downstream_gene:
                parsed_row['downstream_gene'] = downstream_gene

            parsed_row['pubmed'] = row['PUBMEDID']
            if parsed_row['pubmed'] in self.ignore_pmids:
                continue
            parsed_row['chromosome'] = row['Chr_id']
            parsed_row['position'] = None if row['Chr_pos'] == '' else int(row['Chr_pos'])
            parsed_row['mlog_pval'] = None if row['Pvalue_mlog'] == '' else float(row['Pvalue_mlog'])

            # Parse minor allele frequency
            try:
                raf = float(row['Risk Allele Frequency'])
                maf = min(raf, 1.0 - raf)
                parsed_row['MAF'] = maf
            except ValueError:
                pass


            # Parse Odds Ratio
            try:
                parsed_row['odds_ratio'] = float(row['OR or beta'])
            except ValueError:
                pass

            # Parse Odds Ratio confidence interval
            try:
                match = re.search(odds_ratio_ci_pattern, row['95% CI (text)'])
                match_lower, match_upper = match.group(1), match.group(2)
                ci_tuple = float(match_lower), float(match_upper)
                parsed_row['odds_ratio_ci'] = ci_tuple
            except (IndexError, AttributeError, ValueError):
                pass

            # parse the number of SNPs assayed that passed QC
            try:
                snp_number_field = row['Platform [SNPs passing QC]']
                match = re.search(snp_number_pattern, snp_number_field)
                snp_number = match.group(1)
                snp_number = float(snp_number.replace(',', ''))
                if 'million' in snp_number_field:
                    snp_number *= 1e6
                parsed_row['snp_number'] = int(snp_number)
            except (IndexError, AttributeError, ValueError):
                pass

            # parse case and control number
            discovery_size = row['Initial Sample Size']
            for key, pattern in [('cases', case_pattern), ('controls', control_pattern)]:
                try:
                    match_str = re.search(pattern, discovery_size).group(1)
                    parsed_row[key] = int(match_str.replace(',', ''))
                except (IndexError, AttributeError, ValueError):
                    pass
            if 'cases' not in parsed_row and 'controls' not in parsed_row:
                try:
                    match_str = re.search(sample_size_pattern, discovery_size).group(0)
                    sample_size = int(match_str.replace(',', ''))
                    parsed_row['sample_size'] = sample_size
                    parsed_row['merged_sample_size'] = sample_size
                except (IndexError, AttributeError, ValueError):
                    parsed_row['merged_sample_size'] = 0
            else:
                parsed_row['merged_sample_size'] = parsed_row.get('cases', 0) + parsed_row.get('controls', 0)

            catalog_rows.append(parsed_row)
        read_file.close()

        self.catalog_rows = catalog_rows
        return catalog_rows


    def get_filtered_rows(self, mlog_pval_threshold=None):
        """

        """
        if hasattr(self, 'filtered_rows'):
            return self.filtered_rows()

        filtered_rows = list()
        catalog_rows = self.get_catalog_rows()
        for catalog_row in catalog_rows:
            try:
                catalog_row['rsid'],  = catalog_row['rsids']
            except ValueError:
                # association not defined by a single lead SNP
                continue
            if mlog_pval_threshold is not None and catalog_row['mlog_pval'] < mlog_pval_threshold:
                continue
            filtered_rows.append(catalog_row)

        self.filtered_rows = filtered_rows
        return filtered_rows

    def annotate_snap_wingspans(self):
        wingspan_path = os.path.join(self.directory, 'SNAP', 'wingspans_0.1-centimorgans.txt')
        with open(wingspan_path) as wingspan_file:
            reader = csv.DictReader(wingspan_file, delimiter='\t')
            rsid_to_wingspan = {wingspan['rsid']: wingspan for wingspan in reader}

        for catalog_row in self.get_catalog_rows():
            try:
                rsid, = catalog_row['rsids']
            except ValueError:
                continue
            wingspan = rsid_to_wingspan.get(rsid)
            if not wingspan:
                continue
            catalog_row['snap_wingspan'] = wingspan['lower'], wingspan['upper']

    def annotate_dapple_genes(self):
        wingspan_path = os.path.join(self.directory, 'dapple', 'WSinput.txt')
        with open(wingspan_path) as wingspan_file:
            fieldnames = 'rsid', 'chromosome', 'lower', 'upper'
            reader = csv.DictReader(wingspan_file, delimiter='\t', fieldnames=fieldnames)
            rsid_to_wingspan = {wingspan['rsid']: wingspan for wingspan in reader}

        genelist_path = os.path.join(self.directory, 'dapple', 'genelist.tmp')
        rsid_to_genes = dict()
        ensembl_to_gene = data.Data().hgnc.get_ensembl_to_gene()
        symbol_to_gene = data.Data().hgnc.get_symbol_to_gene()
        with open(genelist_path) as genelist_file:
            reader = csv.reader(genelist_file, delimiter='\t')
            for row in reader:
                gene = ensembl_to_gene.get(row[1])
                if not gene:
                    gene = symbol_to_gene.get(row[11])
                if not gene:
                    # cannot map gene to HGNC
                    continue
                rsid_to_genes.setdefault(row[17], list()).append(gene)

        for catalog_row in self.get_catalog_rows():
            try:
                rsid, = catalog_row['rsids']
            except ValueError:
                continue
            wingspan = rsid_to_wingspan.get(rsid)
            if not wingspan:
                continue
            catalog_row['dapple_wingspan'] = wingspan['lower'], wingspan['upper']
            catalog_row['dapple_genes'] = rsid_to_genes.get(rsid, list())

    def annotate_efo(self, path=None):
        """
        Read the mapping file available at the EBI's GWAS Diagram Browser:
        http://www.ebi.ac.uk/fgpt/gwas/#downloadstab

        association receives an item with key 'efo_id' only if one efo_id is mapped
        to the association.
        """
        if path is None:
            file_names = os.listdir(self.directory)
            file_names = filter(lambda s: re.match(r"GWAS-EFO-Mappings[0-9\-]*\.txt$", s), file_names)
            file_names.sort()
            file_name = file_names[-1]
            path = os.path.join(self.directory, file_name)

        efo_graph = data.Data().efo.get_graph()
        trait_tuple_to_efo_ids = dict()
        read_file = open(path)
        unmatched_efo_terms = set()
        reader = csv.DictReader(read_file, delimiter='\t')
        for row in reader:
            catalog_term = row['DISEASETRAIT']
            efo_id = row['EFOURI']
            efo_id = efo_id.rsplit('/', 1)[-1]
            efo_id = efo_id.replace('rdfns#', '')
            efo_id = efo_id.replace('CL#', '')
            efo_id = ':'.join(efo_id.rsplit('_', 1))
            if efo_id not in efo_graph.node:
                unmatched_efo_terms.add(efo_id)
                continue
            trait_tuple = row['PUBMEDID'], catalog_term
            trait_tuple_to_efo_ids.setdefault(trait_tuple, list()).append(efo_id)
        for efo_id in unmatched_efo_terms:
            logging.warning(efo_id + " from EBI's gwas_catalog_to_EFO_mappings not found in EFO.")
        read_file.close()

        # Trait tuple is (pubmed_id, trait)
        trait_tuple_to_rows = dict()
        for association in self.get_catalog_rows():
            trait_tuple = association['pubmed'], association['phenotype']
            trait_tuple_to_rows.setdefault(trait_tuple, list()).append(association)

        for trait_tuple, efo_ids in trait_tuple_to_efo_ids.iteritems():
            for association in trait_tuple_to_rows.get(trait_tuple, list()):
                if len(efo_ids) == 1:
                    association['efo_id'] = efo_ids[0]
                else:
                    association['efo_ids'] = efo_ids

    def annotate_doid(self, ontoprocess_path=None):
        self.annotate_efo()
        doid_graph = data.Data().doid.get_graph()
        efo_to_doid_ids = data.Data().doid.get_xref_to_doids('EFO', 'EFO:')

        remove_doids, pop_doids = self.read_ontprocess_info(ontoprocess_path)
        for association in self.get_catalog_rows():
            efo_id = association.get('efo_id')
            if not efo_id:
                continue
            doid_ids = efo_to_doid_ids.get(efo_id)
            if not doid_ids:
                continue
            if len(doid_ids) > 1:
                print association, doid_ids
                continue
            doid_id, = doid_ids
            if doid_id in remove_doids:
                continue
            doid_id = pop_doids.get(doid_id, doid_id)
            association['doid_id'] = doid_id
            association['doid_name'] = doid_graph.node[doid_id]['name']


    def write_all_snps(self):
        """Write a tab delimited file of SNPs. columns represent rsid and chromosome."""
        all_snps = set()
        for row in self.get_catalog_rows():
            chromosome = row['chromosome']
            if not chromosome:
                continue
            chromosome = int(chromosome)
            all_snps |= {(snp, chromosome) for snp in row['rsids']}
        all_snps = sorted(all_snps, key=lambda x: (x[1], x[0]))
        path = os.path.join(self.directory, 'all_snps.txt')
        with open(path, 'w') as write_file:
            write_file.write('\n'.join('{}\t{}'.format(*snp) for snp in all_snps))

    @staticmethod
    def max_counter_keys(counter, key_subset=None):
        if key_subset is not None:
            counter = collections.Counter({k: v for
                k, v in counter.iteritems() if k in key_subset})
        max_value = max(counter.values() + [0])
        max_keys = [k for k, v in counter.iteritems() if v == max_value]
        return max_keys

    def get_merged_associations(self, doidprocess_path=None, wingspan_key='dapple_wingspan',
                                mlog_cutoff = -math.log10(5e-8), sample_size_cutoff=1000):
        """
        association is a list with the first association representing the
        study reporting the strongest association for that region.
        Default negative log10 p-value threshold = math.log10(5e-8)
        """
        self.annotate_doid(ontoprocess_path=doidprocess_path)
        self.annotate_dapple_genes()
        self.annotate_snap_wingspans()
        associations = self.get_filtered_rows()
        associations = [a for a in associations if a.get('doid_id')]
        associations = [a for a in associations if a.get(wingspan_key)]

        # Calculate posterior probabilities of association
        for a in associations:
            try:
                a['ppa'] = posterior_probability(
                    odds_ratio=a['odds_ratio'], mlog10_pvalue=a['mlog_pval'],
                    prior_sd=0.2, prior_prob=1e-5)
            except (KeyError, ValueError):
                pass
        # Filters associations without posterior probabilities of association
        # associations = [a for a in associations if 'ppa' in a]

        for key in ['MAF', 'ppa', 'odds_ratio', 'odds_ratio_ci', 'snp_number', 'cases', 'controls', 'sample_size']:
            print key, sum(key in a for a in associations) * 1. / len(associations)

        path = os.path.join(self.processed_dir, 'catalog-rows.txt')
        write_file = open(path, 'w')
        fieldnames = ['doid_id', 'doid_name', 'rsid', 'pubmed', 'mlog_pval',
                      'odds_ratio', 'ppa', 'cases', 'controls', 'sample_size',
                      'merged_sample_size', 'snp_number']
        writer = csv.DictWriter(write_file, fieldnames=fieldnames,
                                delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(associations)
        write_file.close()

        for a in associations:
            a['is_HC'] = (a['mlog_pval'] >= mlog_cutoff and
                          a.get('merged_sample_size', 0) >= sample_size_cutoff)

        associations.sort(key=lambda a: (a['is_HC'], a['mlog_pval'], a['date']), reverse=True)
        doid_to_associations = dict()
        for association in associations:
            doid_id = association['doid_id']
            chromosome = association['chromosome']
            wingspan = association[wingspan_key]
            included_associations = doid_to_associations.setdefault(doid_id, list())
            add = True
            for included_association in included_associations:
                incl_wing = included_association[0][wingspan_key]
                included_chromosome = included_association[0]['chromosome']
                if (chromosome == included_chromosome and
                    (incl_wing[0] <= wingspan[0] and wingspan[0] <= incl_wing[1]) or
                    (incl_wing[0] <= wingspan[1] and wingspan[1] <= incl_wing[1])):
                    # overlapping association with greater significance already included
                    add = False
                    included_association.append(association)
                    break
            if add:
                included_associations.append([association])

        merged_associations = list()
        disease_to_statuses = dict()
        for doid, associations in doid_to_associations.items():
            for association in associations:

                # Create a counter of gene reports
                report_counter = collections.Counter()
                for a in association:
                    report_counter.update(a['reported_genes'])
                for gene in report_counter.keys():
                    if not gene.coding:
                        del report_counter[gene]
                mode_reported = GwasCatalog.max_counter_keys(report_counter)

                # Mapped genes
                all_mapped_genes = set()
                sequentially_mapped_gene = None

                for a in association:
                    mapped_genes = filter(operator.attrgetter('coding'), a.get('mapped_genes', []))
                    upstream_gene = a.get('upstream_gene')
                    downstream_gene = a.get('downstream_gene')
                    stream_genes = {upstream_gene, downstream_gene}
                    stream_genes.discard(None)
                    stream_genes = set(filter(operator.attrgetter('coding'), stream_genes))
                    mapped_gene = None
                    # GWAS catalog finds that rsid is within a single gene
                    if len(mapped_genes) == 1:
                        mapped_gene, = mapped_genes
                    # GWAS catalog finds that rsid is within multiple genes
                    elif len(mapped_genes) > 1:
                        mode_mapped_genes = GwasCatalog.max_counter_keys(report_counter, key_subset=mapped_genes)
                        if len(mode_mapped_genes) == 1:
                            mapped_gene, = mode_mapped_genes
                    # GWAS catalog finds an upstream and downstream gene
                    elif len(mapped_genes) == 0:
                        mode_stream_genes = GwasCatalog.max_counter_keys(report_counter, key_subset=stream_genes)
                        if len(mode_stream_genes) == 1:
                            mapped_gene, = mode_stream_genes

                    a['mapped_gene'] = mapped_gene
                    if not sequentially_mapped_gene and mapped_gene:
                        sequentially_mapped_gene = mapped_gene

                    all_mapped_genes |= set(mapped_genes)
                    all_mapped_genes |= stream_genes


                candidate_genes = all_mapped_genes | set(report_counter)
                resolved_gene = mode_reported[0] if len(mode_reported) == 1 else sequentially_mapped_gene

                high_confidence = any(a['is_HC'] for a in association)

                merged = {'disease_code': doid, 'disease_name': association[0]['doid_name'],
                          'resolved_gene': resolved_gene,
                          'candidate_genes': '|'.join(map(str, candidate_genes)),
                          'report_counter': report_counter,
                          'mapped_genes': '|'.join(str(a['mapped_gene']) for a in association),
                          'dapple_genes': '|'.join(map(str, association[0].get('dapple_genes', []))),
                          'studies': '|'.join(a['pubmed'] for a in association),
                          'snps': '|'.join(a['rsid'] for a in association),
                          'mlog_pvals': '|'.join('{0:.3f}'.format(a['mlog_pval']) for a in association),
                          'confidence': 'HC' if high_confidence else 'LC'}
                merged_associations.append(merged)

                disease_tuple = merged['disease_code'], merged['disease_name']
                statuses = disease_to_statuses.setdefault(disease_tuple, dict())
                confidence = merged['confidence']
                primary_set = statuses.setdefault('{}_primary'.format(confidence), set())
                secondary_set = statuses.setdefault('{}_secondary'.format(confidence), set())
                primary_set.add(resolved_gene)
                secondary_set.update(candidate_genes)

        status_rows = list()
        for (disease_code, disease_name), statuses in disease_to_statuses.items():
            HC_primary = statuses.get('HC_primary', set())
            HC_secondary = statuses.get('HC_secondary', set())
            LC_primary = statuses.get('LC_primary', set())
            LC_secondary = statuses.get('LC_secondary', set())
            for gene_set in [HC_primary, HC_secondary, LC_primary, LC_secondary]:
                gene_set.discard(None)

            # pecking order: HC_primary, HC_secondary, LC_primary, LC_secondary
            LC_secondary.difference_update(HC_primary | HC_secondary | LC_primary)
            LC_primary.difference_update(HC_primary | HC_secondary)
            HC_secondary.difference_update(HC_primary)
            for status_key, gene_set in statuses.items():
                for gene in gene_set:
                    status_row = {'disease_code': disease_code, 'disease_name': disease_name,
                     'gene_code': gene.hgnc_id, 'gene_symbol': gene.symbol,
                     'status': status_key}
                    status_rows.append(status_row)
        status_rows.sort(key=operator.itemgetter('disease_name', 'status', 'gene_symbol'))
        assert len(status_rows) == len(set(map(operator.itemgetter('disease_code', 'gene_code'), status_rows)))


        #write files
        fieldnames = ['disease_code', 'disease_name', 'resolved_gene',
                      'confidence', 'candidate_genes', 'report_counter',
                      'mapped_genes', 'dapple_genes', 'studies', 'snps', 'mlog_pvals']
        path = os.path.join(self.processed_dir, 'associations-processing.txt')
        write_file = open(path, 'w')
        writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(merged_associations)
        write_file.close()

        # write statuses
        fieldnames = ['disease_code', 'disease_name', 'gene_code',
                      'gene_symbol', 'status']
        path = os.path.join(self.processed_dir, 'association-statuses.txt')
        write_file = open(path, 'w')
        writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(status_rows)
        write_file.close()


        # Compute the number of associations per disease
        disease_to_counter = dict()
        for status_row in status_rows:
            disease = status_row['disease_code'], status_row['disease_name']
            status = status_row['status']
            counter = disease_to_counter.setdefault(disease, collections.Counter())
            counter[status] += 1

        disease_count_rows = list()
        for (disease_code, disease_name), counter in disease_to_counter.items():
            row = {'disease_code': disease_code, 'disease_name': disease_name,
                   'HC_primary': counter['HC_primary'],
                   'HC_secondary': counter['HC_secondary'],
                   'LC_primary': counter['LC_primary'],
                   'LC_secondary': counter['LC_secondary'],
                   'excluded': sum(counter.values()) - counter['HC_primary']}
            disease_count_rows.append(row)
        disease_count_rows.sort(key=operator.itemgetter('HC_primary', 'LC_primary'), reverse=True)

        fieldnames = ['disease_code', 'disease_name', 'HC_primary', 'HC_secondary', 'LC_primary', 'LC_secondary', 'excluded']
        path = os.path.join(self.processed_dir, 'associations-per-disease.txt')
        with open(path, 'w') as write_file:
            writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(disease_count_rows)

        return merged_associations


    @staticmethod
    def read_ontprocess_info(path):
        """
        Terms to omit should be encoded like:
        remove TERM_ID # comment can go here

        Terms to omit with annotations transferred to a second node:
        pop TERM_ID --> TERM_ID # comment can go here
        """
        if path is None:
            return set(), dict()

        with open(path) as read_file:
            lines = list(read_file)

        remove_terms = set()
        pop_terms = dict()
        for line in lines:
            line = line.split('#', 1)[0].rstrip()
            line_list = line.split(' ')
            command = line_list.pop(0)
            if command == 'remove':
                remove_terms.add(line_list[0])
            elif command == 'pop':
                pop_terms[line_list[0]] = line_list[2]
            else:
                assert False
        return remove_terms, pop_terms



if __name__ =='__main__':
    import pprint

    doidprocess_path = '/home/dhimmels/Documents/serg/gene-disease-hetnet/data-integration/doid-ontprocess-info.txt'

    gcat = GwasCatalog()
    merged_associations = gcat.get_merged_associations(doidprocess_path=doidprocess_path)

    gcat_no_wtccc2 = GwasCatalog(processed_dirname='processed-no-wtccc2', ignore_pmids={'21833088'})
    merged_associations_no_wtccc2 = gcat_no_wtccc2.get_merged_associations(doidprocess_path=doidprocess_path)


    #print len(gcat.get_merged_associations('snap_wingspan'))
    #print sum(len(v) for v in gcat.get_merged_associations('snap_wingspan').values())

    #gcat.annotate_dapple_genes()
    #gcat.annotate_snap_wingspans()
    #rows = gcat.get_filtered_rows()
    #print 'Total', len(rows)
    #print 'dapple', sum('dapple_wingspan' in row for row in rows)
    #print 'snap', sum('snap_wingspan' in row for row in rows)

    #catalog_rows = gcat.get_catalog_rows()
    #gcat.annotate_dapple_genes()
    #filtered_rows = gcat.get_filtered_rows()
    #sum(int(bool(row.get('dapple_wingspan'))) for row in filtered_rows)
    #gcat.write_all_snps()
    #merged_associations = gcat.get_merged_associations()


    """
    disease_to_genes = dict()
    for a in merged_associations:
        genes = a['gene_list']
        if not len(genes) == 1:
            continue
        gene, = genes
        if gene.locus_group != 'protein-coding gene':
            print 'Noncoding', gene
            continue
        disease_to_genes.setdefault(a['doid_name'], set()).add(gene)
    pprint.pprint(disease_to_genes)
    """
    #gcat_no_wtccc2 = GwasCatalog(processed_dirname='processed-no-wtccc2', ignore_pmids={'21833088'})
    #merged_associations_no_wtccc2 = gcat_no_wtccc2.get_doid_to_associations()




    #associations = reduce(operator.add, doid_to_associations.values())
    #print collections.Counter([len(a['reported_genes']) for a in associations])
    #pprint.pprint(doid_to_associations['DOID:2377'])