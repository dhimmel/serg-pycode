import os
import csv
import collections
import re
import logging
import datetime
import operator

import data



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

        symbol_to_gene = data.Data().hgnc.get_symbol_to_gene()

        # Parsing Parameters
        rsid_pattern = re.compile(r"rs[0-9]+")
        gene_split_pattern = re.compile(r"[, ]+")
        invalid_genes = {'Intergenic', 'NR', 'Pending'}
        date_format = '%m/%d/%Y'

        read_file = open(path)
        reader = csv.DictReader(read_file, delimiter='\t')
        catalog_rows = list()
        for row in reader:
            parsed_row = collections.OrderedDict()
            parsed_row['phenotype'] = row['Disease/Trait']
            parsed_row['rsids'] = re.findall(rsid_pattern, row['SNPs'])
            parsed_row['date'] = datetime.datetime.strptime(row['Date'], date_format)
            reported_symbols = re.split(gene_split_pattern, row['Reported Gene(s)'])
            reported_genes = {symbol_to_gene.get(symbol)
                              for symbol in set(reported_symbols) - invalid_genes}
            reported_genes.discard(None)
            parsed_row['reported_symbols'] = reported_symbols
            parsed_row['reported_genes'] = reported_genes
            parsed_row['pubmed'] = row['PUBMEDID']
            if parsed_row['pubmed'] in self.ignore_pmids:
                continue
            parsed_row['chromosome'] = row['Chr_id']
            parsed_row['position'] = None if row['Chr_pos'] == '' else int(row['Chr_pos'])
            parsed_row['mlog_pval'] = None if row['Pvalue_mlog'] == '' else float(row['Pvalue_mlog'])
            catalog_rows.append(parsed_row)
        read_file.close()

        self.catalog_rows = catalog_rows
        return catalog_rows


    def get_filtered_rows(self):
        if hasattr(self, 'filtered_rows'):
            return self.filtered_rows()

        mlog_pval_threshold = 8.0

        filtered_rows = list()
        catalog_rows = self.get_catalog_rows()
        for catalog_row in catalog_rows:
            try:
                catalog_row['rsid'],  = catalog_row['rsids']
            except ValueError:
                # association not defined by a single lead SNP
                continue
            if catalog_row['mlog_pval'] < mlog_pval_threshold:
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

    def annotate_doid(self):
        self.annotate_efo()
        doid_graph = data.Data().doid.get_graph()
        efo_to_doid_ids = data.Data().doid.get_xref_to_doids('EFO', 'EFO:')
        remove_doids, pop_doids = self.read_doid_processing_info()
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

    def get_merged_associations(self, wingspan_key='dapple_wingspan'):
        """association is a list with the first association representing the
        study reporting the strongest association for that region.
        """
        self.annotate_doid()
        self.annotate_dapple_genes()
        self.annotate_snap_wingspans()
        associations = self.get_filtered_rows()
        associations = [a for a in associations if a.get('doid_id')]
        associations = [a for a in associations if a.get(wingspan_key)]
        associations.sort(key=lambda a: (a['mlog_pval'], a['date']))
        associations.reverse()
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
        for doid, associations in doid_to_associations.items():
            for association in associations:
                gene_counter = collections.Counter()
                for a in association:
                    gene_counter.update(a['reported_genes'])
                minp_reported_genes = association[0]['reported_genes']
                max_gene_count = max(gene_counter.values() + [0])
                max_genes = [k for k, v in gene_counter.items() if v == max_gene_count]
                genes = minp_reported_genes if len(minp_reported_genes) == 1 else max_genes
                merged = {'doid_code': doid, 'doid_name': association[0]['doid_name'],
                          'genes': '|'.join(map(str, genes)),
                          'gene_list': genes,
                          'gene_counter': gene_counter,
                          'dapple_genes': '|'.join(map(str, association[0].get('dapple_genes', []))),
                          'studies': '|'.join(a['pubmed'] for a in association),
                          'snps': '|'.join(a['rsid'] for a in association),
                          'mlog_pvals': '|'.join('{0:.3f}'.format(a['mlog_pval']) for a in association)}
                if len(genes) == 1:
                    merged['gene'], = genes
                merged_associations.append(merged)

        #write files
        path_singular = os.path.join(self.processed_dir, 'associations-singular.txt')
        path_multiple = os.path.join(self.processed_dir, 'associations-multiple.txt')
        fieldnames_singular = ['doid_code', 'doid_name', 'gene', 'dapple_genes', 'studies', 'snps', 'mlog_pvals']
        fieldnames_multiple = ['doid_code', 'doid_name', 'genes', 'dapple_genes', 'studies', 'snps', 'mlog_pvals']
        file_singular = open(path_singular, 'w')
        file_multiple = open(path_multiple, 'w')
        writer_singular = csv.DictWriter(file_singular, delimiter='\t', fieldnames=fieldnames_singular, extrasaction='ignore')
        writer_multiple = csv.DictWriter(file_multiple, delimiter='\t', fieldnames=fieldnames_multiple, extrasaction='ignore')
        writer_singular.writeheader()
        writer_multiple.writeheader()
        for merged_association in merged_associations:
            if 'gene' in merged_association:
                writer_singular.writerow(merged_association)
            else:
                writer_multiple.writerow(merged_association)
        file_singular.close()
        file_multiple.close()

        association_tuples = set()
        for association in merged_associations:
            if 'gene' not in association:
                continue
            gene = association['gene']
            if gene.locus_group != 'protein-coding gene':
                continue
            atup = association['doid_code'], association['doid_name'], association['gene']
            association_tuples.add(atup)

        association_tuples = sorted(association_tuples)
        path = os.path.join(self.processed_dir, 'associations.txt')
        with open(path, 'w') as write_file:
            writer = csv.writer(write_file, delimiter='\t')
            writer.writerow(['doid_code', 'doid_name', 'symbol'])
            writer.writerows(association_tuples)

        path = os.path.join(self.processed_dir, 'associations-per-disease.txt')
        disease_counter = collections.Counter((a[0], a[1]) for a in association_tuples)
        counts = list()
        for (code, name), count in disease_counter.items():
            counts.append([code, name, count])
        counts.sort(key=lambda x: -x[2])
        with open(path, 'w') as write_file:
            writer = csv.writer(write_file, delimiter='\t')
            writer.writerow(['doid_code', 'doid_name', 'count'])
            writer.writerows(counts)

        return merged_associations

    def read_doid_processing_info(self, path=None):
        if path is None:
            path = '/home/dhimmels/Documents/serg/gene-disease-hetnet/data-integration/doid-processing-info.txt'
        with open(path) as read_file:
            lines = list(read_file)

        remove_doids = set()
        pop_doids = dict()
        for line in lines:
            line = line.split('#', 1)[0].rstrip()
            line_list = line.split(' ')
            command = line_list.pop(0)
            if command == 'remove':
                remove_doids.add(line_list[0])
            elif command == 'pop':
                pop_doids[line_list[0]] = line_list[2]
            else:
                assert False
        return remove_doids, pop_doids



if __name__ =='__main__':
    import pprint
    gcat = GwasCatalog()
    print len(gcat.get_merged_associations('snap_wingspan'))
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