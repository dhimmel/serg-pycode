import csv
import os
import re
import itertools
import collections

import utilities.omictools
import bioparser.data

    
class MorbidMap(object):

    def __init__(self, omim_dir=None):
        if not omim_dir:
            omim_dir = bioparser.data.current_path('omim')
        self.omim_dir = omim_dir

    def read_morbidmap(self):
        path = os.path.join(self.omim_dir, 'morbidmap')
        read_file = open(path)
        fieldnames = 'disorder_info', 'symbols', 'omim_gene_id', 'locus'
        reader = csv.DictReader(read_file, delimiter='|', fieldnames=fieldnames)
        for row in reader:
            row['symbols'] = row['symbols'].split(', ')
            yield row
        read_file.close()
    
    def get_associations(self, exclude_mimless_disorders=True):
        """http://www.omim.org/help/faq"""
        disorder_name_pattern = r'^[{[?]*(.*?)([}\]]|, [0-9]{6}| \([1-4]\)$)'
        mim_pattern = re.compile(r"[0-9]{6}")
        association_pattern = re.compile(r"\(([1-4])\)$")
        hgnc = bioparser.data.Data().hgnc
        symbol_to_gene = hgnc.get_symbol_to_gene()
        id_gene_tuples = set()
        associations = list()
        for row in self.read_morbidmap():
            disorder_info = row['disorder_info']
            disorder_name = re.search(disorder_name_pattern, disorder_info).group(1)
            mim_match = re.search(mim_pattern, disorder_info)
            mim_number = mim_match.group(0) if mim_match else None
            if exclude_mimless_disorders and not mim_number:
                continue
            if '[' in disorder_info:
                disorder_type = 'nondisease'
            elif '{' in disorder_info:
                disorder_type = 'multifactorial'
            else:
                disorder_type = 'mendelian'

            confirmed = 0 if '?' in disorder_info else 1

            association_type = re.search(association_pattern, disorder_info).group(1)
            
            genes = set(symbol_to_gene.get(symbol) for symbol in row['symbols'])
            genes.discard(None)
            for gene in genes:
                association= {'mim_number': mim_number, 'gene': gene,
                              'disorder_name': disorder_name,
                              'disorder_type': disorder_type,
                              'association_type': association_type,
                              'confirmed': confirmed}
                associations.append(association)
        id_gene_tuples = list(id_gene_tuples)
        return associations

    def get_disorders(self):
        """ """
        associations = self.get_associations(exclude_mimless_disorders=True)
        disorder_attributes = ('mim_number', 'disorder_name', 'disorder_type')
        disorder_tuples = set()
        for association in associations:
            disorder_tuple = tuple(association[key] for key in disorder_attributes)
            disorder_tuples.add(disorder_tuple)
        
        # For mim disorder IDs with multiple names take the first name alphabetically
        disorder_tuples = sorted(disorder_tuples)
        mim_to_disorder = dict()
        for disorder_tuple in disorder_tuples:
            disorder = dict(zip(disorder_attributes, disorder_tuple))
            mim_number = disorder['mim_number']
            if mim_number not in mim_to_disorder:
                mim_to_disorder[mim_number] = disorder
        disorders = mim_to_disorder.values()
        return disorders


if __name__ == '__main__':
    
    mm = MorbidMap()
    disorders = mm.get_disorders()
    import pprint
    pprint.pprint(disorders)
        
    