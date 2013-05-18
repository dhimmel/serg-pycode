import csv
import collections
import os
import re
import gzip

import data

class iRefIndex(object):
    
    def __init__(self, iref_dir=None):
        if iref_dir is None:
            iref_dir = data.current_path('iref')
        self.iref_dir = iref_dir
        
    def row_generator(self):
        file_names = os.listdir(self.iref_dir)
        file_name = filter(lambda s: re.match(r"9606\.mitab\.[0-9]*\.txt.gz$", s), file_names)[0]
        path = os.path.join(self.iref_dir, file_name)
        with gzip.open(path) as iref_file:
            fieldnames = iref_file.next().rstrip().lstrip('#').split('\t')
        iref_file = gzip.open(path)
        reader = csv.DictReader(iref_file, delimiter='\t', fieldnames=fieldnames)
        reader.next() # skip header
        for row in reader:
            for fieldname, field in row.items():
                row[fieldname] = field.split('|')
            
            for field in 'altA', 'altB':
                row[field] = {alt_id.split(':')[0]: alt_id.split(':')[1] for alt_id in row[field]}
                
            yield row
        iref_file.close()
        
    def row_to_interaction(self, row):
        entrez_a = row['altA'].get('entrezgene/locuslink')
        entrez_b = row['altB'].get('entrezgene/locuslink')
        if not entrez_a or not entrez_b:
            return None
        if entrez_a == entrez_b:
            # Exclude self-interaction
            return None
        return frozenset([entrez_a, entrez_b])
    
    def all_interactions(self, symbol=True):
        entrez_to_gene = data.Data().hgnc.get_entrez_to_gene()
        hgnc_entrez_ids = set(entrez_to_gene)
        rows = self.row_generator()
        all_interactions = {self.row_to_interaction(row) for row in rows}
        all_interactions.remove(None)
        all_interactions = filter(lambda fset: fset <= hgnc_entrez_ids, all_interactions)
        all_interactions = map(tuple, all_interactions)
        if symbol:
            all_interactions = [(entrez_to_gene[interaction[0]].symbol,
                                 entrez_to_gene[interaction[1]].symbol)
                                 for interaction in all_interactions]
        return all_interactions                
                
        
    
if __name__ =='__main__':
    iref = iRefIndex()
    rows = iref.row_generator()
    for i, row in enumerate(rows):
        if i % 1000 == 0:
            print i, row['altA'].get('entrezgene/locuslink'), row['altA'].get('entrezgene/locuslink')
    #print len(iref.all_interactions())
    

    
    