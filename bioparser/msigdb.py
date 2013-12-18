import os
import csv

import data

# http://www.broadinstitute.org/gsea/downloads.jsp
abbrev_to_name = {
    'c1.all': 'positional gene sets',
    'c2.all': 'all curated gene sets',
    'c2.cgp': 'chemical and genetic perturbations',
    'c2.cp': 'all canonical pathways',
    'c2.cp.biocarta': 'BioCarta gene sets',
    'c2.cp.kegg': 'KEGG gene sets',
    'c2.cp.reactome': 'Reactome gene sets',
    'c3.all': 'all motif gene sets',
    'c3.mir': 'microRNA targets',
    'c3.tft': 'transcription factor targets',
    'c4.all': 'all computational gene sets',
    'c4.cgn': ' cancer gene neighborhoods',
    'c4.cm': 'cancer modules',
    'c5.all': 'all GO gene sets',
    'c5.bp': 'GO biological processes',
    'c5.cc': 'GO cellular components',
    'c5.mf': 'GO molecular functions',
    'c6.all': 'all oncogenic signatures gene sets',
    'c7.all': 'all immunologic signatures gene sets'
}


class MSigDB(object):

    def __init__(self, directory=None):
        if directory is None:
            directory = data.current_path('msigdb')
        self.directory = directory
        self.abbrev_to_name = abbrev_to_name

    @staticmethod
    def gmt_row_generator(path):
        """
        (name, description, members)
        """
        read_file = open(path)
        reader = csv.reader(read_file, delimiter='\t')
        for row in reader:
            yield row[0], row[1], row[2:]
        read_file.close()

    def get_gmt_path(self, abbrev, version='v4.0', id_type='entrez'):
        filename = '{}.{}.{}.gmt'.format(abbrev, version, id_type)
        path = os.path.join(self.directory, filename)
        return path

    def gene_set_generator(self, abbrev):
        identifiers_to_genes = data.Data().hgnc.identifiers_to_genes
        path = self.get_gmt_path(abbrev, id_type='entrez')
        for name, description, members in self.gmt_row_generator(path):
            genes = identifiers_to_genes(members, id_type='entrez')
            yield name, description, genes



if __name__ == '__main__':
    msigdb = MSigDB()
    for abbrev in abbrev_to_name.keys():
        list(msigdb.gene_set_generator(abbrev))
