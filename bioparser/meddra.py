import os
import csv

import networkx

import data

class MedDRA(object):
    
    def __init__(self, directory=None):
        if not directory:
            directory = data.current_path('meddra')
        self.directory = directory
    
    def get_networkx(self):
        """
        Creates a networx DiGraph to represent the MedDRA ontology. Each
        node is identified using it's meddra code. Additional attributes are
        stored for each node. Lower level terms are excluded since one LLT
        shares the same code as a preferred term.
        """
        graph = networkx.DiGraph(source=self.directory)
        self.graph = graph
        
        # Add nodes excluding llt
        levels = ['soc', 'hlgt', 'hlt', 'pt']
        for level in levels:
            read_function = getattr(self, 'read_{}'.format(level))
            for row in read_function():
                code = row['{}_code'.format(level)]
                name = row['{}_name'.format(level)]
                row['name'] = name
                row['level'] = level
                graph.add_node(code, row)
        
        # Add edges
        relation_files = ['soc_hlgt', 'hlgt_hlt', 'hlt_pt']
        fieldnames = ['parent_code', 'child_code']
        for name in relation_files:
            for relation in self.read_asc(name, fieldnames):
                graph.add_edge(relation['parent_code'], relation['child_code'])
        
        assert networkx.is_directed_acyclic_graph(graph)
        return self.graph        
    
    def read_asc(self, name, fieldnames):
        """
        Returns a generator of dictionaries which each a row of the file
        specified using the name argument. Do not include the .asc extension in
        name. Fieldnames are the labels for the columns. For fieldname
        documentation refer to dist_file_format_xx_x_English.pdf.
        """
        path = os.path.join(self.directory, 'MedAscii', '{}.asc'.format(name))
        read_file = open(path)
        reader = csv.DictReader(read_file, fieldnames=fieldnames, delimiter='$')
        for row in reader:
            yield row
        read_file.close()

    def read_soc(self):
        """System Organ Class"""
        fieldnames = ['soc_code', 'soc_name', 'soc_abbrev',
            'soc_whoart_code', 'soc_harts_code', 'soc_costart_sym', 
            'soc_icd9_code', 'soc_icd9cm_code', 'soc_icd10_code', 'soc_jart_code']
        return self.read_asc('soc', fieldnames)

    def read_hlgt(self):
        """High Level Group Term"""
        fieldnames = ['hlgt_code', 'hlgt_name', 'hlgt_whoart_code',
            'hlgt_harts_code', 'hlgt_costart_sym', 'hlgt_icd9_code', 
            'hlgt_icd9cm_code', 'hlgt_icd10_code', 'hlgt_jart_code']
        return self.read_asc('hlgt', fieldnames)

    def read_hlt(self):
        """High Level Term"""
        fieldnames = ['hlt_code', 'hlt_name', 'hlt_whoart_code',
            'hlt_harts_code', 'hlt_costart_sym', 'hlt_icd9_code', 
            'hlt_icd9cm_code', 'hlt_icd10_code', 'hlt_jart_code']
        return self.read_asc('hlt', fieldnames)

    def read_pt(self):
        """Preferred Term"""
        fieldnames = ['pt_code', 'pt_name', 'null_field', 'pt_soc_code',
            'pt_whoart_code', 'pt_harts_code', 'pt_costart_sym', 'pt_icd9_code',
            'pt_icd9cm_code', 'pt_icd10_code', 'pt_jart_code']
        return self.read_asc('pt', fieldnames)

    def read_llt(self):
        """Lowest Level Term"""
        fieldnames = ['llt_code', 'llt_name', 'pt_code', 'llt_whoart_code',
            'llt_harts_code', 'llt_costart_sym', 'llt_icd9_code', 'llt_icd9cm_code',
            'llt_icd10_code', 'llt_currency', 'llt_jart_code']
        return self.read_asc('llt', fieldnames)


    
if __name__ == '__main__':
    meddra = MedDRA()
    graph = meddra.get_networkx()
    print graph.predecessors('10002043')
