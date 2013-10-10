# Parse Human Disease Ontology
import os
import re

import networkx

import networkx_ontology
import obo
import data

class DO(obo.OBO):
    
    def __init__(self, directory=None):
        if directory is None:
            directory = data.current_path('doid')
        obo_filename = 'doid.obo'
        keep_attributes = ['name', 'xref', 'def', 'synonym']
        super(DO, self).__init__(directory, obo_filename, keep_attributes)

    def get_graph(self):
        """
        """
        graph = super(DO, self).get_graph()
        pattern = re.compile(r'^["](.*?)["]')
        for node, data in graph.nodes_iter(data=True):
            data['id_'] = int(node.split(':')[1])
            # takes the first definition although oftentimes multiple exist
            definition = data['def_'] or ''
            match = re.search(pattern, definition)
            data['def_'] = match.group(1) if match else None
            
            xref_dict = dict()
            #if isinstance(data['xref'], str):
            #    data['xref'] = [data['xref']]
            for xref in data['xref'] or list():
                key, value = xref.split(':')
                xref_dict.setdefault(key, list()).append(value)
            data['xref'] = xref_dict
        return graph


if __name__ =='__main__':
    do = DO()
    do.get_graph()
