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
        obo_filename = 'HumanDO.obo'
        keep_attributes = ['name', 'xref', 'def', 'synonym']
        super(DO, self).__init__(directory, obo_filename, keep_attributes)

    def get_graph(self):
        """
        """
        graph = super(DO, self).get_graph()
        for node, node_data in graph.nodes_iter(data=True):
            node_data['int_id'] = int(node.split(':')[1])
            # takes the first definition although oftentimes multiple exist
            """
            xref_dict = dict()
            for xref in node_data.get('xref', list()):
                key, value = xref.split(':')
                xref_dict.setdefault(key, list()).append(xref)
            node_data['xref_dict'] = xref_dict
            """
        return graph

    def annotate_categories(self):
        """
        annotate nodes attributes to include category_codes and category_names.
        """
        graph = self.get_graph()
        category_codes = set()
        category_codes |= set(graph.successors('DOID:7')) # disease of anatomical entity
        category_codes.add('DOID:14566') # disease of cellular proliferation
        category_codes.add('DOID:150') # disease of mental health
        category_codes.add('DOID:0014667') # disease of metabolism
        for node, node_data in graph.nodes_iter(data=True):
            node_category_codes = category_codes & networkx.ancestors(graph, node)
            node_category_names = {graph.node[code]['name'] for code in node_category_codes}
            node_data['category_codes'] = node_category_codes
            node_data['category_names'] = node_category_names

    def get_doid_to_xrefs(self, code, prepend=''):
        """
        code is the vocabulary code to extract
        """
        graph = self.get_graph()
        doid_to_xrefs = dict()
        for node, node_data in graph.nodes_iter(data=True):
            xref_dict = node_data.get('xref', dict())
            xrefs = xref_dict.get(code, list())
            xrefs = {'{}{}'.format(prepend, xref) for xref in xrefs}
            doid_to_xrefs[node] = xrefs
        return doid_to_xrefs

    def get_xref_to_doids(self, code, prepend=''):
        doid_to_xrefs = self.get_doid_to_xrefs(code, prepend)
        xref_to_doids = dict()
        for doid, xrefs in doid_to_xrefs.iteritems():
            for xref in xrefs:
                xref_to_doids.setdefault(xref, set()).add(doid)
        return xref_to_doids


if __name__ =='__main__':
    do = DO()
    #onto = do.get_ontology()

    #onto.annotate_information_content()
    #g = onto.graph
    import pprint
    pprint.pprint(do.get_doid_to_xrefs('OMIM', 'OMIM:'))


