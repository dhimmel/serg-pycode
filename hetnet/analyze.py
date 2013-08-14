import os

import hetnet

class Analyzer(object):
    
    def __init__(self, graph):
        self.graph = graph

    def get_pairs(self):
        """ """
        
    def set_metapaths(self):
        """ """
    
    def select_negatives(self, positives):
        """ """
    
    def calculate_paths(self):
        """ """


class GraphAnalysis(object):
    
    def __init__(self, graph_dir):
        """ """
        if not os.path.exists(anal_dir):
            os.mkdir(anal_dir)
        self.anal_dir = anal_dir
    
    

class LearningEdges(object):
    
    def __init__(self, graph, source_metanode, target_metanode):
        self.graph = graph
    
    def set_eligible_sources(self, sources):
        """ """

    def set_eligible_targets(self, targets):
        """ """
    
    def set_positives(self):
        """ """
        
    def set_negatives(self):
        """compute from positives"""
        
    def save(self):
        """save to disk"""