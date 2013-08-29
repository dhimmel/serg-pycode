import ast
import ConfigParser
import os

import hetnet
import readwrite


class GraphAgent(object):
    
    def __init__(self, network_dir):
        """ """
        self.network_dir = network_dir
        self.path = os.path.join(network_dir, 'graph', 'graph.json.gz')
    
    def get(self):
        if not hasattr(self, 'graph'):
            self.graph = self._read()
        return self.graph
    
    def _read(self):
        """ """
        return readwrite.graph.read_json(self.path)
    
    def write(self, graph):
        """ """
        readwrite.graph.write_json(graph, self.path)


class MetaEdgeAgent(object):
    
    def __init__(self, graph_agent, identifier):
        """ """
        self.graph_agent = graph_agent
        self.identifier = identifier
        self.directory = os.path.join(self.graph_agent.network_dir, identifier)
        self.path = os.path.join(self.directory, 'metaedge.txt')

    def set(self, metaedge):
        """Set learning metaedge"""
        self.metaedge = metaedge
        self._write()
        
    def _write(self):
        """ """
        if not os.path.isdir(self.directory):
            os.mkdir(self.directory)
        readwrite.metaedge.write_text(self.metaedge, self.path)

    def get(self):
        if not hasattr(self, 'metaedge'):
            self.metaedge = self._read()
        return self.metaedge
    
    def _read(self):
        """ """
        metagraph = self.graph_agent.get().metagraph
        return readwrite.metaedge.read_text(metagraph, self.path)


class MetaPathsAgent(object):
    
    def __init__(self, metaedge_agent, identifier):
        self.metaedge_agent = metaedge_agent
        self.identifier = identifier
        directory = os.path.join(metaedge_agent.directory, 'metapaths')
        if not os.path.isdir(directory):
            os.mkdir(directory)
        self.path = os.path.join(directory, '{}.txt'.format(identifier))
    
    def set(self, metapaths):
        """ """
        self.metapaths = metapaths
        self._write()
        
    def _write(self):
        """ """
        readwrite.metapaths.write_text(self.metapaths, self.path)
    
    def get(self):
        """ """
        if not hasattr(self, 'metapaths'):
            self.metapaths = self._read()
        return self.metapaths
        
    def _read(self):
        metagraph = self.metaedge_agent.graph_agent.get().metagraph
        return readwrite.metapaths.read_text(self.path, metagraph)

class FeaturesAgent(object):
    
    def __init__(self, network_dir, identifier):
        self.network_dir = network_dir
        self.identifier = identifier
        directory = os.path.join(network_dir, 'features')
        if not os.path.isdir(directory):
            os.mkdir(directory)
        self.path = os.path.join(directory, '{}.txt'.format(identifier))

    def set(self, features):
        self.features = features
        self._write()
        
    def _write(self):
        readwrite.features.write_text(self.features, self.path)
    
    def get(self):
        if not hasattr(self, 'features'):
            self.features = self._read()
        return self.features
    
    def _read(self):
        features = readwrite.features.read_text(self.path)

    
class LearningEdgesAgent(object):

    def __init__(self, metaedge_agent, identifier):
        self.metaedge_agent = metaedge_agent
        self.identifier = identifier
        directory = os.path.join(metaedge_agent.directory, 'learning-edges')
        if not os.path.isdir(directory):
            os.mkdir(directory)
        self.path = os.path.join(directory, '{}.txt'.format(identifier))

    def set(self, learning_edges):
        self.learning_edges = learning_edges
        self._write()
    
    def _write(self):
        readwrite.learning_edges.write_text(self.learning_edges, self.path)

    def get(self):
        """ """
        if not hasattr(self, 'learning_edges'):
            self.learning_edges = self._read()
        return self.learning_edges
        
    def _read(self):
        graph = self.metaedge_agent.graph_agent.get()
        return readwrite.learning_edges.read_text(self.path, graph)







if __name__ == '__main__':
    import readwrite
    network_dir = '/home/dhimmels/Documents/serg/ashg13/130814-1'
    graph_json_path = os.path.join(network_dir, 'graph', 'graph.json.gz')
    graph = readwrite.read_json(graph_json_path)
    
    metapath_agent = MetapathAgent(graph, network_dir)
    print metapath_agent.get_learning_metaedge()
    
    
    
    