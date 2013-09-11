import ast
import ConfigParser
import os

import hetnet
import readwrite


class GraphAgent(object):
    
    def __init__(self, network_dir):
        """ """
        self.network_dir = network_dir
        self.graph_dir = os.path.join(network_dir, 'graph')
        for directory in (self.network_dir, self.graph_dir):
            if not os.path.isdir(directory):
                os.mkdir(directory)
        self.path = os.path.join(self.graph_dir, 'graph.pkl.gz')
    
    def get(self):
        if not hasattr(self, 'graph'):
            self.graph = self._read()
        return self.graph
    
    def _read(self):
        """ """
        return readwrite.graph.read_pickle(self.path)
    
    def set(self, graph):
        self.graph = graph
        self._write()
    
    def _write(self):
        """ """
        readwrite.graph.write_pickle(self.graph, self.path)
        #print datetime.datetime.now().time().isoformat()

    def write_additional_formats(self):
        yaml_path = os.path.join(self.graph_dir, 'graph.yaml.gz')
        readwrite.graph.write_yaml(self.graph, yaml_path)
        json_path = os.path.join(self.graph_dir, 'graph.json.gz')
        readwrite.graph.write_json(self.graph, json_path)


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
    
    def __init__(self, graph_agent, identifier):
        self.graph_agent = graph_agent
        self.identifier = identifier
        directory = os.path.join(graph_agent.network_dir, 'features')
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
        return readwrite.features.read_text(self.path)

    
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

class AnalysisAgent(object):
    
    def __init__(self, metaedge_agent, identifier):
        """ """
        self.metaedge_agent = metaedge_agent
        self.graph_agent = metaedge_agent.graph_agent
        #self.network_dir = self.graph_agent.network_dir
        self.identifier = identifier
        
        # analyses directory
        directory = os.path.join(metaedge_agent.directory, 'analyses')
        if not os.path.isdir(directory):
            os.mkdir(directory)
        
        # analyses/identifier directory
        directory = os.path.join(directory, identifier)
        if not os.path.isdir(directory):
            os.mkdir(directory)
        self.directory = directory
        self.path = os.path.join(directory, 'arguments.cfg')
    
    def set(self, features_id, metapaths_id, learning_edges_id):
        self.features_agent = hetnet.agents.FeaturesAgent(self.graph_agent, features_id)
        self.metapaths_agent = hetnet.agents.MetaPathsAgent(self.metaedge_agent, metapaths_id)
        self.learning_edges_agent = hetnet.agents.LearningEdgesAgent(self.metaedge_agent, learning_edges_id)
        self._write()
    
    def write_config(self, features_id='', metapaths_id='', learning_edges_id=''):
        config = ConfigParser.SafeConfigParser()

        section = 'redundant'
        config.add_section(section)
        config.set(section, 'network_dir', self.graph_agent.network_dir)
        config.set(section, 'metaedge_id', self.metaedge_agent.identifier)

        section = 'analysis'
        config.add_section(section)
        config.set(section, 'features_id', features_id)
        config.set(section, 'metapaths_id', metapaths_id)
        config.set(section, 'learning_edges_id', learning_edges_id)
        
        with open(self.path, 'w') as write_file:
            config.write(write_file)
    
    def _write(self):
        features_id = self.features_agent.identifier
        metapaths_id = self.metapaths_agent.identifier
        learning_edges_id = self.learning_edges_agent.identifier
        self.write_config(features_id, metapaths_id, learning_edges_id)
    
    def get(self):
        self._read()
        return #TODO
    
    def _read(self):
        config = ConfigParser.SafeConfigParser()
        config.read(self.path)
        section = 'analysis'
        features_id = config.get(section, 'features_id')
        self.features_agent = FeaturesAgent(self.graph_agent, features_id)
        metapaths_id = config.get(section, 'metapaths_id')
        self.metapaths_agent = MetaPathsAgent(self.metaedge_agent, metapaths_id)
        learning_edges_id = config.get(section, 'learning_edges_id')
        self.learning_edges_agent = LearningEdgesAgent(self.metaedge_agent, learning_edges_id)
    

if __name__ == '__main__':
    import readwrite
    network_dir = '/home/dhimmels/Documents/serg/ashg13/130814-1'
    graph_json_path = os.path.join(network_dir, 'graph', 'graph.json.gz')
    graph = readwrite.read_json(graph_json_path)
    
    metapath_agent = MetapathAgent(graph, network_dir)
    print metapath_agent.get_learning_metaedge()
    
    
    
    