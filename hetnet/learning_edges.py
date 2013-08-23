import ast
import ConfigParser
import os

import hetnet

        
class MetapathAgent(object):
    
    def __init__(self, graph, network_dir):
        self.graph = graph
        self.metapaths_dir = os.path.join(network_dir, 'metapaths')
    
    def write_learning_metaedge_config(self, metaedge):
        """ """
        metaedge_id = metaedge.get_id()
        config.add_section(section)
        config_path = os.path.join(self.metapaths_dir, 'learning-metaedge.cfg')
        config = ConfigParser.SafeConfigParser()
        config.add_section('learning metaedge')
        config.set('learning metaedge', 'edge_id', metaedge_id)
        with open(config_path, 'w') as config_file:
            config.write(config_path)
    
    def get_learning_metaedge(self):
        """ """
        if not hasattr(self, 'learning_metaedge'):
            config_path = os.path.join(self.metapaths_dir, 'learning-metaedge.cfg')
            config = ConfigParser.SafeConfigParser()
            with open(config_path) as config_file:
                config.read(config_path)
            edge_id_str = config.get('learning metaedge', 'edge_id')
            learning_edge_id = ast.literal_eval(edge_id_str)
            self.learning_metaedge = graph.metagraph.edge_dict[learning_edge_id]
        return self.learning_metaedge
        

class LearningEdgeAgent(object):
    
    def __init__(self, graph, network_dir):
        self.graph = graph
        os.path.join(network_dir, 'metapaths')
        
        self.positives = dict()
        self.negatives = dict()
        
    def set_positives(self, edges=None):
        """
        edges=None specifies all edges of the specified type should be positive
        """
        
        if edges is None:
            self.graph


if __name__ == '__main__':
    import readwrite
    network_dir = '/home/dhimmels/Documents/serg/ashg13/130814-1'
    graph_json_path = os.path.join(network_dir, 'graph', 'graph.json.gz')
    graph = readwrite.read_json(graph_json_path)
    
    metapath_agent = MetapathAgent(graph, network_dir)
    print metapath_agent.get_learning_metaedge()
    
    
    
    