# Parse Human Disease Ontology

class DO(object):
    
    def __init__(self, doid_dir=None):
        if doid_dir is None:
            doid_dir = data.current_path('doid')
        self.doid_dir = doid_dir
        self.obo_path = os.path.join(doid_dir, 'HumanDO.obo.txt')
        self.graph_json_path = os.path.join(doid_dir, 'efo.networkx.json')
        self.graph_pkl_path = os.path.join(doid_dir, 'efo.networkx.pkl')
