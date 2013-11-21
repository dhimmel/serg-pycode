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


if __name__ =='__main__':
    do = DO()
    onto = do.get_ontology()

    onto.annotate_information_content()
    g = onto.graph
    
    """
    import itertools
    for n0, n1 in itertools.combinations(g.nodes(), 2):
        onto.intrinsic_similarity(n0, n1)
    """
    doids = ['DOID:1686', 'DOID:11830', 'DOID:9008', 'DOID:13189', 'DOID:585', 'DOID:0050156', 'DOID:5379', 'DOID:678', 'DOID:5408', 'DOID:1909', 'DOID:12236', 'DOID:3083', 'DOID:3324', 'DOID:6364', 'DOID:824', 'DOID:11949', 'DOID:1067', 'DOID:9074', 'DOID:2914', 'DOID:9970', 'DOID:10608', 'DOID:8398', 'DOID:7147', 'DOID:11476', 'DOID:2043', 'DOID:289', 'DOID:13099', 'DOID:7148', 'DOID:12306', 'DOID:10763', 'DOID:684', 'DOID:8552', 'DOID:2841', 'DOID:769', 'DOID:2355', 'DOID:3963', 'DOID:3620', 'DOID:3459', 'DOID:2986', 'DOID:14330', 'DOID:3455', 'DOID:11829', 'DOID:1107', 'DOID:8893', 'DOID:11119', 'DOID:9744', 'DOID:594', 'DOID:1040', 'DOID:3905', 'DOID:3908', 'DOID:4450', 'DOID:216', 'DOID:4905', 'DOID:9352', 'DOID:4481', 'DOID:1024', 'DOID:6536', 'DOID:11714', 'DOID:1380', 'DOID:9296', 'DOID:14227', 'DOID:14221', 'DOID:5517', 'DOID:2468', 'DOID:13378', 'DOID:1936', 'DOID:3393', 'DOID:11555', 'DOID:12995', 'DOID:557', 'DOID:5844', 'DOID:12361', 'DOID:12365', 'DOID:13241', 'DOID:2377', 'DOID:4960', 'DOID:12849', 'DOID:635', 'DOID:3369', 'DOID:1094', 'DOID:8778', 'DOID:11782', 'DOID:8567', 'DOID:1826', 'DOID:10652', 'DOID:0050425', 'DOID:10286', 'DOID:332', 'DOID:13641', 'DOID:4001', 'DOID:12930', 'DOID:4007', 'DOID:10976', 'DOID:1485', 'DOID:10487', 'DOID:9835', 'DOID:12185', 'DOID:11335', 'DOID:0050742', 'DOID:3312', 'DOID:9206', 'DOID:1312', 'DOID:0050589', 'DOID:13608', 'DOID:8577', 'DOID:10941', 'DOID:13223', 'DOID:7693', 'DOID:2154', 'DOID:10871', 'DOID:14268', 'DOID:9538', 'DOID:399', 'DOID:3910', 'DOID:986', 'DOID:1459', 'DOID:10892', 'DOID:11612', 'DOID:0050741', 'DOID:1595', 'DOID:706', 'DOID:5419', 'DOID:9952', 'DOID:90', 'DOID:8986', 'DOID:784', 'DOID:3310', 'DOID:418']

    
    import itertools
    for n0, n1 in itertools.combinations(doids, 2):
        print '\t'.join([g.node[n0]['name'], g.node[n1]['name'],
              str(round(onto.intrinsic_similarity(n0, n1)['resnik_similarity'], 3)),
              str(round(onto.intrinsic_similarity(n0, n1)['lin_similarity'], 3)),
              str(round(onto.intrinsic_similarity(n0, n1)['jen_similarity'], 3))
              ])
