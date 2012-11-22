import collections

import networkx

import bioparser.data

pkl_path = '/home/dhimmels/Documents/serg/ipanet/ipanet.pkl'
g = networkx.read_gpickle(pkl_path)
#print networkx.info(g)
#node_kind_counts = collections.Counter(data['kind'] for node, data in g.nodes_iter(data=True))
#print node_kind_counts


#print 'Number of connected components:', networkx.number_connected_components(g)


"""
for node, data in g.nodes_iter(data=True):
    print node
    print data
    print g.neighbors(node)
"""    







################################################################################
################################# Network Stats ################################
#print 'calculating statistics'


if __name__ == '__main__':


    # pagerank
    #pr = networkx.pagerank(g)
    #pr_sorted = sorted(pr.items(), key=lambda x: x[1])
    #print pr_sorted[:5]
    pass