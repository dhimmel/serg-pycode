import os

import scipy.stats
import statsmodels.distributions.empirical_distribution
import matplotlib.pyplot as plt
import networkx
import numpy

import nxutils

def node_degree(g, node, edge_kind):
    """ """
    edges = g.edges(node, keys=True)
    edges = filter(lambda e: e[2] == edge_kind, edges)
    return len(edges)

def node_degrees(g, node_kind, edge_kind):
    """Return a dictionary with nodes of the specified kind as keys and the
    number of incident edges of the specified kind as values.
    """
    kind_to_nodes = nxutils.get_kind_to_nodes(g)
    #kind_to_edges = nxutils.get_kind_to_edges(g)
    nodes = kind_to_nodes[node_kind]
    node_to_degree = {node: node_degree(g, node, edge_kind) for node in nodes}
    return node_to_degree

def node_degrees_from_edges(g, edge_kind, edges):
    """ """
    node_kind_to_degrees = dict()
    for edge in edges:
        for node, edge_kind in [(edge[0], edge_kind), (edge[1], (edge_kind[1], edge_kind[0], edge_kind[2]))]:
            node_kind = g.node[node]['kind']
            edge_kind_counter = nxutils.node_degree_counter(g, node)
            degree = edge_kind_counter[edge_kind]
            node_kind_to_degrees.setdefault(node_kind, list()).append(degree)
    return node_kind_to_degrees

def plot_ecdf(values, **kwargs):
    ecdf = statsmodels.distributions.empirical_distribution.ECDF(values)
    x = numpy.arange(0, max(values))
    y = ecdf(x)
    plt.step(x, y, **kwargs)

def cdf_degree_plot(g):
    """For g.graph['edge_kind'] edges."""
    edge_kind = g.graph['edge_kind']
    kind_to_edges = networkx_extensions.get_kind_to_edges(g)
    remaining_edges = kind_to_edges[edge_kind]
    negatives = g.graph['negatives']
    positives = g.graph['positives']
    
    
    
    for edges, subset_label in ([positives, 'positives'],
                                [negatives, 'negatives'],
                                [remaining_edges, 'remaining']):
        node_kind_to_degrees = node_degrees_from_edges(g, edge_kind, edges)
        for node_kind, degrees in node_kind_to_degrees.items():
            label = subset_label + ': ' + node_kind
            plot_ecdf(degrees, label=label)

    
    source_kind = g.graph['source_kind']
    target_kind = g.graph['target_kind']
    source_degrees = node_degrees(g, source_kind, edge_kind).values()
    target_degrees = node_degrees(g, target_kind, edge_kind).values()
    plot_ecdf(source_degrees, label=source_kind, linestyle='--')
    plot_ecdf(target_degrees, label=target_kind, linestyle='--')
    
    plt.legend(loc=4, title='Edges: Nodes')   

    plt.title('Node Degree Distriubtion for ' + edge_kind + ' Edges')
    plt.xlabel('Node Degree')
    plt.ylabel('Cumulative Density')
    

def cdf_degree_plot2(g):
    edge_kind = g.graph['source_kind'], g.graph['target_kind'], g.graph['edge_key']
    kind_to_edges = nxutils.get_kind_to_edges(g)
    negatives = g.graph['negatives']
    positives = g.graph['positives']
    
    
    
    for edges, subset_label in ([positives, 'positives'],
                                [negatives, 'negatives']):
        node_kind_to_degrees = node_degrees_from_edges(g, edge_kind, edges)
        for node_kind, degrees in node_kind_to_degrees.items():
            label = subset_label + ': ' + node_kind
            print label, degrees
            plot_ecdf(degrees, label=label)
    
    plt.legend(loc=4, title='Edges: Nodes')   

    plt.title('Node Degree Distriubtion for ' + edge_kind[2] + ' Edges')
    plt.xlabel('Node Degree')
    plt.ylabel('Cumulative Density')
    

def cdf_degree_plot3(g):
    source_kind = g.graph['source_kind']
    target_kind = g.graph['target_kind']
    #edge_key = g.graph['edge_key']
    negatives = g.graph['negatives']
    positives = g.graph['positives']
    edge_kind = g.graph['source_kind'], g.graph['target_kind'], g.graph['edge_key']
    edge_kind_rev = g.graph['target_kind'], g.graph['source_kind'], g.graph['edge_key']
    #nxutils.node_degree_counter(g, node)
    source_degrees = list()
    target_degrees = list()
    for edges, subset_label in ([positives, 'positives'],
                                [negatives, 'negatives']):
        for node, neighbor, edge_key in edges:
            source_degrees.append(nxutils.node_degree_counter(g, node)[edge_kind])
            target_degrees.append(nxutils.node_degree_counter(g, neighbor)[edge_kind_rev])
        label = subset_label + ': ' + source_kind
        plot_ecdf(source_degrees, label=label)
        label = subset_label + ': ' + target_kind
        plot_ecdf(target_degrees, label=label)
    
    plt.legend(loc=4, title='Edges: Nodes')   

    plt.title('Node Degree Distriubtion for ' + edge_kind[2] + ' Edges')
    plt.xlabel('Node Degree')
    plt.ylabel('Cumulative Density')
    
                




#node_degrees(g, 'drug', 'indication').itervalues()
#node_degrees(g, 'disease', 'indication').itervalues()


if __name__ == '__main__':
    """
    ipanet_dir = '/home/dhimmels/Documents/serg/ipanet/'
    network_id = '130116-1'
    pkl_path_prepared = os.path.join(ipanet_dir, 'networks', network_id, 'prepared-graph.pkl')
    pkl_path = os.path.join(ipanet_dir, 'networks', network_id, 'graph.pkl')
    g = networkx.read_gpickle(pkl_path_prepared)
    print 'read gpickle'
    
    cdf_degree_plot(g)
    #figure = plt.figure(figsize=(5, 4))
    fig_path = os.path.join(ipanet_dir, 'networks', network_id, 'cumulative-degree-plot.eps')
    plt.savefig(fig_path, format="eps")
    #plt.show()
    """
    path = '/home/dhimmels/Documents/serg/ipanet/networks/130222-1/graphs/learning-graph.pkl'
    g = networkx.read_gpickle(path)
    print 'read network pickle'
    cdf_degree_plot3(g)
    fig_path = os.path.join('/home/dhimmels/Dropbox/serg/lab-meeting-130318', 'degree_plot_130222-1.eps')
    plt.savefig(fig_path, format="eps")
    plt.close()
    
    
    path = '/home/dhimmels/Documents/serg/ipanet/networks/130225-1/graphs/learning-graph.pkl'
    g = networkx.read_gpickle(path)
    print 'read network pickle'
    print len(g.graph['positives']), len(g.graph['negatives'])
    cdf_degree_plot3(g)
    fig_path = os.path.join('/home/dhimmels/Dropbox/serg/lab-meeting-130318', 'degree_plot_130225-1.eps')
    plt.savefig(fig_path, format="eps")
    plt.close()

    path = '/home/dhimmels/Documents/serg/networks/130301-1/graphs/learning-graph.pkl'
    g = networkx.read_gpickle(path)
    print 'read network pickle'
    cdf_degree_plot3(g)
    fig_path = os.path.join('/home/dhimmels/Dropbox/serg/lab-meeting-130318', 'degree_plot_130301-1.eps')
    plt.savefig(fig_path, format="eps")
    plt.close()



