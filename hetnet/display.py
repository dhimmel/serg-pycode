

def graph_metrics(graph):
    lines = list()

    metanode_to_nodes = graph.get_metanode_to_nodes()
    total_nodes = len(graph.node_dict)
    connected_nodes = sum(any(node.edges.values()) for node in graph.node_dict.values())
    lines.append('total nodes: {}, connected nodes: {}'.format(total_nodes, connected_nodes))
    lines.append('metanode counts:')
    for metanode, nodes in metanode_to_nodes.items():
        line = '{}: {}'.format(metanode, len(nodes))
        lines.append(line)
    lines.append('')

    metaedge_to_edges = graph.get_metaedge_to_edges(exclude_inverts=True)
    lines.append('edge count: {}'.format(sum(map(len, metaedge_to_edges.values()))))
    lines.append('metaedge counts:')
    for metaedge, edges in metaedge_to_edges.items():
        line = '{}: {}'.format(metaedge, len(edges))
        lines.append(line)
    lines.append('')

    return '\n'.join(lines)