import random

import hetnet.graph

def permute_graph(graph, seed=0, verbose=False):
    """
    Shuffle edges within metaedge category. Preserves node degree but randomizes
    edges.
    """

    if verbose: print 'Creating permuted graph template'
    permuted_graph = hetnet.graph.Graph(graph.metagraph)
    for node_id, node in graph.node_dict.items():
        permuted_graph.add_node(node_id, node.metanode.id_, data=node.data)

    if verbose: print 'Retrieving graph edges'
    metaedge_to_edges = graph.get_metaedge_to_edges(exclude_inverts=True)

    if verbose: print 'Adding permuted edges'
    for metaedge, edges in metaedge_to_edges.items():
        if verbose: print metaedge
        pair_list = [(edge.source.id_, edge.target.id_) for edge in edges]
        permuted_pair_list = permute_pair_list(pair_list, seed=seed)

        kind = metaedge.kind
        direction = metaedge.direction
        for pair in permuted_pair_list:
            permuted_graph.add_edge(pair[0], pair[1], kind, direction)

    return permuted_graph


def permute_pair_list(pair_list, n_perm=None, seed=0):
    """
    If n_perm is not specific, perform 10 times the number of edges of permutations
    """
    random.seed(seed)

    pair_set = set(pair_list)
    edge_number = len(pair_list)
    if n_perm is None:
        n_perm = edge_number * 10
    for i in xrange(n_perm):

        # Same two random edges
        i_0 = random.randrange(edge_number)
        i_1 = random.randrange(edge_number)

        # Same edge selected twice
        if i_0 == i_1:
            continue
        pair_0 = pair_list.pop(i_0)
        pair_1 = pair_list.pop(i_1 - 1 if i_0 < i_1 else i_1)
        new_pair_0 = pair_0[0], pair_1[1]
        new_pair_1 = pair_1[0], pair_0[1]

        # Prevent self-connections
        self_loop = False
        duplicate = False
        for pair in new_pair_0, new_pair_1:
            if pair[0] == pair[1]:
                self_loop = True
            if pair in pair_set:
                duplicate = True

        # If edge already exists
        if duplicate or self_loop:
            for pair in pair_0, pair_1:
                pair_list.append(pair)

        # If edge is a novel permutation
        else:
            for pair in pair_0, pair_1:
                pair_set.remove(pair)
            for pair in new_pair_0, new_pair_1:
                pair_set.add(pair)
                pair_list.append(pair)

    assert len(pair_list) == edge_number
    return pair_list




