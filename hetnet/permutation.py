import random

import hetnet.graph

def permute_graph(graph, multiplier=10, seed=0, metaedge_to_excluded=dict(), verbose=False):
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

        excluded_pair_set = metaedge_to_excluded.get(metaedge, set())
        pair_list = [(edge.source.id_, edge.target.id_) for edge in edges]
        directed = metaedge.direction != 'both'
        permuted_pair_list = permute_pair_list(
            pair_list, directed=directed, multiplier=multiplier,
            excluded_pair_set=excluded_pair_set, seed=seed, verbose=verbose)

        kind = metaedge.kind
        direction = metaedge.direction
        for pair in permuted_pair_list:
            permuted_graph.add_edge(pair[0], pair[1], kind, direction)

    return permuted_graph


def permute_pair_list(pair_list, directed=False, multiplier=10, excluded_pair_set=set(), seed=0, verbose=False):
    """
    If n_perm is not specific, perform 10 times the number of edges of permutations
    May not work for directed edges
    """
    random.seed(seed)

    pair_set = set(pair_list)
    assert len(pair_set) == len(pair_list)

    edge_number = len(pair_list)
    n_perm = int(edge_number * multiplier)

    if verbose:
        orig_pair_set = pair_set.copy()
        print '{} edges, {} permutations (seed = {}, directed = {}, {} excluded_edges)'.format(
            edge_number, n_perm, seed, directed, len(excluded_pair_set))
        print_at = range(0, n_perm, n_perm / 10) + [n_perm - 1]

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
        excluded = False
        for pair in new_pair_0, new_pair_1:
            if pair[0] == pair[1]:
                self_loop = True
            if pair in pair_set:
                duplicate = True
            if not directed and (pair[1], pair[0]) in pair_set:
                duplicate = True
            if pair in excluded_pair_set:
                excluded = True

        # If edge is invalid
        if duplicate or self_loop or excluded:
            for pair in pair_0, pair_1:
                pair_list.append(pair)

        # If edge is a novel permutation
        else:
            for pair in pair_0, pair_1:
                pair_set.remove(pair)
            for pair in new_pair_0, new_pair_1:
                pair_set.add(pair)
                pair_list.append(pair)

        # print updates
        if verbose and i in print_at:
            percent_done = 100.0 * float(i) / n_perm
            percent_same = 100.0 * float(len(orig_pair_set & pair_set)) / len(pair_set)
            print '{:.1f}% complete: {:.1f}% unchanged'.format(percent_done, percent_same)

    assert len(pair_set) == edge_number
    return pair_list




