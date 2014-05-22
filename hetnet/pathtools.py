import resource
import sys
import collections
import itertools
import random

import hetnet


cache_gets = 0
cache_sets = 0
max_MB = 50 * 1024
prune_fraction = 0.2
memcheck_frequency = 0.01
cache = collections.OrderedDict()


def memory_usage_resource():
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom *= rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem

def cache_get(key):
    global cache_gets
    cache_gets += 1
    value = cache.pop(key)
    cache[key] = value
    return value

def cache_set(key, value):
    global cache_sets
    cache_sets += 1
    cache[key] = value
    if random.random() <= memcheck_frequency and memory_usage_resource() > max_MB:
        n_remove = int(len(cache) * prune_fraction)
        remove_keys = tuple(itertools.tee(cache.iterkeys(), n_remove))
        for key in remove_keys:
            del cache[key]
        memory = memory_usage_resource()
        hitrate = cache_hit_rate()
        print('Deleted {} cached items. Using {:.1f}GB for {:.2f} item cache. Hitrate {:.5f}'.format(
            n_remove, memory / 1024.0, hitrate))

def cache_hit_rate():
    return float(cache_gets) / cache_sets


def crdfs_paths_from(node, metapath):
    """
    Cached recursive depth-first-search: computes all paths from
    source_node of kind metapath. Paths with duplicate nodes are excluded.
    Returns a tuple of tuple paths where the elements of the tuple path are
    hetnet.Edge() objects. Refer to the cache_get and cache_set functions for
    the specifics of the caching algorithm.
    """
    if not metapath:
        return tuple(),
    args = node, metapath
    if args in cache:
        return cache_get(args)
    paths = list()
    metapath_tail = metapath.sub
    for edge in node.edges[metapath[0]]:
        for tail in crdfs_paths_from(edge.target, metapath_tail):
            if node in (e.target for e in tail):
                continue
            paths.append((edge, ) + tail)
    paths = tuple(paths)
    cache_set(args, paths)
    return paths

def crdfs_paths_fromto(source_node, target_node, metapath, exclude_nodes=set(), exclude_edges=set()):
    """
    Cached recursive depth-first-search: computes all paths from
    source_node to target_node of kind metapath. Paths with duplicate
    nodes, with nodes in exclude_nodes, or edges in exclude_edges are excluded.
    Returns of tuple of hetnet.Path() objects.
    """
    paths = list()
    for edge_list in crdfs_paths_from(source_node, metapath):
        if not edge_list[-1].target == target_node:
            continue
        if exclude_edges and exclude_edges & set(edge_list):
            continue
        path = hetnet.Path(edge_list)
        if exclude_nodes and exclude_nodes & set(path.get_nodes()):
            continue
        paths.append(path)
    return tuple(paths)




def rdfs_paths_from(node, metapath):
    """
    Recursive depth-first-search to find all paths corresponding to metapath
    and starting at node. Paths with duplicate nodes are NOT excluded currently.
    Returns a tuple of paths (as tuples of hetnet.Edge() objects
    rather than a hetnet.Path() objects).
    """
    if not metapath:
        return (),
    paths = list()
    metapath_tail = metapath.sub
    for edge in node.edges[metapath[0]]:
        for tail in rdfs_paths_from(edge.target, metapath_tail):
            if node in (e.target for e in tail):
                continue
            paths.append((edge, ) + tail)
    return tuple(paths)


def rdfs_paths_fromto(source_node, target_node, metapath, exclude_nodes=set(), exclude_edges=set()):
    """
    Recursive depth-first-search to find all paths corresponding to metapath,
    originating on source_node and terminating on target_node. Paths with
    duplicate nodes ARE excluded currently. Returns a generator of
    hetnet.Path() objects.
    """
    paths = list()
    for edge_list in rdfs_paths_from(source_node, metapath):
        if not edge_list[-1].target == target_node:
            continue
        if exclude_edges and exclude_edges & set(edge_list):
            continue
        path = hetnet.Path(edge_list)
        if exclude_nodes and exclude_nodes & set(path.get_nodes()):
            continue
        paths.append(path)
    return tuple(paths)


def NPC(paths_s, paths_t):
    if len(paths_t) == 0:
        paths = list()
    else:
        target = paths_t[0].source()
        paths = [path for path in paths_s if path.target() == target]
    denom = len(paths_s) + len(paths_t)
    if denom:
        return 2.0 * len(paths) / denom
    else:
        return None

