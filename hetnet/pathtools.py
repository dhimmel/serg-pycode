import hetnet

import cachetools

"""
def recursive_path_search(graph, source_node, metapath):
    sub_metapath = metapath.sub
    if sub_metapath is None:
        return list()
    edge_lists = list()
    for edge in source_node.edges(metapath[0]):
        subpaths = recursive_path_search(graph, edge.target, sub_metapath)
        for subpath in subpaths:
            edge_list = [edge] + subpath
            edge_lists.append(edge_list)
    return edge_lists





def recursive_path_search(graph, source_node, metapath):
    sub_metapath = metapath.sub
    if sub_metapath is not None:
        for edge in source_node.edges[metapath[0]]:
            subpaths = recursive_path_search(graph, edge.target, sub_metapath)
            for subpath in subpaths:
                yield [edge] + subpath
"""






def rdfs_paths_from(node, metapath, preceding_edges=tuple()):
    """
    Recursive depth-first-search to find all paths corresponding to metapath
    and starting at node. Paths with duplicate nodes are NOT excluded currently.
    Returns a generator that yields paths (as tuples of hetnet.Edge() objects
    rather than a hetnet.Path() objects).
    """
    if not metapath:
        yield preceding_edges
    else:
        sub_metapath = metapath.sub
        for edge in node.edges[metapath[0]]:
            succeeding_paths = rdfs_paths_from(
                edge.target, sub_metapath, preceding_edges + (edge, ))
            for succeeding_path in succeeding_paths:
                yield succeeding_path

def rdfs_paths_fromto(source_node, target_node, metapath):
    """
    Recursive depth-first-search to find all paths corresponding to metapath,
    originating on source_node and terminating on target_node. Paths with
    duplicate nodes ARE excluded currently. Returns a generator of
    hetnet.Path() objects.
    """
    for edge_list in rdfs_paths_from(source_node, metapath):
        if not edge_list[-1] == target_node:
            continue
        path = hetnet.Path(edge_list)
        nodes = path.get_nodes()
        if len(set(nodes)) > len(nodes):
            continue
        yield path


@cachetools.LFUCache
def cached_rdfs_paths_from(get_metapath, node, metapath, preceding_edges=tuple()):
    """
    CACHED: Recursive depth-first-search to find all paths corresponding to metapath
    and starting at node. Paths with duplicate nodes are NOT excluded currently.
    Returns a generator that yields paths (as tuples of hetnet.Edge() objects
    rather than a hetnet.Path() objects).
    """
    if metapath:
        metapath_to_paths = node.metapaths
        metapath_head = metapath.max_overlap(metapath_to_paths.keys())
        if metapath_head:
            for path in metapath_to_paths[metapath_head]:
                head_len = len(metapath_head)
                metapath_tail = get_metapath(metapath[:head_len]) if head_len else None
                for succeeding_path in rdfs_paths_from(
                        path.target(), metapath_tail, preceding_edges + path.edges):
                    yield succeeding_path
        else:
            sub_metapath = metapath.sub
            for edge in node.edges[metapath[0]]:
                succeeding_paths = rdfs_paths_from(
                    edge.target, sub_metapath, preceding_edges + (edge, ))
                for succeeding_path in succeeding_paths:
                    yield succeeding_path
    else:
        yield preceding_edges


def compute_metapaths(graph, metapaths):










#import psutil
import resource
import sys
import collections
import itertools
import functools

"""
def memory_usage_psutil():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.get_memory_info()[0] / float(2 ** 20)
    return mem
"""

def memory_usage_resource():
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem



class memoize_lru(object):
    '''Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    https://wiki.python.org/moin/PythonDecoratorLibrary#Memoize
    '''
    def __init__(self, func):
        max_GB = 100.0
        remove = 0.2
        self.func = func
        self.cache = collections.OrderedDict()
        self.max_MB = max_GB * 1024.0
        self.remove = remove

    def __call__(self, *args):
        if args in self.cache:
            value = self.cache.pop(args)
            self.cache[args] = value
            return value

        else:
            value = self.func(*args)
            cache = self.cache
            cache[args] = value
            if memory_usage_resource() > self.max_MB:
                n_remove = int(len(cache) * self.remove)
                remove_keys = list(itertools.tee(cache.iterkeys(), n_remove))
                for key in remove_keys:
                    del cache[key]
            return value



cache = collections.OrderedDict()
cache_gets = 0
cache_sets = 0
max_MB = 50 * 1024
prune_fraction = 0.2

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
    if memory_usage_resource() > max_MB:
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


def cached_rdfs_paths_from(node, metapath):
    """
    Recursive depth-first-search to find all paths corresponding to metapath
    and starting at node. Paths with duplicate nodes are NOT excluded currently.
    Returns a generator that yields paths (as tuples of hetnet.Edge() objects
    rather than a hetnet.Path() objects).
    """
    if not metapath:
        return tuple(),
    args = node, metapath
    if args in cache:
        return cache_get(args)
    paths = list()
    metapath_tail = metapath.sub
    for edge in node.edges[metapath[0]]:
        for tail in rdfs_paths_from(edge.target, metapath_tail):
            if node in (e.target for e in tail):
                continue
            paths.append((edge, ) + tail)
    paths = tuple(paths)
    cache_set(args, paths)
    return paths



def rdfs_paths_from(node, metapath):
    """
    Recursive depth-first-search to find all paths corresponding to metapath
    and starting at node. Paths with duplicate nodes are NOT excluded currently.
    Returns a generator that yields paths (as tuples of hetnet.Edge() objects
    rather than a hetnet.Path() objects).
    """
    if not metapath:
        return tuple(),
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
    return paths



graph = hetnet.readwrite.graph.read_pickle('/home/dhimmels/Documents/serg/gene-disease-hetnet/networks/140518-metricsweep/graph/graph.pkl.gz' )
metagraph = graph.metagraph

metapaths = metagraph.extract_metapaths('gene', 'disease', 5)
node_IRF1 = graph.node_dict['IRF1']
node_CXCR4 = graph.node_dict['CXCR4']
node_MS = graph.node_dict['DOID:2377']

rdfs_paths_fromto(node_IRF1, node_MS, metapaths[1])

print('Initial Memory Usage: {:.1f}. Max Memory Usage: {:.1f}'.format(
    memory_usage_resource() / 1024.0, max_MB / 1024.0))
for metapath in metapaths:
    paths = rdfs_paths_from(node_CXCR4, metapath)
    print metapath, len(paths)










