import numpy



def random_walk(r, seed_vector, adj_matrix, iteration_cutoff = 10e-10):
    """
    Random walk with restart
    http://dx.doi.org/10.1016/j.ajhg.2008.02.013
    """
    # column normalize adj_matrix
    W = adj_matrix / adj_matrix.sum(axis=0)
    p0 = seed_vector / sum(seed_vector)
    p0 = p0[:, None]
    pt = p0
    #l1_norm_t = numpy.linalg.norm(pt, 1)
    pt1 = None
    steps = 0
    while True:
        pt1 = random_walk_step(r, p0, W, pt)
        #l1_norm_t1 = numpy.linalg.norm(pt1, 1)
        change = max(abs(pt1 - pt))
        #change = abs(l1_norm_t - l1_norm_t1)
        #print change
        pt = pt1
        #l1_norm_t = l1_norm_t1
        steps += 1
        if change < iteration_cutoff:
            return pt, steps


def random_walk_step(r, p0, W, pt):
    """
    r = probability of restart
    p0 = initial probability vector
    W = column-normalized adjacency matrix
    pt = probability vector at time step t
    pt = probability vector at time step t + 1
    """
    #print 'r', r
    #print 'W', W
    #print 'pt', pt
    #print 'p0', p0
    #pt1 = (1 - r) * W.dot(pt) + r * p0
    pt1 = (1 - r) * W * pt + r * p0
    return pt1

