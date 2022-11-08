from itertools import islice

import numpy as np

from .search import best_first, bfs


def get_nth(generator, n=0):
    return next(islice(generator, n), None)


def get_neighbors_fn(extend_fn, targets=[], scale=1e3, bias=False): # TODO: could also include diagonal
    # https://github.mit.edu/caelan/lis-openrave/blob/master/manipulation/motion/cspace.py#L171
    def neighbors_fn(current):
        d = len(current)
        for target in targets:
            new = get_nth(extend_fn(current, target), n=1)
            if bias or (new is target):
                yield new
        for k in range(d):
            direction = np.zeros(d)
            direction[k] = scale
            for sign in [-1, +1]:
                # TODO: hash the confs
                target = tuple(np.array(current) + sign * direction)
                #yield target
                new = get_nth(extend_fn(current, target), n=1)
                yield new
    return neighbors_fn


def lattice(start, goal, extend_fn, collision_fn, distance_fn=None, **kwargs):
    """
    :param start: Start configuration - conf
    :param goal: End configuration - conf
    :param distance_fn: Distance function - distance_fn(q1, q2)->float
    :param extend_fn: Extension function - extend_fn(q1, q2)->[q', ..., q"]
    :param collision_fn: Collision function - collision_fn(q)->bool
    :param kwargs: Keyword arguments
    :return: Path [q', ..., q"] or None if unable to find a solution
    """
    #collision_fn = lambda q: False
    neighbors_fn = get_neighbors_fn(extend_fn, targets=[goal])
    if distance_fn is None:
        return bfs(start, goal, neighbors_fn, collision_fn, **kwargs)
    return best_first(start, goal, distance_fn, neighbors_fn, collision_fn, **kwargs)
