from itertools import takewhile

import numpy as np

from .rrt import TreeNode
from .trajectory.linear import get_default_limits, solve_linear
from .trajectory.retime import spline_duration
from .utils import argmin, negate, circular_difference, UNBOUNDED_LIMITS, get_distance, get_delta

ASYMETRIC = True


def asymmetric_extend(q1, q2, extend_fn, backward=False):
    if backward and ASYMETRIC:
        return reversed(list(extend_fn(q2, q1))) # Forward model
    return extend_fn(q1, q2)


def extend_towards(tree, target, distance_fn, extend_fn, collision_fn, swap=False, tree_frequency=1, **kwargs):
    assert tree_frequency >= 1
    last = argmin(lambda n: distance_fn(n.config, target), tree)
    extend = list(asymmetric_extend(last.config, target, extend_fn, backward=swap))
    safe = list(takewhile(negate(collision_fn), extend))
    for i, q in enumerate(safe):
        if (i % tree_frequency == 0) or (i == len(safe) - 1):
            last = TreeNode(q, parent=last)
            tree.append(last)
    success = len(extend) == len(safe)
    return last, success

##################################################

def calculate_radius(d=2):
    # TODO: unify with get_threshold_fn
    # Sampling-based Algorithms for Optimal Motion Planning
    # http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.419.5503&rep=rep1&type=pdf
    # https://en.wikipedia.org/wiki/Volume_of_an_n-ball
    interval = (1 - 0)
    vol_free = interval ** d
    radius = 1./2
    vol_ball = np.pi * (radius ** d)
    gamma = 2 * ((1 + 1. / d) * (vol_free / vol_ball)) ** (1. / d)
    # threshold = gamma * (math.log(n) / n) ** (1. / d)
    return gamma


def default_weights(conf, weights=None, scale=1.):
    if weights is not None:
        return weights
    d = len(conf)
    weights = scale*np.ones(d)
    return weights


def get_embed_fn(weights):
    weights = np.array(weights)
    return lambda q: weights * q


def get_distance_fn(weights, p_norm=2):
    embed_fn = get_embed_fn(weights)
    return lambda q1, q2: np.linalg.norm(embed_fn(q2) - embed_fn(q1), ord=p_norm)


def distance_fn_from_extend_fn(extend_fn):
    # TODO: can compute cost between waypoints from extend_fn
    def distance_fn(q1, q2):
        path = list(extend_fn(q1, q2)) # TODO: cache
        return len(path) # TODO: subtract endpoints?
    return distance_fn

##################################################

def get_difference_fn(circular={}):
    def fn(q2, q1):
        return tuple(circular_difference(v2, v1, interval=circular.get(i, UNBOUNDED_LIMITS))
                     for i, (v2, v1) in enumerate(zip(q2, q1)))
    return fn


def get_cost_fn(distance_fn=get_distance, constant=0., coefficient=1.):
    def fn(q1, q2):
        return constant + coefficient*distance_fn(q1, q2)
    return fn


def get_duration_fn(difference_fn=get_delta, t_constant=0., t_min=0., **kwargs):
    v_max, a_max = get_default_limits(d=None, **kwargs)
    def fn(q1, q2):
        # TODO: be careful that not colinear with other waypoints
        difference = difference_fn(q1, q2)
        t_transit = 0.
        if not np.allclose(np.zeros(len(difference)), difference, atol=1e-6, rtol=0):
            t_transit = solve_linear(difference, v_max, a_max, only_duration=True)
            assert t_transit is not None
            #curve = solve_linear(difference, v_max, a_max)
            #t_transit = spline_duration(curve)
        t = t_constant + t_transit
        return max(t_min, t) # TODO: clip function
    return fn