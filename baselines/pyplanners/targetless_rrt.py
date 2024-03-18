from random import random
import time
import numpy as np

from .utils import irange, argmin, RRT_ITERATIONS, INF, elapsed_time
from .rrt import TreeNode, configs
from .bbox import compute_nearest_bbox_goal


def targetless_rrt(start, vertices_move, vertices_still, distance_fn, sample_fn, extend_fn, collision_fn, goal_test,
        goal_probability=.2, max_iterations=RRT_ITERATIONS, max_time=INF):
    """
    :param start: Start configuration - conf
    :param distance_fn: Distance function - distance_fn(q1, q2)->float
    :param sample_fn: Sample function - sample_fn()->conf
    :param extend_fn: Extension function - extend_fn(q1, q2)->[q', ..., q"]
    :param collision_fn: Collision function - collision_fn(q)->bool
    :param max_iterations: Maximum number of iterations - int
    :param max_time: Maximum runtime - float
    :return: Path [q', ..., q"] or None if unable to find a solution
    """
    start_time = time.time()
    if collision_fn(start):
        return None
    nodes = [TreeNode(start)]

    for i in irange(max_iterations):
        if elapsed_time(start_time) >= max_time:
            break
        goal = random() < goal_probability or i == 0
        if goal:
            s, last = compute_nearest_bbox_goal(nodes, vertices_move, vertices_still)
        else:
            s = sample_fn()
            last = argmin(lambda n: distance_fn(n.config, s), nodes)
        for q in extend_fn(last.config, s):
            if collision_fn(q):
                break
            last = TreeNode(q, parent=last)
            nodes.append(last)
            if goal_test(last.config):
                return configs(last.retrace())
        else:
            if goal:
                return configs(last.retrace())
    return None
