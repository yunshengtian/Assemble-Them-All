from random import random
import time

from .utils import irange, argmin, RRT_ITERATIONS, apply_alpha, RED, INF, elapsed_time


class TreeNode(object):

    def __init__(self, config, parent=None):
        self.config = config
        self.parent = parent

    #def retrace(self):
    #    if self.parent is None:
    #        return [self]
    #    return self.parent.retrace() + [self]

    def retrace(self):
        sequence = []
        node = self
        while node is not None:
            sequence.append(node)
            node = node.parent
        return sequence[::-1]

    def clear(self):
        self.node_handle = None
        self.edge_handle = None

    def draw(self, env, color=apply_alpha(RED, alpha=0.5)):
        # https://github.mit.edu/caelan/lis-openrave
        from manipulation.primitives.display import draw_node, draw_edge
        self.node_handle = draw_node(env, self.config, color=color)
        if self.parent is not None:
            self.edge_handle = draw_edge(
                env, self.config, self.parent.config, color=color)

    def __str__(self):
        return 'TreeNode(' + str(self.config) + ')'
    __repr__ = __str__


def configs(nodes):
    if nodes is None:
        return None
    return list(map(lambda n: n.config, nodes))


def rrt(start, goal_sample, distance_fn, sample_fn, extend_fn, collision_fn, goal_test=lambda q: False,
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
    if not callable(goal_sample):
        g = goal_sample
        goal_sample = lambda: g
    nodes = [TreeNode(start)]
    for i in irange(max_iterations):
        if elapsed_time(start_time) >= max_time:
            break
        goal = random() < goal_probability or i == 0
        s = goal_sample() if goal else sample_fn()

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
