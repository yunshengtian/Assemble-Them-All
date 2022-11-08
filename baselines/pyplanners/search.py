from collections import deque, namedtuple
from heapq import heappop, heappush

import numpy as np
import time

from .utils import INF, elapsed_time

# https://github.mit.edu/caelan/lis-openrave/tree/master/manipulation/motion
# https://github.mit.edu/caelan/lis-openrave/commit/4d8683407ec79a7c39dab62d6779804730ff598d

Node = namedtuple('Node', ['g', 'parent'])


def retrace(visited, q):
    if q is None:
        return []
    return retrace(visited, visited[tuple(q)].parent) + [q]


def bfs(start, goal, neighbors_fn, collision_fn, max_iterations=INF, max_time=INF):
    """
    :param start: Start configuration - conf
    :param goal: End configuration - conf
    :param collision_fn: Collision function - collision_fn(q)->bool
    :param max_time: Maximum runtime - float
    :return: Path [q', ..., q"] or None if unable to find a solution
    """
    start_time = time.time()
    if collision_fn(start) or collision_fn(goal):
        return None
    iterations = 0
    visited = {tuple(start): Node(g=0, parent=None)}
    queue = deque([start])
    while queue and (iterations < max_iterations) and (elapsed_time(start_time) < max_time):
        iterations += 1
        current = queue.popleft()
        if goal is not None and tuple(current) == tuple(goal):
            return retrace(visited, current)
        for new in neighbors_fn(current):
            # TODO - make edges for real (and store bad edges)
            if (tuple(new) not in visited) and not collision_fn(new):
                visited[tuple(new)] = Node(visited[tuple(current)].g + 1, current)
                queue.append(new)
    return None

##################################################

def weighted(weight=1.):
    if weight == INF:
        return lambda g, h: h
    return lambda g, h: g + weight*h

uniform = weighted(0)
astar = weighted(1)
wastar2 = weighted(2)
wastar3 = weighted(2)
greedy = weighted(INF)
lexicographic = lambda g, h: (h, g)

def best_first(start, goal, distance_fn, neighbors_fn, collision_fn,
               max_iterations=INF, max_time=INF, priority=lexicographic):  # TODO - put start and goal in neighbors_fn
    """
    :param start: Start configuration - conf
    :param goal: End configuration - conf
    :param distance_fn: Distance function - distance_fn(q1, q2)->float
    :param collision_fn: Collision function - collision_fn(q)->bool
    :param max_time: Maximum runtime - float
    :return: Path [q', ..., q"] or None if unable to find a solution
    """
    start_time = time.time()
    if collision_fn(start) or collision_fn(goal):
        return None
    queue = [(priority(0, distance_fn(start, goal)), 0, start)]
    visited = {tuple(start): Node(g=0, parent=None)}
    iterations = 0
    while queue and (iterations < max_iterations) and (elapsed_time(start_time) < max_time):
        _, current_g, current = heappop(queue)
        current = np.array(current)
        if visited[tuple(current)].g != current_g:
            continue
        # TODO: lazy collision_fn checking
        iterations += 1
        if tuple(current) == tuple(goal):
            return retrace(visited, current)
        for new in neighbors_fn(current):
            new_g = current_g + distance_fn(current, new)
            if (tuple(new) not in visited or new_g < visited[tuple(new)].g) and not collision_fn(new):
                visited[tuple(new)] = Node(new_g, current)
                # ValueError: The truth value of an array with more than one
                # element is ambiguous.
                heappush(queue, (priority(new_g, distance_fn(new, goal)), new_g, new))
    return None
