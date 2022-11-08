from heapq import heappush, heappop
from collections import namedtuple, defaultdict

from .nearest import BruteForceNeighbors, KDNeighbors
from .primitives import default_weights, get_embed_fn, get_distance_fn
from .utils import INF, elapsed_time, get_pairs, default_selector, irange, \
    merge_dicts, compute_path_cost, get_length, is_path, flatten

import time
import numpy as np

Metric = namedtuple('Metric', ['p_norm', 'weights'])
Node = namedtuple('Node', ['g', 'parent'])
Solution = namedtuple('PRMSolution', ['path', 'samples', 'edges', 'colliding_vertices', 'colliding_edges'])

unit_cost_fn = lambda v1, v2: 1.
zero_heuristic_fn = lambda v: 0

ORDINAL = 1e3
REVERSIBLE = True
ROADMAPS = [] # TODO: not ideal

def sample_until(sample_fn, num_samples, max_time=INF):
    # TODO: is this actually needed?
    # TODO: compute number of rejected samples
    start_time = time.time()
    samples = []
    while (len(samples) < num_samples) and (elapsed_time(start_time) < max_time):
        conf = sample_fn()  # TODO: include
        # TODO: bound function based on distance_fn(start, conf) and individual distances
        # if (max_cost == INF) or (distance_fn(start, conf) + distance_fn(conf, goal)) < max_cost:
        # TODO: only keep edges that move toward the goal
        samples.append(conf)
    return samples

def retrace_path(visited, vertex):
    if vertex is None:
        return []
    return retrace_path(visited, visited[vertex].parent) + [vertex]

def dijkstra(start_v, neighbors_fn, cost_fn=unit_cost_fn):
    # Update the heuristic over time
    # TODO: overlap with discrete
    # TODO: all pairs shortest paths
    # TODO: max_cost
    start_g = 0.
    visited = {start_v: Node(start_g, None)}
    queue = [(start_g, start_v)]
    while queue:
        current_g, current_v = heappop(queue)
        if visited[current_v].g < current_g:
            continue
        for next_v in neighbors_fn(current_v):
            next_g = current_g + cost_fn(current_v, next_v)
            if (next_v not in visited) or (next_g < visited[next_v].g):
                visited[next_v] = Node(next_g, current_v)
                heappush(queue, (next_g, next_v))
    return visited

def get_priority_fn(w=1.):
    assert 0. <= w <= INF
    def priority_fn(g, h):
        if w == 0.:
            return g
        if w == INF:
            return (h, g)
        return g + w*h
    return priority_fn

def wastar_search(start_v, end_v, neighbors_fn, cost_fn=unit_cost_fn,
                  heuristic_fn=zero_heuristic_fn, max_cost=INF, max_time=INF, **kwargs):
    # TODO: lazy wastar to get different paths
    # TODO: multi-start / multi-goal
    # TODO: use previous search tree as heuristic
    #heuristic_fn = lambda v: cost_fn(v, end_v)
    priority_fn = get_priority_fn(**kwargs)
    goal_test = lambda v: v == end_v

    start_time = time.time()
    start_g = 0.
    start_h = heuristic_fn(start_v)
    visited = {start_v: Node(start_g, None)}
    queue = [(priority_fn(start_g, start_h), start_g, start_v)]
    while queue and (elapsed_time(start_time) < max_time):
        _, current_g, current_v = heappop(queue)
        if visited[current_v].g < current_g:
            continue
        if goal_test(current_v):
            return retrace_path(visited, current_v)
        for next_v in neighbors_fn(current_v): # TODO: lazily compute neighbors
            next_g = current_g + cost_fn(current_v, next_v)
            if (next_v not in visited) or (next_g < visited[next_v].g):
                visited[next_v] = Node(next_g, current_v)
                next_h = heuristic_fn(next_v)
                if (next_g + next_h) < max_cost: # Assumes admissible
                    next_p = priority_fn(next_g, next_h)
                    heappush(queue, (next_p, next_g, next_v))
    return None

##################################################

class Roadmap(object):
    def __init__(self, extend_fn, weights=None, distance_fn=None, cost_fn=None,
                 p_norm=2, max_degree=5, max_distance=INF, approximate_eps=0., **kwargs):
        # TODO: custom cost_fn
        assert (weights is not None) or (distance_fn is not None)
        self.distance_fn = distance_fn
        self.extend_fn = extend_fn
        self.weights = weights
        self.cost_fn = cost_fn
        self.p_norm = p_norm
        self.max_degree = max_degree
        self.max_distance = max_distance
        self.approximate_eps = approximate_eps
        if self.weights is None:
            self.nearest = BruteForceNeighbors(self.distance_fn)
        else:
            self.nearest = KDNeighbors(embed_fn=get_embed_fn(self.weights), **kwargs)
            #self.nearest = BruteForceNeighbors(get_distance_fn(weights, p_norm=p_norm))
        self.edges = set()
        self.outgoing_from_vertex = defaultdict(set)
        self.edge_costs = {}
        self.edge_paths = {}
        self.colliding_vertices = {}
        self.colliding_edges = {}
        self.colliding_intermediates = {}
    @property
    def samples(self):
        return self.nearest.data
    @property
    def vertices(self):
        return list(range(len(self.samples))) # TODO: inefficient
    def index(self, x):
        #return self.samples.index(x)
        for i, q in enumerate(self.samples):
            if x is q: #x == q:
                return i
        return ValueError(x)
    def __iter__(self):
        return iter([self.samples, self.vertices, self.edges])
    def add_edge(self, v1, v2):
        assert REVERSIBLE # TODO
        edges = {(v1, v2), (v2, v1)}
        self.edges.update(edges)
        self.outgoing_from_vertex[v1].add(v2)
        self.outgoing_from_vertex[v2].add(v1)
        return edges
    def add_samples(self, samples):
        new_edges = set()
        for v1, sample in self.nearest.add_data(samples):
        #for v1, sample in enumerate(self.vertices):
            # TODO: could dynamically compute distances
            # if len(self.outgoing_from_vertex[v1]) >= self.max_degree:
            #     raise NotImplementedError()
            max_degree = min(self.max_degree + 1, len(self.samples))
            for d, v2, _ in self.nearest.query_neighbors(sample, k=max_degree, eps=self.approximate_eps,
                                                         p=self.p_norm, distance_upper_bound=self.max_distance):
                if (v1 != v2): # and (d <= self.max_distance):
                    new_edges.update(self.add_edge(v1, v2))
        return new_edges
    def is_colliding(self, v1, v2):
        edge = (v1, v2)
        return self.colliding_vertices.get(v1, False) or \
               self.colliding_vertices.get(v2, False) or \
               self.colliding_edges.get(edge, False)
    def is_safe(self, v1, v2):
        edge = (v1, v2)
        return not self.colliding_edges.get(edge, True)
    def neighbors_fn(self, v1):
        for v2 in self.outgoing_from_vertex[v1]:
            if not self.is_colliding(v1, v2):
                yield v2
    def check_vertex(self, v, collision_fn):
        x = self.samples[v]
        if v not in self.colliding_vertices:
            # TODO: could update the colliding adjacent edges as well
            self.colliding_vertices[v] = collision_fn(x)
        return not self.colliding_vertices[v]
    def check_intermediate(self, v1, v2, index, collision_fn):
        if (v1, v2, index) not in self.colliding_intermediates:
            x = self.get_path(v1, v2)[index]
            self.colliding_intermediates[v1, v2, index] = collision_fn(x)
            if self.colliding_intermediates[v1, v2, index]:
                # TODO: record when all safe
                self.colliding_edges[v1, v2] = self.colliding_intermediates[v1, v2, index]
                if REVERSIBLE:
                    self.colliding_edges[v2, v1] = self.colliding_edges[v1, v2]
        return not self.colliding_intermediates[v1, v2, index]
    def check_edge(self, v1, v2, collision_fn):
        if (v1, v2) not in self.colliding_edges:
            segment = default_selector(self.get_path(v1, v2)) # TODO: check_intermediate
            self.colliding_edges[v1, v2] = any(map(collision_fn, segment))
            if REVERSIBLE:
                self.colliding_edges[v2, v1] = self.colliding_edges[v1, v2]
        return not self.colliding_edges[v1, v2]
    def check_path(self, path, collision_fn):
        for v in default_selector(path):
            if not self.check_vertex(v, collision_fn):
                return False
        # for v1, v2 in default_selector(get_pairs(path)):
        #     if not self.check_edge(v1, v2, collision_fn):
        #         return False
        # return True
        intermediates = []
        for v1, v2 in get_pairs(path):
            intermediates.extend((v1, v2, index) for index in range(len(self.get_path(v1, v2))))
        for v1, v2, index in default_selector(intermediates):
           if not self.check_intermediate(v1, v2, index, collision_fn):
               return False
        return True
    def check_roadmap(self, collision_fn):
        for vertex in self.vertices:
            self.check_vertex(vertex, collision_fn)
        for vertex1, vertex2 in self.edges:
            self.check_edge(vertex1, vertex2, collision_fn)
    def get_cost(self, v1, v2):
        edge = (v1, v2)
        if edge not in self.edge_costs:
            self.edge_costs[edge] = self.cost_fn(self.samples[v1], self.samples[v2])
            if REVERSIBLE:
                self.edge_costs[edge[::-1]] = self.edge_costs[edge]
        return self.edge_costs[edge]
    def get_path(self, v1, v2):
        edge = (v1, v2)
        if edge not in self.edge_paths:
            path = list(self.extend_fn(self.samples[v1], self.samples[v2]))
            self.edge_paths[edge] = path
            if REVERSIBLE:
                self.edge_paths[edge[::-1]] = path[::-1]
        return self.edge_paths[edge]
    def augment(self, sample_fn, num_samples=100):
        n = len(self.samples)
        if n >= num_samples:
            return self
        samples = sample_until(sample_fn, num_samples=num_samples - n)
        self.add_samples(samples)
        return self

##################################################

def get_metrics(conf, weights=None, p_norm=2, distance_fn=None, cost_fn=None):
    # TODO: can embed pose and/or points on the robot for other distances
    if (weights is None) and (distance_fn is None):
        weights = default_weights(conf, weights=weights)
        #distance_fn = distance_fn_from_extend_fn(extend_fn)
    if cost_fn is None:
        if distance_fn is None:
            cost_fn = get_distance_fn(weights, p_norm=p_norm) # TODO: additive cost, acceleration cost
        else:
            cost_fn = distance_fn
    return weights, distance_fn, cost_fn

def lazy_prm(start, goal, sample_fn, extend_fn, collision_fn, distance_fn=None, cost_fn=None, roadmap=None, num_samples=100,
             weights=None, circular={}, p_norm=2, lazy=True, max_cost=INF, max_time=INF, w=1., meta=False, verbose=True, **kwargs): #, max_paths=INF):
    """
    :param start: Start configuration - conf
    :param goal: End configuration - conf
    :param sample_fn: Sample function - sample_fn()->conf
    :param extend_fn: Extension function - extend_fn(q1, q2)->[q', ..., q"]
    :param collision_fn: Collision function - collision_fn(q)->bool
    :param max_time: Maximum runtime - float
    :param kwargs: Keyword arguments
    :return: Path [q', ..., q"] or None if unable to find a solution
    """
    # TODO: compute hyperparameters using start, goal, and sample_fn statistics
    # TODO: scale default parameters based on
    # TODO: precompute and store roadmap offline
    # TODO: multiple collision functions to allow partial reuse
    # TODO: multi-query motion planning
    start_time = time.time()
    weights, distance_fn, cost_fn = get_metrics(start, weights=weights, p_norm=p_norm, distance_fn=distance_fn, cost_fn=cost_fn)
    if roadmap is None:
        roadmap = Roadmap(extend_fn, weights=weights, distance_fn=distance_fn, cost_fn=cost_fn, circular=circular)
                          #leafsize=10, compact_nodes=True, copy_data=False, balanced_tree=True, boxsize=None, **kwargs)
        roadmap.add_samples([start, goal] + sample_until(sample_fn, num_samples))
    roadmap = roadmap.augment(sample_fn, num_samples=num_samples)

    samples, vertices, edges = roadmap
    start_vertex = roadmap.index(start)
    end_vertex = roadmap.index(goal)
    degree = np.average(list(map(len, roadmap.outgoing_from_vertex.values())))

    # TODO: update collision occupancy based on proximity to existing colliding (for diversity as well)
    # TODO: minimize the maximum distance to colliding
    if not lazy:
        roadmap.check_roadmap(collision_fn)

    weight_fn = roadmap.get_cost
    if meta:
        lazy_fn = lambda v1, v2: int(not roadmap.is_safe(v1, v2)) # TODO: score by length
        #weight_fn = lazy_fn
        #weight_fn = lambda v1, v2: (lazy_fn(v1, v2), cost_fn(samples[v1], samples[v2])) # TODO:
        weight_fn = lambda v1, v2: ORDINAL*lazy_fn(v1, v2) + roadmap.get_cost(v1, v2)

    visited = dijkstra(end_vertex, roadmap.neighbors_fn, weight_fn)
    heuristic_fn = lambda v: visited[v].g if (v in visited) else INF # TODO: lazily apply costs
    #heuristic_fn = zero_heuristic_fn
    #heuristic_fn = lambda v: weight_fn(v, end_vertex)
    path = None
    while (elapsed_time(start_time) < max_time) and (path is None): # TODO: max_attempts
        lazy_path = wastar_search(start_vertex, end_vertex, neighbors_fn=roadmap.neighbors_fn,
                                  cost_fn=weight_fn, heuristic_fn=heuristic_fn,
                                  max_cost=max_cost, max_time=max_time-elapsed_time(start_time), w=w)
        if lazy_path is None:
            break
        if verbose:
            print('Candidate | Length: {} | Cost: {:.3f} | Vertices: {} | Samples: {} | Degree: {:.3f} | Time: {:.3f}'.format(
                len(lazy_path), compute_path_cost(lazy_path, cost_fn=weight_fn),
                len(roadmap.colliding_vertices), len(roadmap.colliding_intermediates),
                degree, elapsed_time(start_time)))
        if roadmap.check_path(lazy_path, collision_fn):
            path = lazy_path

    if path is None:
        forward_visited = set(dijkstra(start_vertex, roadmap.neighbors_fn))
        backward_visited = set(dijkstra(end_vertex, roadmap.neighbors_fn))
        # for v in roadmap.vertices:
        #     if not roadmap.colliding_vertices.get(v, False):
        #         # TODO: add edges if the collision-free degree drops
        #         num_colliding = sum(not roadmap.colliding_vertices.get(v2, False) for v2 in roadmap.outgoing_from_vertex[v])

        if verbose:
            print('Failure | Forward: {} | Backward: {} | Vertices: {} | Samples: {} | Degree: {:.3f} | Time: {:.3f}'.format(
                len(forward_visited), len(backward_visited),
                len(roadmap.colliding_vertices), len(roadmap.colliding_intermediates),
                degree, elapsed_time(start_time)))
        return Solution(path, samples, edges, roadmap.colliding_vertices, roadmap.colliding_edges)

    if verbose:
        print('Solution | Length: {} | Cost: {:.3f} | Vertices: {} | Samples: {} | Degree: {:.3f} | Time: {:.3f}'.format(
            len(path), compute_path_cost(path, cost_fn=weight_fn),
            len(roadmap.colliding_vertices), len(roadmap.colliding_intermediates),
            degree, elapsed_time(start_time)))
    #waypoints = [samples[v] for v in path]
    #solution = [start] + refine_waypoints(waypoints, extend_fn)
    solution = [start] + list(flatten(roadmap.get_path(v1, v2) for v1, v2 in get_pairs(path)))
    return Solution(solution, samples, edges, roadmap.colliding_vertices, roadmap.colliding_edges)

##################################################

def create_param_sequence(initial_samples=100, step_samples=100, **kwargs):
    # TODO: iteratively increase the parameters
    # TODO: generalize to degree, distance, cost
    return (merge_dicts(kwargs, {'num_samples': num_samples})
            for num_samples in irange(start=initial_samples, stop=INF, step=step_samples))

def lazy_prm_star(start, goal, sample_fn, extend_fn, collision_fn, distance_fn=None, cost_fn=None, max_cost=INF, success_cost=INF,
                  param_sequence=None, resuse=True, weights=None, circular={}, p_norm=2, max_time=INF, verbose=False, **kwargs):
    # TODO: bias to stay near the (past/hypothetical) path
    # TODO: proximity pessimistic collision checking
    # TODO: roadmap reuse in general
    start_time = time.time()
    weights, distance_fn, cost_fn = get_metrics(start, weights=weights, p_norm=p_norm, distance_fn=distance_fn, cost_fn=cost_fn)
    #print(weights, distance_fn, cost_fn)
    #input()
    if param_sequence is None:
        param_sequence = create_param_sequence()

    roadmap = None
    best_path = None
    best_cost = max_cost
    for i, params in enumerate(param_sequence):
        remaining_time = (max_time - elapsed_time(start_time))
        if remaining_time <= 0.:
            break
        if verbose:
            print('\nIteration: {} | Cost: {:.3f} | Elapsed: {:.3f} | Remaining: {:.3f} | Params: {}'.format(
                i, best_cost, elapsed_time(start_time), remaining_time, params))
        if (roadmap is None) or not resuse:
            roadmap = Roadmap(extend_fn, weights=weights, distance_fn=distance_fn, cost_fn=cost_fn, circular=circular)
            roadmap.add_samples([start, goal] + sample_until(sample_fn, params['num_samples']))

        new_path = lazy_prm(start, goal, sample_fn, extend_fn, collision_fn, roadmap=roadmap,
                            cost_fn=cost_fn, weights=weights, circular=circular, p_norm=p_norm,
                            max_time=remaining_time, max_cost=best_cost,
                            verbose=verbose, **params, **kwargs)[0]
        new_cost = compute_path_cost(new_path, cost_fn=cost_fn)
        if verbose:
            print('Path: {} | Cost: {:.3f} | Length: {}'.format(is_path(new_path), new_cost, get_length(new_path)))
        if new_cost < best_cost:
            best_path = new_path
            best_cost = new_cost
        if best_cost < success_cost:
            break
    if roadmap is not None:
        ROADMAPS.append(roadmap)
    return best_path
