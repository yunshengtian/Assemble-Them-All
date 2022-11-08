from random import shuffle
from itertools import islice
from collections import deque, defaultdict, namedtuple
import time
import contextlib
import pstats
import cProfile
import random

import numpy as np

INF = float('inf')
PI = np.pi


# TODO: deprecate these defaults
RRT_ITERATIONS = 20
RRT_RESTARTS = 2
RRT_SMOOTHING = 20


try:
   user_input = raw_input
except NameError:
   user_input = input


RED = (1, 0, 0)
GREEN = (0, 1, 0)
BLUE = (0, 0, 1)


Interval = namedtuple('Interval', ['lower', 'upper']) # AABB
UNIT_LIMITS = Interval(0., 1.)
CIRCULAR_LIMITS = Interval(-PI, PI)
UNBOUNDED_LIMITS = Interval(-INF, INF)

##################################################

def apply_alpha(color, alpha=1.):
   return tuple(color[:3]) + (alpha,)


def irange(start, stop=None, step=1):  # np.arange
    if stop is None:
        stop = start
        start = 0
    while start < stop:
        yield start
        start += step


def negate(test):
    return lambda *args, **kwargs: not test(*args, **kwargs)


def clip(value, min_value=-INF, max_value=+INF):
    return min(max(min_value, value), max_value)


def argmin(function, sequence):
    # TODO: use min
    values = list(sequence)
    scores = [function(x) for x in values]
    return values[scores.index(min(scores))]


def get_pairs(sequence):
    sequence = list(sequence)
    return list(zip(sequence[:-1], sequence[1:]))


def merge_dicts(*args):
    result = {}
    for d in args:
        result.update(d)
    return result
    # return dict(reduce(operator.add, [d.items() for d in args]))


def flatten(iterable_of_iterables):
    return (item for iterables in iterable_of_iterables for item in iterables)


def randomize(sequence):
    sequence = list(sequence)
    shuffle(sequence)
    return sequence


def is_even(num):
    return num % 2 == 0


def is_odd(num):
    return num % 2 == 1


def bisect(sequence):
    sequence = list(sequence)
    indices = set()
    queue = deque([(0, len(sequence)-1)])
    while queue:
        lower, higher = queue.popleft()
        if lower > higher:
            continue
        index = int((lower + higher) / 2.)
        assert index not in indices
        #if is_even(higher - lower):
        yield sequence[index]
        queue.extend([
            (lower, index-1),
            (index+1, higher),
        ])


def take(iterable, n=INF):
    if n == INF:
        n = None  # NOTE - islice takes None instead of INF
    elif n is None:
        n = 0  # NOTE - for some of the uses
    return islice(iterable, n)


def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    enums['names'] = sorted(enums.keys(), key=lambda k: enums[k])
    return type('Enum', (), enums)


def elapsed_time(start_time):
    return time.time() - start_time


@contextlib.contextmanager
def profiler(field='tottime', num=10):
    pr = cProfile.Profile()
    pr.enable()
    yield
    pr.disable()
    pstats.Stats(pr).sort_stats(field).print_stats(num) # cumtime | tottime


def inf_sequence():
    return iter(int, 1)


def find(test, sequence):
    for item in sequence:
        if test(item):
            return item
    raise RuntimeError()


def get_sign(x):
    if x > 0:
        return +1
    if x < 0:
        return -1
    return x


def strictly_increasing(sequence):
    return all(x2 > x1 for x1, x2 in get_pairs(sequence))

##################################################

def get_delta(q1, q2):
    return np.array(q2) - np.array(q1)


def get_difference(q2, q1):
    return get_delta(q1, q2)


def get_distance(q1, q2):
    return np.linalg.norm(get_delta(q1, q2))


def get_unit_vector(vec):
    norm = np.linalg.norm(vec)
    if norm == 0:
        return vec
    return np.array(vec) / norm


def is_path(path):
    return path is not None


def compute_path_cost(path, cost_fn=get_distance):
    if not is_path(path):
        return INF
    #path = waypoints_from_path(path)
    return sum(cost_fn(*pair) for pair in get_pairs(path))


def get_length(path):
    if not is_path(path):
        return INF
    return len(path)


def remove_redundant(path, tolerance=1e-3):
    assert path
    new_path = [path[0]]
    for conf in path[1:]:
        difference = get_difference(new_path[-1], np.array(conf))
        if not np.allclose(np.zeros(len(difference)), difference, atol=tolerance, rtol=0):
            new_path.append(conf)
    return new_path


def waypoints_from_path(path, difference_fn=None, tolerance=1e-3):
    if difference_fn is None:
        difference_fn = get_difference
    path = remove_redundant(path, tolerance=tolerance)
    if len(path) < 2:
        return path
    waypoints = [path[0]]
    last_conf = path[1]
    last_difference = get_unit_vector(difference_fn(last_conf, waypoints[-1]))
    for conf in path[2:]:
        difference = get_unit_vector(difference_fn(conf, waypoints[-1]))
        if not np.allclose(last_difference, difference, atol=tolerance, rtol=0):
            waypoints.append(last_conf)
            difference = get_unit_vector(difference_fn(conf, waypoints[-1]))
        last_conf = conf
        last_difference = difference
    waypoints.append(last_conf)
    return waypoints


def refine_waypoints(waypoints, extend_fn):
    #if len(waypoints) <= 1:
    #    return waypoints
    return list(flatten(extend_fn(q1, q2) for q1, q2 in get_pairs(waypoints))) # [waypoints[0]] +

##################################################

def convex_combination(x, y, w=0.5):
    return (1-w)*np.array(x) + w*np.array(y)


def uniform_generator(d):
    while True:
        yield np.random.uniform(size=d)


def halton_generator(d, seed=None):
    # TODO: randomly sample an initial point and then wrap around
    # TODO: apply random noise on top
    # https://ghalton.readthedocs.io/en/latest/
    import ghalton
    if seed is None:
        seed = random.randint(0, 100-1) # ghalton.EA_PERMS[d-1]
    #ghalton.PRIMES, ghalton.EA_PERMS
    #sequencer = ghalton.Halton(d)
    #sequencer = ghalton.GeneralizedHalton(d, seed) # TODO: seed not working
    sequencer = ghalton.GeneralizedHalton(ghalton.EA_PERMS[:d])
    #sequencer.reset()
    #sequencer.seed(seed) # TODO: seed not working
    sequencer.get(seed) # Burn this number of values
    while True:
        [weights] = sequencer.get(1)
        yield np.array(weights)


def unit_generator(d, use_halton=False, **kwargs):
    # TODO: mixture generator
    if use_halton:
        try:
            import ghalton
        except ImportError:
            print('ghalton is not installed (https://pypi.org/project/ghalton/)')
            use_halton = False
    return halton_generator(d, **kwargs) if use_halton else uniform_generator(d)


def interval_generator(lower, upper, **kwargs):
    assert len(lower) == len(upper)
    assert np.less_equal(lower, upper).all()
    if np.equal(lower, upper).all():
        return iter([lower])
    return (convex_combination(lower, upper, w=weights)
            for weights in unit_generator(d=len(lower), **kwargs))

##################################################

def forward_selector(path):
    return path


def backward_selector(path):
    return reversed(list(path))


def random_selector(path):
    return randomize(path)


def bisect_selector(path):
    return bisect(path)


default_selector = bisect_selector # random_selector

##################################################

def is_hashable(value):
    #return isinstance(value, Hashable) # TODO: issue with hashable and numpy 2.7.6
    try:
        hash(value)
    except TypeError:
        return False
    return True

def value_or_id(value):
    if is_hashable(value):
        return value
    return id(value)

##################################################

def incoming_from_edges(edges):
    incoming_vertices = defaultdict(set)
    for v1, v2 in edges:
        incoming_vertices[v2].add(v1)
    return incoming_vertices


def outgoing_from_edges(edges):
    # neighbors_from_index = {v: set() for v in vertices}
    # for v1, v2 in edges:
    #     neighbors_from_index[v1].add(v2)
    outgoing_vertices = defaultdict(set)
    for v1, v2 in edges:
        outgoing_vertices[v1].add(v2)
    return outgoing_vertices


def adjacent_from_edges(edges):
    undirected_edges = defaultdict(set)
    for v1, v2 in edges:
        undirected_edges[v1].add(v2)
        undirected_edges[v2].add(v1)
    return undirected_edges

##################################################

def normalize_interval(value, interval=UNIT_LIMITS):
    # TODO: move more out of pybullet-planning
    lower, upper = interval
    assert lower <= upper
    return (value - lower) / (upper - lower)


def rescale_interval(value, old_interval=UNIT_LIMITS, new_interval=UNIT_LIMITS):
    lower, upper = new_interval
    return convex_combination(lower, upper, w=normalize_interval(value, old_interval))


def wrap_interval(value, interval=UNIT_LIMITS):
    lower, upper = interval
    if (lower == -INF) and (+INF == upper):
        return value
    assert -INF < lower <= upper < +INF
    return (value - lower) % (upper - lower) + lower


def interval_distance(value1, value2, interval=UNIT_LIMITS):
    value1 = wrap_interval(value1, interval)
    value2 = wrap_interval(value2, interval)
    if value1 > value2:
        value1, value2 = value2, value1
    lower, upper = interval
    return min(value2 - value1, (value1 - lower) + (upper - value2))


def circular_difference(theta2, theta1, interval=UNIT_LIMITS):
    extent = get_interval_extent(interval)
    diff_interval = Interval(-extent/2, +extent/2)
    return wrap_interval(theta2 - theta1, interval=diff_interval)

##################################################

def get_interval_center(interval):
    lower, upper = interval
    return np.average([lower, upper], axis=0)


def get_interval_extent(interval):
    lower, upper = interval
    return get_delta(lower, upper)


def even_space(start, stop, step=1, endpoint=True):
    sequence = np.arange(start, stop, step=step)
    if not endpoint:
        return sequence
    return np.append(sequence, [stop])