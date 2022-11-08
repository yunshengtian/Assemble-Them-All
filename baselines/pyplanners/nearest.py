from scipy.spatial import KDTree
from itertools import product
from .utils import get_interval_extent, UNBOUNDED_LIMITS, INF

import numpy as np

def expand_circular(x, circular={}):
    domains = []
    for k in range(len(x)):
        interval = circular.get(k, UNBOUNDED_LIMITS)
        extent = get_interval_extent(interval)
        if extent != INF:
            domains.append([
                -extent, 0., +extent, # TODO: choose just one
            ])
        else:
            domains.append([0.])
    for dx in product(*domains):
        wx = x + np.array(dx)
        yield wx

##################################################

class NearestNeighbors(object):
    def __init__(self):
        pass
        # self.data = []
        # self.add_data(data)
    def add_data(self, new_data):
        raise NotImplementedError()
    def query_neighbors(self, x, k=1, **kwargs):
        raise NotImplementedError()

##################################################

class KDNeighbors(NearestNeighbors):
    # https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KDTree.html
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html
    # TODO: approximate KDTrees
    # https://github.com/lmcinnes/pynndescent
    # https://github.com/spotify/annoy
    # https://github.com/flann-lib/flann
    def __init__(self, data=[], circular={}, embed_fn=lambda x: x, **kwargs): # [0, 1]
        super(NearestNeighbors, self).__init__()
        # TODO: maintain tree and brute-force list
        self.data = [] # TODO: self.kd_tree.data
        self.embedded = []
        self.kd_tree = None
        self.circular = circular
        self.embed_fn = embed_fn
        self.kwargs = kwargs
        self.stale = True
        self.add_data(data)
    def update(self):
        if not self.stale:
            return
        self.stale = False
        if self.embedded:
            self.kd_tree = KDTree(self.embedded,
                                  #leafsize=10, compact_nodes=True, copy_data=False, balanced_tree=True, boxsize=None
                                  **self.kwargs)
    def add_data(self, new_data):
        indices = []
        for x in new_data:
            index = len(self.data)
            indices.append(index)
            self.data.append(x)
            self.embedded.append(self.embed_fn(x))
        self.stale |= bool(indices)
        self.update()
        return zip(indices, new_data)
    def remove_data(self, new_data):
        raise NotImplementedError() # TODO: need to keep track of data indices (using id?)
    def query_neighbors(self, x, k=1, **kwargs):
        # TODO: class **kwargs
        closest_neighbors = {}
        for wx in expand_circular(x, circular=self.circular):
            embedded = self.embed_fn(wx)
            #print(x, embedded)
            # k=1, eps=0, p=2, distance_upper_bound=inf, workers=1
            for d, i in zip(*self.kd_tree.query(embedded, k=k, **kwargs)):
                if d < closest_neighbors.get(i, INF):
                    closest_neighbors[i] = d
        return [(d, i, self.data[i]) for i, d in sorted(
            closest_neighbors.items(), key=lambda pair: pair[1])][:k] # TODO: filter

##################################################

class BruteForceNeighbors(NearestNeighbors):
    def __init__(self, distance_fn, data=[], **kwargs):
        super(BruteForceNeighbors, self).__init__()
        self.distance_fn = distance_fn
        self.data = []
        self.add_data(data)
    def add_data(self, new_data):
        indices = []
        for x in new_data:
            index = len(self.data)
            indices.append(index)
            self.data.append(x)
        return zip(indices, new_data)
    def query_neighbors(self, x, k=1, **kwargs):
        # TODO: store pairwise distances
        neighbors = []
        for i, x2 in enumerate(self.data):
            d = self.distance_fn(x, x2)
            neighbors.append((d, i, x2))
        # TODO: heapq
        return sorted(neighbors, key=lambda pair: pair[0])[:k]
