try:
    from Tkinter import Tk, Canvas, Toplevel, LAST
    #import TKinter as tk
except ModuleNotFoundError:
    from tkinter import Tk, Canvas, Toplevel, LAST
    #import tkinter as tk
from _tkinter import TclError

import numpy as np
import traceback

from collections import namedtuple

from ..utils import get_pairs, get_delta, INF, get_interval_center, get_interval_extent

Box = namedtuple('Box', ['lower', 'upper'])
Circle = namedtuple('Circle', ['center', 'radius'])

class PRMViewer(object):
    def __init__(self, width=500, height=500, title='PRM', background='tan'):
        # TODO: matplotlib viewer
        tk = Tk()
        tk.withdraw()
        top = Toplevel(tk)
        top.wm_title(title)
        top.protocol('WM_DELETE_WINDOW', top.destroy)
        self.width = width
        self.height = height
        self.canvas = Canvas(top, width=self.width, height=self.height, background=background)
        self.canvas.pack()

    def pixel_from_point(self, point):
        (x, y) = point
        # return (int(x*self.width), int(self.height - y*self.height))
        return (x * self.width, self.height - y * self.height)

    def draw_point(self, point, radius=5, color='black'):
        (x, y) = self.pixel_from_point(point)
        self.canvas.create_oval(x - radius, y - radius, x + radius, y + radius, fill=color, outline='')

    def draw_line(self, segment, color='black'):
        (point1, point2) = segment
        (x1, y1) = self.pixel_from_point(point1)
        (x2, y2) = self.pixel_from_point(point2)
        self.canvas.create_line(x1, y1, x2, y2, fill=color, width=1)

    def draw_arrow(self, point1, point2, color='black'):
        (x1, y1) = self.pixel_from_point(point1)
        (x2, y2) = self.pixel_from_point(point2)
        self.canvas.create_line(x1, y1, x2, y2, fill=color, width=2, arrow=LAST)

    def draw_rectangle(self, box, width=2, color='brown'):
        (point1, point2) = box
        (x1, y1) = self.pixel_from_point(point1)
        (x2, y2) = self.pixel_from_point(point2)
        self.canvas.create_rectangle(x1, y1, x2, y2, outline='black', fill=color, width=width)

    def draw_circle(self, center, radius, width=2, color='black'):
        (x1, y1) = self.pixel_from_point(np.array(center) - radius * np.ones(2))
        (x2, y2) = self.pixel_from_point(np.array(center) + radius * np.ones(2))
        self.canvas.create_oval(x1, y1, x2, y2, outline='black', fill=color, width=width)

    def clear(self):
        self.canvas.delete('all')


#################################################################

STEP_SIZE = 1e-2
MIN_PROXIMITY = 1e-3

def contains_box(point, box, buffer=0.):
    (lower, upper) = box
    lower = lower - buffer*np.ones(len(lower))
    upper = upper + buffer*np.ones(len(upper))
    return np.less_equal(lower, point).all() and \
           np.less_equal(point, upper).all()
    #return np.all(point >= lower) and np.all(upper >= point)

def contains_circle(point, circle, buffer=0.):
    center, radius = circle
    return np.linalg.norm(np.array(point) - np.array(center)) <= (radius + buffer)

def contains(point, shape, **kwargs):
    if isinstance(shape, Box):
        return contains_box(point, shape, **kwargs)
    if isinstance(shape, Circle):
        return contains_circle(point, shape, **kwargs)
    raise NotImplementedError(shape)

def point_collides(point, obstacles, buffer=MIN_PROXIMITY, **kwargs):
    return any(contains(point, obst, buffer=buffer, **kwargs) for obst in obstacles)

def sample_line(segment, step_size=STEP_SIZE):
    (q1, q2) = segment
    diff = get_delta(q1, q2)
    dist = np.linalg.norm(diff)
    for l in np.arange(0., dist, step_size):
        yield tuple(np.array(q1) + l * diff / dist)
    yield q2

def line_collides(line, obst, step_size=STEP_SIZE, **kwargs):  # TODO - could also compute this exactly
    return any(point_collides(point, obstacles=[obst], **kwargs) for point in sample_line(line, step_size=step_size))

def is_collision_free(line, obstacles, **kwargs):
    return not any(line_collides(line, obst, **kwargs) for obst in obstacles)

def create_box(center, extents):
    (x, y) = center
    (w, h) = extents
    lower = (x - w / 2., y - h / 2.)
    upper = (x + w / 2., y + h / 2.)
    return Box(np.array(lower), np.array(upper))

def create_cylinder(center, radius):
    return Circle(center, radius)

get_box_center = get_interval_center
get_box_extent = get_interval_extent

def sample_box(box):
    (lower, upper) = box
    return np.random.random(len(lower)) * get_box_extent(box) + lower

def sample_circle(circle):
    center, radius = circle
    theta = np.random.uniform(low=0, high=2*np.pi)
    direction = np.array([np.cos(theta), np.sin(theta)])
    return np.array(center) + direction

#################################################################

def draw_shape(viewer, shape, **kwargs):
    if viewer is None:
        return
    if isinstance(shape, Box):
        return viewer.draw_rectangle(shape, **kwargs)
    if isinstance(shape, Circle):
        center, radius = shape
        return viewer.draw_circle(center, radius, **kwargs)
    raise NotImplementedError(shape)

def draw_environment(obstacles, regions, **kwargs):
    try:
        viewer = PRMViewer(**kwargs)
    except TclError:
        traceback.print_exc()
        return None

    for obstacle in obstacles:
        draw_shape(viewer, obstacle, color='brown')
    for name, region in regions.items():
        if name == 'env':
            continue
        draw_shape(viewer, region, color='green')
    return viewer

def add_segments(viewer, segments, step_size=INF, **kwargs):
    if (viewer is None) or (segments is None):
        return
    for line in segments:
        viewer.draw_line(line, **kwargs)
        #for p in [p1, p2]:
        for p in sample_line(line, step_size=step_size):
            viewer.draw_point(p, radius=2, **kwargs)

def add_path(viewer, path, **kwargs):
    segments = list(get_pairs(path))
    return add_segments(viewer, segments, **kwargs)

def hex_from_8bit(rgb):
    assert all(0 <= v <= 2**8-1 for v in rgb)
    return '#%02x%02x%02x' % tuple(rgb)

def hex_from_rgb(rgb):
    assert all(0. <= v <= 1.for v in rgb)
    return hex_from_8bit([int(v*(2**8-1)) for v in rgb])

def spaced_colors(n, s=1, v=1):
    import colorsys
    return [colorsys.hsv_to_rgb(h, s, v) for h in np.linspace(0, 1, n, endpoint=False)]

def add_timed_path(viewer, times, path, **kwargs):
    if viewer is None:
        return
    # TODO: color based on velocity
    import colorsys

    min_value = min(times)
    max_value = max(times)

    def get_color(t):
        fraction = (t - min_value) / (max_value - min_value)
        rgb = colorsys.hsv_to_rgb(h=(1-fraction)*(2./3), s=1., v=1.)
        return hex_from_rgb(rgb)

    for t, p in zip(times, path):
        viewer.draw_point(p, radius=2, color=get_color(t), **kwargs)
    for (t1, p1), (t2, p2) in get_pairs(list(zip(times, path))):
        t = (t1 + t2) / 2
        line = (p1, p2)
        viewer.draw_line(line, color=get_color(t), **kwargs)

def draw_solution(segments, obstacles, regions):
    viewer = draw_environment(obstacles, regions)
    add_segments(viewer, segments)

def add_roadmap(viewer, roadmap, **kwargs):
    if viewer is None:
        return
    for line in roadmap:
        viewer.draw_line(line, **kwargs)

def draw_roadmap(roadmap, obstacles, regions):
    viewer = draw_environment(obstacles, regions)
    add_roadmap(viewer, roadmap)

def add_points(viewer, points, **kwargs):
    if viewer is None:
        return
    for sample in points:
        viewer.draw_point(sample, **kwargs)
