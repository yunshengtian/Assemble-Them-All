import numpy as np

from .retime import get_interval, spline_duration
from .limits import find_max_velocity, find_max_acceleration
from ..utils import get_distance, INF, even_space

V_MAX = 0.8*np.ones(2)
A_MAX = abs(V_MAX - 0.) / abs(0.2 - 0.)
#V_MAX = INF*np.ones(2)
#A_MAX = 1e6*np.ones(2)
DEFAULT_RESOLUTION = 1e-2

##################################################

def filter_proximity(times, positions, resolutions):
    assert len(times) == len(positions)
    if len(times) <= 2:
        return times, positions
    new_times = [times[0]]
    new_positions = [positions[0]]
    for idx in range(1, len(times)-1):
        # TODO: search a list of existing samples (such as knot points)
        current_delta = np.absolute(get_distance(positions[idx], new_positions[-1]) / resolutions)
        next_delta = np.absolute(get_distance(positions[idx+1], new_positions[-1]) / resolutions)
        if (current_delta >= 1).any() or (next_delta >= 1).any():
            new_times.append(times[idx])
            new_positions.append(positions[idx])
    new_times.append(times[-1])
    new_positions.append(positions[-1])
    return new_times, new_positions

##################################################

def inf_norm(vector):
    #return max(map(abs, vector))
    return np.linalg.norm(vector, ord=INF)

def time_discretize_curve(positions_curve, max_velocities=None,
                          resolution=DEFAULT_RESOLUTION, verbose=True, **kwargs): # TODO: min_time?
    start_t, end_t = get_interval(positions_curve, **kwargs)
    norm = INF
    d = len(positions_curve(start_t))
    resolutions = resolution
    if np.isscalar(resolution):
        resolutions = resolution*np.ones(d)
    if max_velocities is None:
        # TODO: adjust per trajectory segment
        v_max_t, max_v = find_max_velocity(positions_curve, start_t=start_t, end_t=end_t, norm=norm)
        a_max_t, max_a = find_max_acceleration(positions_curve, start_t=start_t, end_t=end_t, norm=norm)
        #v_max_t, max_v = INF, np.linalg.norm(V_MAX)
        time_step = min(np.divide(resolutions, max_v))
        #time_step = 0.1*time_step
        if verbose:
            print('Max velocity: {:.3f}/{:.3f} (at time {:.3f}) | Max accel: {:.3f}/{:.3f} (at time {:.3f}) | '
                  'Step: {:.3f} | Duration: {:.3f}'.format(
                max_v, np.linalg.norm(V_MAX, ord=norm), v_max_t, max_a, np.linalg.norm(A_MAX, ord=norm), a_max_t,
                time_step, spline_duration(positions_curve))) # 2 | INF
    else:
        time_step = np.min(np.divide(resolutions, max_velocities))

    times = even_space(start_t, end_t, step=time_step)
    positions = [positions_curve(t) for t in times]
    times, positions = filter_proximity(times, positions, resolution)
    return times, positions

    # TODO: bug here (just use knot points instead?)
    # times.extend(np.hstack(positions_curve.derivative().roots(discontinuity=True))) # TODO: make these points special within filter proximity
    # times = sorted(set(times))
    # positions = [positions_curve(t) for t in times]
    # return times, positions


def derivative_discretize_curve(positions_curve, resolution=DEFAULT_RESOLUTION, time_step=1e-3, **kwargs):
    d = positions_curve.c.shape[-1]
    resolutions = resolution*np.ones(d)
    start_t, end_t = get_interval(positions_curve, **kwargs)
    velocities_curve = positions_curve.derivative()
    #acceleration_curve = velocities_curve.derivative()
    times = [start_t]
    while True:
        velocities = velocities_curve(times[-1])
        dt = min(np.divide(resolutions, np.absolute(velocities)))
        #dt = min(dt, time_step) # TODO: issue if an infinite derivative
        new_time = times[-1] + dt
        if new_time > end_t:
            break
        times.append(new_time)
    times.append(end_t)
    positions = [positions_curve(control_time) for control_time in times]
    # TODO: distance between adjacent positions
    return times, positions

##################################################

def sample_discretize_curve(positions_curve, resolutions, time_step=1e-2, **kwargs):
    start_t, end_t = get_interval(positions_curve, **kwargs)
    times = [start_t]
    samples = [positions_curve(start_t)]
    for t in even_space(start_t, end_t, step=time_step):
        positions = positions_curve(t)
        if np.less_equal(samples[-1] - resolutions, positions).all() and  \
                np.less_equal(positions, samples[-1] + resolutions).all():
            continue
        times.append(t)
        samples.append(positions)
    return times, samples

def distance_discretize_curve(curve, resolution=DEFAULT_RESOLUTION, **kwargs):
    # TODO: could compute for the full interval and then sort by proximity for speed purposes
    # TODO: could sample at a small timestep and prune
    d = curve.c.shape[-1]
    resolutions = resolution*np.ones(d)
    start_t, end_t = get_interval(curve, **kwargs)
    times = [start_t]
    while True:
        t1 = times[-1]
        position = curve(t1)
        candidates = []
        for sign in [-1, +1]:
            target = position + sign*resolutions
            for i in range(len(position)):
                candidates.extend(t for t in curve.solve(target[i])[i] if t1 < t < end_t)
                # curve.roots()
        if not candidates:
            break
        times.append(min(candidates))
    times.append(end_t)
    positions = [curve(control_time) for control_time in times]
    # TODO: record distance from adjacent positions
    return times, positions


def integral_discretize_curve(positions_curve, resolution=DEFAULT_RESOLUTION, **kwargs):
    #from scipy.integrate import quad
    start_t, end_t = get_interval(positions_curve, **kwargs)
    distance_curve = positions_curve.antiderivative()
    #distance = positions_curve.integrate(a, b)
    # TODO: compute a total distance curve
    raise NotImplementedError()
