import random
import time

import numpy as np

from .linear import find_lower_bound, solve_linear, T_MIN
from .limits import check_spline
from .discretize import time_discretize_curve, derivative_discretize_curve, distance_discretize_curve, sample_discretize_curve
from .parabolic import solve_multi_poly, solve_multivariate_ramp
from .retime import EPSILON, trim, spline_duration, append_polys, get_interval, find_extrema
from ..utils import INF, elapsed_time, get_pairs, find, default_selector, waypoints_from_path, irange

def within_velocity_limits(position_curve, max_v=None, **kwargs):
    if max_v is None:
        return True
    velocity_curve = position_curve.derivative(nu=1)
    extrema = find_extrema(velocity_curve, **kwargs)
    return all(np.less_equal(np.absolute(velocity_curve(t)), max_v).all() for t in extrema)


def within_acceleration_limits(position_curve, max_a=None, **kwargs):
    velocity_curve = position_curve.derivative(nu=1)
    return within_velocity_limits(velocity_curve, max_v=max_a, **kwargs)


def within_dynamical_limits(position_curve, max_v=None, max_a=None, **kwargs):
    return within_velocity_limits(position_curve, max_v=max_v, **kwargs) and \
           within_acceleration_limits(position_curve, max_v=max_a, **kwargs)

##################################################

def get_curve_collision_fn(collision_fn=lambda q: False, max_velocities=None, max_accelerations=None): # a_max

    def curve_collision_fn(curve, t0=None, t1=None):
        # TODO: stage the function to check all the easy things like joint limits first
        if curve is None:
            return True
        #if not within_dynamical_limits(curve, max_v=max_velocities, max_a=max_accelerations, start_t=t0, end_t=t1):
        #    return True
        # TODO: can exactly compute limit violations
        # if not check_spline(curve, v_max=max_velocities, a_max=None, verbose=False,
        #                     #start_t=t0, end_t=t1,
        #                     ):
        #     return True
        # _, samples = time_discretize_curve(curve, verbose=False, start_t=t0, end_t=t1, #max_velocities=v_max)
        _, samples = distance_discretize_curve(curve, start_t=t0, end_t=t1)
        if any(map(collision_fn, default_selector(samples))):
           return True
        return False
    return curve_collision_fn

##################################################

def smooth_curve(start_curve, v_max, a_max, curve_collision_fn,
                 sample=True, intermediate=True, cubic=True, refit=True, num=1000, min_improve=0., max_time=INF):
    # TODO: rename smoothing.py to shortcutting.py
    # TODO: default v_max and a_max
    assert (num < INF) or (max_time < INF)
    assert refit or intermediate
    # TODO: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html
    #from scipy.interpolate import UnivariateSpline, LSQUnivariateSpline, LSQBivariateSpline
    from scipy.interpolate import CubicHermiteSpline
    start_time = time.time()

    if curve_collision_fn(start_curve, t0=None, t1=None):
        #return None
        return start_curve
    curve = start_curve
    for iteration in irange(num):
        if elapsed_time(start_time) >= max_time:
            break
        times = curve.x
        durations = [0.] + [t2 - t1 for t1, t2 in get_pairs(times)] # includes start
        positions = [curve(t) for t in times]
        velocities = [curve(t, nu=1) for t in times]

        # ts = [times[0], times[-1]]
        # t1, t2 = curve.x[0], curve.x[-1]
        t1, t2 = np.random.uniform(times[0], times[-1], 2) # TODO: sample based on position
        if t1 > t2: # TODO: minimum distance from a knot
            t1, t2 = t2, t1

        ts = [t1, t2]
        i1 = find(lambda i: times[i] <= t1, reversed(range(len(times)))) # index before t1
        i2 = find(lambda i: times[i] >= t2, range(len(times))) # index after t2
        assert i1 != i2

        local_positions = [curve(t) for t in ts]
        local_velocities = [curve(t, nu=1) for t in ts]
        #print(local_velocities, v_max)
        assert all(np.less_equal(np.absolute(v), v_max + EPSILON).all() for v in local_velocities)
        #if any(np.greater(np.absolute(v), v_max).any() for v in local_velocities):
        #    continue # TODO: do the same with collisions
        x1, x2 = local_positions
        v1, v2 = local_velocities

        #min_t = 0
        min_t = find_lower_bound(x1, x2, v1, v2, v_max=v_max, a_max=a_max)
        #min_t = optimistic_time(x1, x2, v_max=v_max, a_max=a_max)
        current_t = (t2 - t1) - min_improve
        if min_t >= current_t: # TODO: also limit the distance/duration between these two points
            continue

        #best_t = min_t
        if sample:
            max_t = current_t
            ramp_t = solve_multivariate_ramp(x1, x2, v1, v2, v_max, a_max)
            ramp_t = INF if ramp_t is None else ramp_t
            max_t = min(max_t, ramp_t)
            best_t = random.uniform(min_t, max_t)
        else:
            best_t = solve_multivariate_ramp(x1, x2, v1, v2, v_max, a_max)
        if (best_t is None) or (best_t >= current_t):
            continue
        #best_t += 1e-3
        #print(min_t, best_t, current_t)
        local_durations = [t1 - times[i1], best_t, times[i2] - t2]
        #local_times = [0, best_t]
        local_times = [t1, (t1 + best_t)] # Good if the collision function is time sensitive

        if intermediate:
            if cubic:
                local_curve = CubicHermiteSpline(local_times, local_positions, dydx=local_velocities)
            else:
                local_curve = solve_multi_poly(times=local_times, positions=local_positions,
                                               velocities=local_velocities,
                                               v_max=v_max, a_max=a_max)
            if (local_curve is None) or (spline_duration(local_curve) >= current_t) \
                    or curve_collision_fn(local_curve, t0=None, t1=None):
                continue
            # print(new_curve.hermite_spline().c[0,...])
            local_positions = [local_curve(x) for x in local_curve.x]
            local_velocities = [local_curve(x, nu=1) for x in local_curve.x]
            local_durations = [t1 - times[i1]] + [x - local_curve.x[0]
                                                    for x in local_curve.x[1:]] + [times[i2] - t2]

        if refit:
            new_durations = np.concatenate([
                durations[:i1 + 1], local_durations, durations[i2 + 1:]])
            # assert len(new_durations) == (i1 + 1) + (len(durations) - i2) + 2
            new_times = np.cumsum(new_durations)
            # new_times = [new_times[0]] + [t2 for t1, t2 in get_pairs(new_times) if t2 > t1]
            new_positions = positions[:i1 + 1] + local_positions + positions[i2:]
            new_velocities = velocities[:i1 + 1] + local_velocities + velocities[i2:]
            # if not all(np.less_equal(np.absolute(v), v_max).all() for v in new_velocities):
            #    continue
            if cubic:
                # new_curve = CubicSpline(new_times, new_positions)
                new_curve = CubicHermiteSpline(new_times, new_positions, dydx=new_velocities)
            else:
                new_curve = solve_multi_poly(new_times, new_positions, new_velocities, v_max, a_max)
            if (new_curve is None) or (spline_duration(new_curve) >= spline_duration(curve)) \
                    or not check_spline(new_curve, v_max, a_max) or \
                    (not intermediate and curve_collision_fn(new_curve, t0=None, t1=None)):
                continue
        else:
            assert intermediate
            # print(curve.x)
            # print(curve.c[...,0])
            # pre_curve = trim(curve, end=t1)
            # post_curve = trim(curve, start=t1)
            # curve = append_polys(pre_curve, post_curve)
            # print(curve.x)
            # print(curve.c[...,0])

            # print(new_curve.x)
            # print(new_curve.c[...,0])
            pre_curve = trim(curve, end=t1)
            post_curve = trim(curve, start=t2)
            new_curve = append_polys(pre_curve, local_curve, post_curve) # TODO: the numerics are throwing this off?
            # print(new_curve.x)
            # print(new_curve.c[...,0])
            #assert(not curve_collision_fn(new_curve, t0=None, t1=None))
            if (spline_duration(new_curve) >= spline_duration(curve)) or \
                    not check_spline(new_curve, v_max, a_max):
                continue
        print('Iterations: {} | Current time: {:.3f} | New time: {:.3f} | Elapsed time: {:.3f}'.format(
            iteration, spline_duration(curve), spline_duration(new_curve), elapsed_time(start_time)))
        curve = new_curve
    print('Iterations: {} | Start time: {:.3f} | End time: {:.3f} | Elapsed time: {:.3f}'.format(
        num, spline_duration(start_curve), spline_duration(curve), elapsed_time(start_time)))
    check_spline(curve, v_max, a_max)
    return curve

##################################################

def smooth_cubic(path, collision_fn, resolutions, v_max=None, a_max=None, time_step=1e-2,
                 parabolic=True, sample=False, intermediate=True, max_iterations=1000, max_time=INF,
                 min_improve=0., verbose=False):
    start_time = time.time()
    if path is None:
        return None
    assert (v_max is not None) or (a_max is not None)
    assert path and (max_iterations < INF) or (max_time < INF)
    from scipy.interpolate import CubicHermiteSpline

    def curve_collision_fn(segment, t0=None, t1=None):
        #if not within_dynamical_limits(curve, max_v=v_max, max_a=a_max, start_t=t0, end_t=t1):
        #    return True
        _, samples = sample_discretize_curve(segment, resolutions, start_t=t0, end_t=t1, time_step=time_step)
        if any(map(collision_fn, default_selector(samples))):
           return True
        return False

    start_positions = waypoints_from_path(path) # TODO: ensure following the same path (keep intermediate if need be)
    if len(start_positions) == 1:
        start_positions.append(start_positions[-1])

    start_durations = [0] + [solve_linear(np.subtract(p2, p1), v_max, a_max, t_min=T_MIN, only_duration=True)
                             for p1, p2 in get_pairs(start_positions)] # TODO: does not assume continuous acceleration
    start_times = np.cumsum(start_durations) # TODO: dilate times
    start_velocities = [np.zeros(len(start_positions[0])) for _ in range(len(start_positions))]
    start_curve = CubicHermiteSpline(start_times, start_positions, dydx=start_velocities)
    # TODO: directly optimize for shortest spline
    if len(start_positions) <= 2:
        return start_curve

    curve = start_curve
    for iteration in irange(max_iterations):
        if elapsed_time(start_time) >= max_time:
            break
        times = curve.x
        durations = [0.] + [t2 - t1 for t1, t2 in get_pairs(times)]
        positions = [curve(t) for t in times]
        velocities = [curve(t, nu=1) for t in times]

        t1, t2 = np.random.uniform(times[0], times[-1], 2)
        if t1 > t2:
            t1, t2 = t2, t1
        ts = [t1, t2]
        i1 = find(lambda i: times[i] <= t1, reversed(range(len(times)))) # index before t1
        i2 = find(lambda i: times[i] >= t2, range(len(times))) # index after t2
        assert i1 != i2

        local_positions = [curve(t) for t in ts]
        local_velocities = [curve(t, nu=1) for t in ts]
        if not all(np.less_equal(np.absolute(v), np.array(v_max) + EPSILON).all() for v in local_velocities):
            continue

        x1, x2 = local_positions
        v1, v2 = local_velocities

        current_t = (t2 - t1) - min_improve # TODO: percent improve
        #min_t = 0
        min_t = find_lower_bound(x1, x2, v1, v2, v_max=v_max, a_max=a_max)
        if parabolic:
            # Softly applies limits
            min_t = solve_multivariate_ramp(x1, x2, v1, v2, v_max, a_max) # TODO: might not be feasible (soft constraint)
            if min_t is None:
                continue
        if min_t >= current_t:
            continue
        best_t = random.uniform(min_t, current_t) if sample else min_t

        local_durations = [t1 - times[i1], best_t, times[i2] - t2]
        #local_times = [0, best_t]
        local_times = [t1, (t1 + best_t)] # Good if the collision function is time varying

        if intermediate:
            local_curve = CubicHermiteSpline(local_times, local_positions, dydx=local_velocities)
            if curve_collision_fn(local_curve, t0=None, t1=None): # check_spline
                continue
            #local_positions = [local_curve(x) for x in local_curve.x]
            #local_velocities = [local_curve(x, nu=1) for x in local_curve.x]
            local_durations = [t1 - times[i1]] + [x - local_curve.x[0] for x in local_curve.x[1:]] + [times[i2] - t2]

        new_durations = np.concatenate([durations[:i1 + 1], local_durations, durations[i2 + 1:]])
        new_times = np.cumsum(new_durations)
        new_positions = positions[:i1 + 1] + local_positions + positions[i2:]
        new_velocities = velocities[:i1 + 1] + local_velocities + velocities[i2:]

        new_curve = CubicHermiteSpline(new_times, new_positions, dydx=new_velocities)
        if not intermediate and curve_collision_fn(new_curve, t0=None, t1=None):
            continue
        if verbose:
            print('Iterations: {} | Current time: {:.3f} | New time: {:.3f} | Elapsed time: {:.3f}'.format(
                iteration, spline_duration(curve), spline_duration(new_curve), elapsed_time(start_time)))
        curve = new_curve
    if verbose:
        print('Iterations: {} | Start time: {:.3f} | End time: {:.3f} | Elapsed time: {:.3f}'.format(
            max_iterations, spline_duration(start_curve), spline_duration(curve), elapsed_time(start_time)))
    return curve

##################################################

DERIVATIVE_NAMES = [
    'Position',
    'Velocity',
    'Acceleration',
    'Jerk',
]

def plot_curve(positions_curve, derivative=0, dt=1e-3):
    import matplotlib.pyplot as plt
    # test_scores_mean, test_scores_std = estimate_gaussian(test_scores)
    # width = scale * test_scores_std # standard deviation
    # # TODO: standard error (confidence interval)
    # # from learn_tools.active_learner import tail_confidence
    # # alpha = 0.95
    # # scale = tail_confidence(alpha)
    # # width = scale * test_scores_std / np.sqrt(train_sizes)
    # plt.fill_between(train_sizes, test_scores_mean - width, test_scores_mean + width, alpha=0.1)

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color'] # Default color order
    curve = positions_curve.derivative(nu=derivative)

    start_t, end_t = get_interval(positions_curve, start_t=None, end_t=None)
    times = np.append(np.arange(start_t, end_t, step=dt), [end_t])
    for i, coords in enumerate(zip(*[curve(t) for t in times])):
        plt.plot(times, coords, color=colors[i], label='x[{}]'.format(i)) #, marker='o-')

    if derivative == 0:
        #discretized_times = np.append(np.arange(start_t, end_t, step=5e1*dt), [end_t])
        resolution = 5e-2
        #discretized_times, _ = time_discretize_curve(curve, resolution=resolution)
        #discretized_times, _ = derivative_discretize_curve(curve, resolution=resolution)
        discretized_times, _ = distance_discretize_curve(curve, resolution=resolution)
        print('Discretize steps:', np.array(discretized_times))
        for i, coords in enumerate(zip(*[curve(t) for t in discretized_times])):
            plt.plot(discretized_times, coords, color=colors[i], label='x[{}]'.format(i), marker='x') # o | + | x

    for t in curve.x:
        plt.axvline(x=t, color='black')
    extrema = find_extrema(curve)
    print('Extrema:', np.array(extrema))
    #plt.vlines(extrema)
    for t in extrema:
        plt.axvline(x=t, color='green')

    plt.xlabel('Time')
    plt.ylabel(DERIVATIVE_NAMES[derivative])
    ax = plt.subplot()
    ax.autoscale(tight=True)
    plt.legend(loc='best') # 'upper left'
    plt.grid()
    plt.show()
    #return curve
