import math

import numpy as np

from .limits import maximize_curve
from .retime import EPSILON, curve_from_controls, check_time, parabolic_val, append_polys
from ..utils import INF, get_sign, waypoints_from_path, get_pairs

T_MIN = 1e-3 # EPSILON

def quickest_inf_accel(x1, x2, v_max=INF):
    #return solve_zero_ramp(x1, x2, v_max=INF)
    return abs(x2 - x1) / abs(v_max)

def acceleration_cost(p_curve, **kwargs):
    # TODO: minimize max acceleration
    # TODO: minimize sum of squared
    # TODO: minimize absolute value
    # total_accel = p_curve.integrate(p_curve.x[0], p_curve.x[-1]) # TODO: square or abs
    a_curve = p_curve.derivative(nu=2)
    max_t, max_a = maximize_curve(a_curve, **kwargs)
    return max_a

def find_lower_bound(x1, x2, v1=None, v2=None, v_max=None, a_max=None, ord=INF):
    d = len(x1)
    if v_max is None:
        v_max = np.full(d, INF)
    if a_max is None:
        a_max = np.full(d, INF)
    lower_bounds = [
        # Instantaneously accelerate
        np.linalg.norm(np.divide(np.subtract(x2, x1), v_max), ord=ord), # quickest_inf_accel
    ]
    if (v1 is not None) and (v2 is not None):
        lower_bounds.extend([
            np.linalg.norm(np.divide(np.subtract(v2, v1), a_max), ord=ord),
        ])
    return max(lower_bounds)

def check_curve(p_curve, x1, x2, v1, v2, T, v_max=INF, a_max=INF):
    assert p_curve is not None
    end_times = np.append(p_curve.x[:1], p_curve.x[-1:])
    v_curve = p_curve.derivative()
    # print()
    #print(x1, x2, v1, v2, T, v_max, a_max)
    # print([x1, x2], [float(p_curve(t)) for t in end_times])
    # print([v1, v2], [float(v_curve(t)) for t in end_times])
    if not np.allclose([0., T], end_times):
        raise RuntimeError([0., T], end_times)
    if not np.allclose([x1, x2], [float(p_curve(t)) for t in end_times]):
        raise RuntimeError([x1, x2], [float(p_curve(t)) for t in end_times])
    if not np.allclose([v1, v2], [float(v_curve(t)) for t in end_times]):
        raise RuntimeError([v1, v2], [float(v_curve(t)) for t in end_times])
    all_times = p_curve.x
    if not all(abs(v_curve(t)) <= abs(v_max) + EPSILON for t in all_times):
        raise RuntimeError(abs(v_max), [abs(v_curve(t)) for t in all_times])
    a_curve = v_curve.derivative()
    if not all(abs(a_curve(t)) <= abs(a_max) + EPSILON for t in all_times):
        raise RuntimeError(abs(a_max), [abs(a_curve(t)) for t in all_times])
    # TODO: check continuity

##################################################

def zero_one_fixed(x1, x2, T, v_max=INF, dt=EPSILON):
    from scipy.interpolate import PPoly
    #assert 0 < v_max < INF
    sign = get_sign(x2 - x1)
    d = abs(x2 - x1)
    v = d / T
    if v > (v_max + EPSILON):
        return None
    t_hold = T - 2*dt
    assert t_hold > 0
    v = d / t_hold
    #return zero_three_stage(x1, x2, T, v_max=v_max, a_max=1e6) # NOTE: approximation
    coeffs = [[0., x1], [sign*v, x1], [0., x2]]
    times = [0., dt, dt + t_hold, T]
    p_curve = PPoly(c=np.array(coeffs).T, x=times) # Not differentiable
    return p_curve


def zero_two_ramp(x1, x2, T, v_max=INF, a_max=INF):
    sign = get_sign(x2 - x1)
    d = abs(x2 - x1)
    t_accel = T / 2.
    a = d / t_accel ** 2 # Lower accel
    if a > (a_max + EPSILON):
        return None
    if a*t_accel > (v_max + EPSILON):
        return None
    a = min(a, a_max) # Numerical error
    durations = [t_accel, t_accel]
    accels = [sign * a, -sign * a]
    p_curve = curve_from_controls(durations, accels, x0=x1)
    return p_curve


def zero_three_stage(x1, x2, T, v_max=INF, a_max=INF):
    sign = get_sign(x2 - x1)
    d = abs(x2 - x1)
    solutions = np.roots([
        a_max,
        -a_max*T,
        d,
    ])
    solutions = filter(check_time, solutions)
    solutions = [t for t in solutions if (T - 2*t) >= 0]
    if not solutions:
        return None
    t1 = min(solutions)
    if t1*a_max > (v_max + EPSILON):
        return None
    #t1 = min(t1, v_max / a_max)
    t3 = t1
    t2 = T - t1 - t3 # Lower velocity
    durations = [t1, t2, t3]
    accels = [sign * a_max, 0., -sign * a_max]
    p_curve = curve_from_controls(durations, accels, x0=x1)
    return p_curve

##################################################

def opt_straight_line(x1, x2, v_max=INF, a_max=INF, t_min=T_MIN, only_duration=False):
    # TODO: solve for a given T which is higher than the min T
    # TODO: solve for all joints at once using a linear interpolator
    # Can always rest at the start of the trajectory if need be
    # Exploits symmetry
    assert (v_max > 0.) and (a_max > 0.)
    #assert (v_max < INF) or (a_max < INF) or (t_min > 0)
    #v_max = abs(x2 - x1) / abs(v_max)
    d = abs(x2 - x1)
    #assert d > 0
    # if v_max == INF:
    #     raise NotImplementedError()
    # TODO: more efficient version

    if a_max == INF:
        T = d / v_max
        T += 2 * EPSILON
        T = max(t_min, T)
        if only_duration:
            return T
        p_curve = zero_one_fixed(x1, x2, T, v_max=v_max)
        check_curve(p_curve, x1, x2, v1=0., v2=0., T=T, v_max=v_max, a_max=a_max)
        return p_curve

    t_accel = math.sqrt(d / a_max) # 1/2.*a*t**2 = d/2.
    if a_max*t_accel <= v_max:
        T = 2.*t_accel
        #a = a_max
        assert t_min <= T
        T = max(t_min, T)
        if only_duration:
            return T
        p_curve = zero_two_ramp(x1, x2, T, v_max, a_max)
        check_curve(p_curve, x1, x2, v1=0., v2=0., T=T, v_max=v_max, a_max=a_max)
        return p_curve

    t1 = t3 = (v_max - 0.) / a_max
    t2 = (d - 2 * parabolic_val(t1, a=a_max)) / v_max
    T = t1 + t2 + t3

    assert t_min <= T
    T = max(t_min, T)
    if only_duration:
        return T
    p_curve = zero_three_stage(x1, x2, T, v_max=v_max, a_max=a_max)
    check_curve(p_curve, x1, x2, v1=0., v2=0., T=T, v_max=v_max, a_max=a_max)
    return p_curve

##################################################

def get_default_limits(d=None, v_max=None, a_max=None):
    # TODO: unify with find_lower_bound
    # TODO: figure out d from the bounds
    if v_max is None:
        v_max = INF*np.ones(d)
    if a_max is None:
        a_max = INF*np.ones(d)
    assert len(v_max) == len(a_max)
    return v_max, a_max

def solve_linear(difference, v_max, a_max, **kwargs):
    # TODO: careful with circular joints
    # TODO: careful if difference is zero
    unit_v_max = min(np.divide(v_max, np.absolute(difference)))
    unit_a_max = min(np.divide(a_max, np.absolute(difference)))
    return opt_straight_line(x1=0., x2=1., v_max=unit_v_max, a_max=unit_a_max, **kwargs)

def solve_multi_linear(positions, v_max=None, a_max=None, **kwargs):
    from scipy.interpolate import PPoly
    assert positions
    positions = waypoints_from_path(positions, **kwargs)
    if len(positions) == 1:
        positions.append(positions[-1])
    #assert len(positions) >= 2
    d = len(positions[0])
    v_max, a_max = get_default_limits(d, v_max=v_max, a_max=a_max)
    splines = []
    for x1, x2 in get_pairs(positions):
        difference = np.subtract(x2, x1)  # TODO: pass if too small
        curve = solve_linear(difference, v_max, a_max)
        c = np.zeros(curve.c.shape + (d,))
        for k in range(d):
            c[:,:,k] = difference[k]*curve.c
            c[-1,:,k] += x1[k]
        splines.append(PPoly(c=c, x=curve.x))
    assert splines
    return append_polys(*splines)
