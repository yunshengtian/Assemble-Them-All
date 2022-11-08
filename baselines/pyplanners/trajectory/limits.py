import time
import numpy as np

from .retime import EPSILON, get_max_velocity, poly_from_spline, get_interval, filter_extrema, find_extrema
from ..utils import INF, elapsed_time

def old_check_spline(spline, v_max=None, a_max=None, start_idx=None, end_idx=None):
    # TODO: be careful about time vs index (here is index)
    if (v_max is None) and (a_max is None):
        return True
    if start_idx is None:
        start_idx = 0
    if end_idx is None:
        end_idx = len(spline.x) - 1
    signs = [+1, -1]
    for i in range(start_idx, end_idx):
        t0, t1 = spline.x[i], spline.x[i + 1]
        t0, t1 = 0, (t1 - t0)
        boundary_ts = [t0, t1]
        for k in range(spline.c.shape[-1]):
            position_poly = poly_from_spline(spline, i, d=k)
            # print([position_poly(t) for t in boundary_ts])
            # print([spline(t)[k] for t in boundary_ts])

            if v_max is not None:
                vel_poly = position_poly.deriv(m=1)
                if any(abs(vel_poly(t)) > v_max[k] for t in boundary_ts):
                    return False
                if any(not isinstance(r, complex) and (t0 <= r <= t1)
                       for s in signs for r in (vel_poly + s * np.poly1d([v_max[k]])).roots):
                    return False
            # a_max = None

            # TODO: reorder to check endpoints first
            if a_max is not None:  # INF
                accel_poly = position_poly.deriv(m=2)
                # print([(accel_poly(t), a_max[k]) for t in boundary_ts])
                if any(abs(accel_poly(t)) > a_max[k] for t in boundary_ts):
                    return False
                # print([accel_poly(r) for s in signs for r in (accel_poly + s*np.poly1d([a_max[k]])).roots
                #        if isinstance(r, complex) and (t0 <= r <= t1)])
                if any(not isinstance(r, complex) and (t0 <= r <= t1)
                       for s in signs for r in (accel_poly + s * np.poly1d([a_max[k]])).roots):
                    return False
    return True

def check_spline(spline, v_max=None, a_max=None, verbose=False, **kwargs):
    if (v_max is None) and (a_max is None):
        return True
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PPoly.html#scipy.interpolate.PPoly
    # polys = [np.poly1d([spline.c[c, 0, k] for c in range(spline.c.shape[0])]) # Decreasing order
    #         for k in range(spline.c.shape[-1])]
    # from numpy.polynomial.polynomial import Polynomial
    # polys = [Polynomial([spline.c[c, 0, k] for c in range(spline.c.shape[0])][-1])  # Increasing order
    #          for k in range(spline.c.shape[-1])]

    # _, v = find_max_velocity(spline, ord=INF)
    # _, a = find_max_acceleration(spline, ord=INF)
    # print(v, v_max, a, a_max)
    # input()

    if v_max is not None:
        # TODO: maybe the pieces are screwing something
        t, v = find_max_velocity(spline, **kwargs)
        violation = abs(v) > get_max_velocity(v_max) + EPSILON
        if verbose: # or violation:
            print('Violation: {} | Max velocity: {:.3f}/{:.3f} (at time {:.3f})'.format(
                violation, v, get_max_velocity(v_max), t))
        if violation:
            return False
    return True
    # TODO: trusting continuous pieces to respect
    # TODO: solve for the intersection when not at a root
    # if a_max is not None:
    #     t, a = find_max_acceleration(spline, **kwargs)
    #     print('Max accel: {:.3f}/{:.3f} (at time {:.3f})'.format(a, get_max_velocity(a_max), t))
    #     if abs(a) > get_max_velocity(a_max) + EPSILON:
    #         print(t, a)
    #         print(spline.x)
    #         print(t in spline.x)
    #         input()
    #         return False
    # return True

##################################################

def minimize_objective(objective, lower, upper, num=100, max_time=INF, max_iterations=100, verbose=False, **kwargs):
    # https://www.cvxpy.org/examples/basic/socp.html
    from scipy.optimize import minimize #, line_search, brute, basinhopping, minimize_scalar
    start_time = time.time()
    best_t, best_f = None, INF
    bounds = list(zip(lower, upper))
    assert num >= 1
    for iteration in range(num):
        if elapsed_time(start_time) >= max_time:
            break
        x0 = np.random.uniform(lower, upper)
        if max_iterations is None:
            t, f = x0, objective(x0)
        else:
            result = minimize(objective, x0=x0, bounds=bounds, # maxiter=max_iterations,
                              **kwargs) # method=None, jac=None,
            t, f = result.x, result.fun
        if (best_t is None) or (f < best_f):
            best_t, best_f = t, f
            if verbose:
                print(iteration, x0, t, f) # objective(result.x)
    return best_t, best_f

def find_max_curve(curve, start_t=None, end_t=None, norm=INF, **kwargs):
    start_t, end_t = get_interval(curve, start_t=start_t, end_t=end_t)
    # TODO: curve(t) returns a matrix if passed a vector, which is summed by the INF norm
    objective = lambda t: -np.linalg.norm(curve(t[0]), ord=norm) # 2 | INF
    #objective = lambda t: -np.linalg.norm(curve(t[0]), norm=2)**2 # t[0]
    #accelerations_curve = positions_curve.derivative() # TODO: ValueError: failed in converting 7th argument `g' of _lbfgsb.setulb to C/Fortran array
    #grad = lambda t: np.array([-2*sum(accelerations_curve(t))])
    #result = minimize_scalar(objective, method='bounded', bounds=(start_t, end_t)) #, options={'disp': False})

    #print(max(-objective(t) for t in curve.x))
    best_t, best_f = minimize_objective(objective, lower=[start_t], upper=[end_t], **kwargs)
    #best_f = objective(best_t)
    best_t, best_f = best_t[0], -best_f
    return best_t, best_f

##################################################

def find_candidates(curve, **kwargs):
    return filter_extrema(curve, list(curve.x) + find_extrema(curve, **kwargs))

def maximize_curve(curve, ignore=set(), **kwargs): # fn=None
    #d = curve(start_t)
    times = find_candidates(curve, **kwargs)
    times = [t for t in times if t not in ignore] # TODO: filter repeated
    #fn = lambda t: max(np.absolute(curve(t)))
    fn = lambda t: np.linalg.norm(curve(t), ord=INF) if curve(t).shape else float(curve(t))
    #fn = max
    max_t = max(times, key=fn)
    return max_t, fn(max_t)

def exceeds_curve(curve, threshold, start_t=None, end_t=None, **kwargs):
    # TODO: joint limits
    # TODO: solve for the intersection with threshold
    max_t, max_v = maximize_curve(curve, start_t=start_t, end_t=end_t, **kwargs)
    if np.greater(max_v, threshold):
        return max_t
    return None

##################################################

def find_max_velocity(positions_curve, analytical=True, **kwargs):
    velocities_curve = positions_curve.derivative(nu=1)
    if analytical:
        return maximize_curve(velocities_curve)
    else:
        return find_max_curve(velocities_curve, **kwargs)

def find_max_acceleration(positions_curve, **kwargs):
    # TODO: should only ever be quadratic
    accelerations_curve = positions_curve.derivative(nu=2)
    #return find_max_curve(accelerations_curve, max_iterations=None, **kwargs)
    return maximize_curve(accelerations_curve, discontinuity=True,)
                          #ignore=set(positions_curve.x))

def analyze_continuity(curve, epsilon=1e-9, **kwargs):
    # TODO: explicitly check the adjacent curves
    start_t, end_t = get_interval(curve, **kwargs)
    max_t, max_error = start_t, 0.
    for i in range(1, len(curve.x)-1):
        t = curve.x[i]
        if not start_t <= t <= end_t:
            continue
        t1 = t - epsilon # TODO: check i-1
        t2 = t + epsilon # TODO: check i+1
        v1 = curve(t1)
        v2 = curve(t2)
        error = np.linalg.norm(v2 - v1, ord=INF)
        if error > max_error:
            max_t, max_error = t, error
    return max_t, max_error
