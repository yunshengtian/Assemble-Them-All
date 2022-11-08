import numpy as np


def test_spline(best_t, x1, x2, v1, v2):
    observations = [
        (0., x1[0], 0),
        (best_t, x2[0], 0),
        (0., v1[0], 1),
        (best_t, v2[0], 1),
    ]
    degree = len(observations) - 1

    from numpy import poly1d
    terms = []
    for k in range(degree + 1):
        coeffs = np.zeros(degree + 1)
        coeffs[k] = 1.
        terms.append(poly1d(coeffs))
    # series = poly1d(np.ones(degree+1))

    A = []
    b = []
    for t, v, nu in observations:
        A.append([term.deriv(m=nu)(t) for term in terms])
        b.append(v)
    print(A)
    print(b)
    print(np.linalg.solve(A, b))
    # print(polyfit([t for t, _, nu in observations if nu == 0],
    #              [v for _, v, nu in observations if nu == 0], deg=degree))
    # TODO: compare with CubicHermiteSpline