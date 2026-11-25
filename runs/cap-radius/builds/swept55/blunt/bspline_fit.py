"""Small tensor-product B-spline least-squares fitter (NumPy) used by build_cap.py."""
import numpy as np


def clamped_knots(n_poles, degree, breaks=None):
    """Clamped knot vector on [0, 1]; interior knots at 'breaks' (default uniform)."""
    n_int = n_poles - degree - 1
    if breaks is None:
        breaks = np.linspace(0, 1, n_int + 2)[1:-1]
    assert len(breaks) == n_int
    return np.r_[np.zeros(degree + 1), breaks, np.ones(degree + 1)]


def basis_matrix(t, knots, degree):
    """B-spline basis values N_i(t) as an (len(t), n_poles) matrix (Cox - de Boor)."""
    t = np.asarray(t, float)
    n = len(knots) - degree - 1
    m = len(knots) - 1
    N = np.zeros((len(t), m))
    for i in range(m):
        lo, hi = knots[i], knots[i + 1]
        if hi > lo:
            N[:, i] = (t >= lo) & (t < hi)
    # t == 1 belongs to the last non-empty span
    last = max(i for i in range(m) if knots[i + 1] > knots[i])
    N[t >= knots[-1], :] = 0
    N[t >= knots[-1], last] = 1
    for d in range(1, degree + 1):
        Nn = np.zeros((len(t), m - d))
        for i in range(m - d):
            a = knots[i + d] - knots[i]
            b = knots[i + d + 1] - knots[i + 1]
            term = 0
            if a > 0:
                term = term + (t - knots[i]) / a * N[:, i]
            if b > 0:
                term = term + (knots[i + d + 1] - t) / b * N[:, i + 1]
            Nn[:, i] = term
        N = Nn
    return N[:, :n]


def constrained_lstsq(A, Y, C=None, d=None):
    """min |A x - Y| subject to C x = d (column-wise for Y with several columns)."""
    if C is None:
        return np.linalg.lstsq(A, Y, rcond=None)[0]
    n = A.shape[1]
    k = C.shape[0]
    K = np.block([[A.T @ A, C.T], [C, np.zeros((k, k))]])
    rhs = np.vstack([A.T @ Y, d if d is not None else np.zeros((k, Y.shape[1]))])
    sol = np.linalg.solve(K, rhs)
    return sol[:n]
